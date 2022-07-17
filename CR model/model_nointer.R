library(MASS); library(reshape2); library(ggplot2)
library(foreach); library(doParallel)
# library(SpatialFunctionsCR);
library(Rcpp); library(RcppArmadillo)

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
# stopCluster(cl)

setwd("/cluster/home/osr/ad625/Spatial CR")
sourceCpp("code.cpp")
source('functions.r', echo=TRUE)

# REAL DATA ---------------------------------------------------------------

usingTrueData <- T
load("/cluster/home/osr/ad625/Spatial CR/yad_data.rda")

traps$x <- traps$x / 1000
traps$y <- traps$y / 1000

# traps$x <- (traps$x - mean(traps$x)) / sd(traps$x)
# traps$y <- (traps$y - mean(traps$y)) / sd(traps$y)

trapsX <- traps$x
trapsY <- traps$y

trapsX_mean <- mean(trapsX)
trapsY_mean <- mean(trapsY)
trapsX_sd <- sd(trapsX)
trapsY_sd <- sd(trapsY)

traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2

trapsX <- (trapsX - mean(trapsX)) / traps_meansd
trapsY <- (trapsY - mean(trapsY)) / traps_meansd

allTraps <- cbind(trapsX, trapsY)
R_traps <- .2

a1 <- min(trapsX) - R_traps
b1 <- max(trapsX) + R_traps
a2 <- min(trapsY) - R_traps
b2 <- max(trapsY) + R_traps

# CREATE NEW TRAPS -------

a1 <- 0
b1 <- 1
a2 <- 0
b2 <- 1

usingTrueData <- F

numTraps <- 7
trapsX_grid <- seq(a1, b1, length.out = numTraps + 2)[-(numTraps+2)][-1]
trapsY_grid <- seq(a2, b2, length.out = numTraps + 2)[-(numTraps+2)][-1]
traps <- expand.grid(trapsX_grid, trapsY_grid)

ggplot(data = NULL, aes(x = traps[,1], y = traps[,2])) + geom_point()

trapsX <- traps[,1]
trapsY <- traps[,2]

trapsX_mean <- mean(trapsX)
trapsY_mean <- mean(trapsY)
trapsX_sd <- sd(trapsX)
trapsY_sd <- sd(trapsY)

K <- length(trapsX)

S_usage <- rep(1, K)

traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2

trapsX <- (trapsX - mean(trapsX)) / traps_meansd
trapsY <- (trapsY - mean(trapsY)) / traps_meansd

allTraps <- cbind(trapsX, trapsY)
R_traps <- .75

a1 <- min(trapsX) - R_traps
b1 <- max(trapsX) + R_traps
a2 <- min(trapsY) - R_traps
b2 <- max(trapsY) + R_traps


# SIMULATED DATA ----------------------------------------------------------

usingTrueData <- F

# starting params
{
  theta1 <- exp(-100)
  theta2 <- exp(-100)
  theta12 <- exp(-100)
  beta1 <- 250
  beta2 <- 250
  
  beta1_true <- beta1
  beta2_true <- beta2
  theta1_true <- theta1
  theta2_true <- theta2
  theta12_true <- theta12
  
  p0_1 <- .025
  p0_2 <- .025
  sigma_1 <- 3 / trapsX_sd
  sigma_2 <- 3 / trapsX_sd
  
  p0_1_true <- p0_1
  p0_2_true <- p0_2
  
  sigma_1_true <- sigma_1
  sigma_2_true <- sigma_2
  
}

# data simulation
{
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                        Sigma_newpoint = diag(1, nrow = 2), allTraps, .2,
                                        a1, b1, a2, b2)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  
  s1 <- s1[1:N1,]
  s2 <- s2[1:N2,]
  
  print(N1)
  print(N2)
  # ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
  # geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
  
  K <- nrow(traps)
  
  CH_1 <- matrix(NA, nrow = N1, ncol = K)
  for (i in 1:N1) {
    for (k in 1:K) {
      p <- p0_1 * exp(-(1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
      CH_1[i,k] <- rbinom(1, S_usage[k], p)
    }
  }
  indexesCaughtIndividuals1 <- which(apply(CH_1,1,sum) != 0)
  indexesUncaughtIndividuals1 <- which(apply(CH_1,1,sum) == 0)
  D1 <- length(indexesCaughtIndividuals1)
  CH_1 <- CH_1[indexesCaughtIndividuals1,]
  
  CH_2 <- matrix(NA, nrow = N2, ncol = K)
  for (i in 1:N2) {
    for (k in 1:K) {
      p <- p0_2 * exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
      CH_2[i,k] <- rbinom(1, S_usage[k], p)
    }
  }
  indexesCaughtIndividuals2 <- which(apply(CH_2,1,sum) != 0)
  indexesUncaughtIndividuals2 <- which(apply(CH_2,1,sum) == 0)
  D2 <- length(indexesCaughtIndividuals2)
  CH_2 <- CH_2[indexesCaughtIndividuals2,]
  
  CH_leopards <- CH_1
  CH_tigers <- CH_2
  
  N1_0 <- N1
  N2_0 <- N2
  
  s1_0 <- matrix(NA, nrow = N1, ncol = 2)
  s2_0 <- matrix(NA, nrow = N2, ncol = 2)
  
  s1_0[seq_len(D1),] <- s1[indexesCaughtIndividuals1,]
  if(N1 > D1){
    s1_0[(D1 + 1):N1,] <- s1[indexesUncaughtIndividuals1,]
  }
  s2_0[seq_len(D2),] <- s2[indexesCaughtIndividuals2,]
  if(N2 > D2){
    s2_0[(D2 + 1):N2,] <- s2[indexesUncaughtIndividuals2,]  
  }
  
  print(D1)
  print(D2)
}

# CLEAN DATA --------------------------------------------------------------

D1 <- nrow(CH_leopards)
D2 <- nrow(CH_tigers)

K <- ncol(CH_leopards) # number of traps

CH1 <- CH_leopards
CH2 <- CH_tigers

# PRIOR -------------------------------------------------------------------

Nmax <- 1200

# proposal interaction parameters
{
  epsilon_beta <- 70^2
}

# intensity parameter
{
  meanBeta <- N1
  varBeta <- 10000
  b_beta <- meanBeta / varBeta
  a_beta <- meanBeta * b_beta  
  
  # ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,400))
}

# proposal value
sigma_prop <- .01
Sigma_newpoint <- cov(cbind(trapsX, trapsY))
lambda_newPoints1 <- 4
lambda_newPoints2 <- 4

# capture probability params
{
  a_p <- 1
  b_p <- 1
  
  # leopards (.5 to 2.5)
  sigma_0_1 <- 1.5 / traps_meansd
  sd_sigma_1 <- 2 / traps_meansd
  
  # tigers ( 1.5 to 3.5)
  sigma_0_2 <- 2.5 / traps_meansd
  sd_sigma_2 <- 2 / traps_meansd
  
  # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1))
  # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_2, sd = sd_sigma_2))
  
  sigma_p0_1_prop <- .0002
  sigma_sigma1_prop <- .005
  
  sigma_p0_2_prop <- .0001
  sigma_sigma2_prop <- .01
}

# points proposal parameters
{
  probMixture <- c(.95,.95)
}

# MCMC --------------------------------------------------------------------

pointsToAccept <- 1000#nburn <- 0
niter <- 5000
nchain <- 1
nthin <- 1
nburn <- 0

# grid fornew points density
{
  gridLength_x <- seq(min(trapsX) - .05, 
                      max(trapsY) + .05, length.out = 100)
  gridLength_y <- seq(min(trapsY) - .05, 
                      max(trapsY) + .05, length.out = 100)
}

# variables for output
{
  sigma_output <- array(NA, dim = c(nchain, niter, 2))
  # s_output <- array(NA, dim = c(nchain, niter, 2, Nmax, 2))
  N_output <- array(NA, dim = c(nchain, niter, 2))
  p0_output <- array(NA, dim = c(nchain, niter, 2))
  params_output <- array(NA, dim = c(nchain, niter, 5))
  
  # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
  # params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
  
  # acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
}

for(chain in 1:nchain) {
  
  # starting values
  # if(nburn != 0){
  
  # home centers
  {
    if(usingTrueData){
      N1 <- 2 * D1
      N2 <- 2 * D2
    } else {
      N1 <- N1_0
      N2 <- N2_0
    }
    
    CH1 <- rbind(CH1[1:D1,],
                 matrix(0, nrow = N1 - D1, ncol = K))
    
    CH2 <- rbind(CH2[1:D2,],
                 matrix(0, nrow = N2 - D2, ncol = K))
    
    if(usingTrueData){
      
      s1 <- matrix(NA, nrow = Nmax, ncol = 2)
      for (i in 1:N1) {
        if(i <= D1){ # assign the home center to the average location of the traps where the indiv. was captured
          s1[i,1] <- mean(trapsX[CH1[i,] != 1])
          s1[i,2] <- mean(trapsY[CH1[i,] != 1])#apply(trapsX[CH1[i,] == 1, 2:3], 2, mean)
          s1[i,] <- mvrnorm(1, mu = s1[i,], diag(.005, nrow = 2))
        } else { # otherwise a random point
          s1[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
        }
      }
      
      s2 <- matrix(NA, nrow = Nmax, ncol = 2)
      for (i in 1:N2) {
        if(i <= D2){ # assign the home center to the average location of the traps where the indiv. was captured
          s2[i,1] <- mean(trapsX[CH2[i,] != 1])
          s2[i,2] <- mean(trapsY[CH2[i,] != 1])
          s2[i,] <- mvrnorm(1, mu = s2[i,], diag(.005, nrow = 2))
        } else { # otherwise a random point
          s2[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
        }
      }
      
    } else {
      
      s1 <- matrix(NA, nrow = Nmax, ncol = 2)
      s2 <- matrix(NA, nrow = Nmax, ncol = 2)
      
      s1[1:N1,] <- s1_0
      s2[1:N2,] <- s2_0
      
    }
    
  }
  
  # capture probability parameter
  {
    if(usingTrueData){
      
      sigma_1 <- sigma_0_1
      
      c_i <- sapply(1:N1, function(i){
        sapply(1:K, function(k){
          exp(- (1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
        })
      })
      
      p0_1 <- rbeta(1, sum(CH1) + 1, sum(mean(S_usage) - CH1) + 1) / mean(c_i)
      
      sigma_2 <- sigma_0_2
      
      c_i <- sapply(1:N2, function(i){
        sapply(1:K, function(k){
          exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
        })
      })
      
      p0_2 <- rbeta(1, sum(CH2) + 1, sum(mean(S_usage) - CH2) + 1) / mean(c_i)
      
    } else {
      p0_1 <- p0_1_true
      sigma_1 <- sigma_1_true
      p0_2 <- p0_2_true
      sigma_2 <- sigma_2_true
    }
    
  }
  
  # point process parameters
  {
    if(usingTrueData){
      theta1 <- .002
      theta2 <- .002
      theta12 <- .002  
      beta1 <- N1
      beta2 <- N2
    } else {
      theta1 <- theta1_true
      theta2 <- theta2_true
      theta12 <- theta12_true
      beta1 <- N1
      beta2 <- N2
    }
    
  }
  
  # new proposal parameters
  {
    
    numAcceptedPoints <- c(0,0)
    acceptedPoints1 <- matrix(NA, pointsToAccept, 2)
    acceptedPoints2 <- matrix(NA, pointsToAccept, 2)
    
    # compute normalizing constant for uniform
    {
      M <- 10000
      
      normConstUnif <- mean(sapply(1:M, function(i){
        x <- c(runif(1, a1, b1), runif(1, a2, b2))
        checkPointIsInRegionTraps(x, allTraps, R_traps)
      }))
      
      S_area <- (b1 - a1) * (b2 - a2) * normConstUnif
      
      normConst1 <- c(normConstUnif)
      normConst2 <- c(normConstUnif)
      
      mixtureMeans1 <- matrix(0, nrow = 0, ncol = 2)
      mixtureSd1 <- matrix(0, nrow = 0, ncol = 2)
      mixtureMeans2 <- matrix(0, nrow = 0, ncol = 2)
      mixtureSd2 <- matrix(0, nrow = 0, ncol = 2)
      
    }
    
  }
  
  # }
  
  # output variables
  {
    sigma_iter <- array(NA, dim = c(niter, 2))
    # s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
    N_iter <- matrix(NA, nrow = niter, ncol = 2)
    p0_iter <- array(NA, dim = c(niter, 2))
    params_iter <- matrix(NA, nrow = niter, ncol = 5)
    
    # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
    params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
    
    acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
    theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
  }
  
  iter <- 1
  outOfBurnIn <- F
  mixture1Updated <- F
  mixture2Updated <- F
  
  while(iter < niter){
    
    if(!outOfBurnIn){
      print(paste0("Chain = ",chain," / AcceptedPoints1 = ",numAcceptedPoints[1]," / AcceptedPoints2 = ",numAcceptedPoints[2])) 
    } else if(((iter - nburn)/nthin) %% 5 == 0){
      print(paste0("Chain = ",chain," / Iteration = ",(iter - nburn)/nthin))
    }
    
    print(paste0("N1 = ",N1," - N2 = ",N2))
    
    # UPDATE P ----------------------------------------------------------------
    
    list_p <- update_p(p0_1, sigma_1, p0_2, sigma_2,
                       s1, s2, CH1, CH2, trapsX, trapsY,
                       S_usage, N1, N2,
                       sigma_0_1, sd_sigma_1,
                       sigma_0_2, sd_sigma_2,
                       sigma_p0_1_prop, sigma_p0_2_prop,
                       sigma_sigma1_prop, sigma_sigma2_prop,
                       a_p, b_p)
    p0_1 <- list_p$p0_1
    p0_2 <- list_p$p0_2
    sigma_1 <- list_p$sigma_1
    sigma_2 <- list_p$sigma_2
    
    # UPDATE HOME CENTERS -----------------------------------------------
    
    # if(iter < nburn){
    # list_s <-  update_s_cpp(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
    #                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
    #                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
    #                         allTraps, R_traps,
    #                         a1, b1, a2, b2)  
    
    if(outOfBurnIn){
      
      # list_s <-  update_s_cpp(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
      #                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
      #                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
      #                         allTraps, R_traps,
      #                         a1, b1, a2, b2)
      list_s <- update_s_cpp_new(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
                                 CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                                 Sigma_newpoint, lambda_newPoints1, lambda_newPoints2, probMixture[1],
                                 probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
                                 normConst1, normConst2, allTraps, R_traps, S_area, S_usage, trapsX,
                                 trapsY, a1, b1, a2, b2)
    } else {
      list_s <-  update_s_cpp_withacceptance(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
                                             CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                             S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
                                             allTraps, R_traps,
                                             a1, b1, a2, b2)
    }
    
    s1 <- list_s$s1
    s2 <- list_s$s2
    N1 <- list_s$N1
    N2 <- list_s$N2
    
    CH1 <- rbind(CH1[1:D1,],
                 matrix(0, nrow = N1 - D1, ncol = K))
    
    CH2 <- rbind(CH2[1:D2,],
                 matrix(0, nrow = N2 - D2, ncol = K))
    
    if(!outOfBurnIn){
      if(numAcceptedPoints[1] < pointsToAccept & list_s$isPoint1Accepted){
        acceptedPoints1[numAcceptedPoints[1] + 1,] <- list_s$acceptedPoint1
        numAcceptedPoints[1] <- numAcceptedPoints[1] + 1
      }
      
      if(numAcceptedPoints[2] < pointsToAccept & list_s$isPoint2Accepted){
        acceptedPoints2[numAcceptedPoints[2] + 1,] <- list_s$acceptedPoint2
        numAcceptedPoints[2] <- numAcceptedPoints[2] + 1
      }  
    }
    
    # UPDATE BETA ---------
    
    beta1 <- rgamma(1, a_beta + N1, b_beta + 1)
    beta2 <- rgamma(1, a_beta + N2, b_beta + 1)
    
    # WRITE RESULTS ---------------
    
    outOfBurnIn <- numAcceptedPoints[1] >= pointsToAccept & numAcceptedPoints[2] >= pointsToAccept
    
    if(outOfBurnIn){
      
      # trueIter <- (iter - nburn)/nthin
      trueIter <- iter
      
      # s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
      # s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
      p0_iter[trueIter,] <- c(p0_1, p0_2)
      sigma_iter[trueIter,] <- c(sigma_1, sigma_2)
      N_iter[trueIter,] <- c(N1,N2)  
      params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
      # papangelou_density_iter[trueIter,,,] <- computeNewPointsDensity(gridLength_x, gridLength_y,
      # s1, s2, N1, N2,
      # theta1,  theta2,  theta12,
      # beta1,  beta2)
      
      iter <- iter + 1
    }
    
    # update proposals
    {
      if(outOfBurnIn & !mixture1Updated){
        
        numCenters_all <- 5:30
        kmeans_totss <- rep(length(numCenters_all))
        
        for (l in seq_along(numCenters_all)) {
          
          numCenters <- numCenters_all[l]
          
          kmeansfit <- kmeans(acceptedPoints1, centers = numCenters)
          
          kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints1)) 
        }
        
        best_k <- numCenters_all[which.min(kmeans_totss)]
        kmeansfit <- kmeans(acceptedPoints1, centers = best_k)
        
        mixtureMeans1 <- matrix(NA, nrow = best_k, ncol = 2)
        mixtureSd1 <- matrix(NA, nrow = best_k, ncol = 2)
        normConstNormals1 <- rep(NA, best_k)
        
        for (i in 1:best_k) {
          mixtureMeans1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, mean)
          mixtureSd1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, sd)
          
          normConstNormals1[i] <- mean(sapply(1:M, function(j){
            x <- c(rnorm(1, mixtureMeans1[i,1], mixtureSd1[i,1]), 
                   rnorm(1, mixtureMeans1[i,2], mixtureSd1[i,2]))
            checkPointIsInRegionTraps(x, allTraps, R_traps)
          }))
        }
        
        normConst1 <- c(normConstUnif, normConstNormals1)
        
        mixture1Updated <- T
      }
      
      if(outOfBurnIn & !mixture2Updated){
        
        numCenters_all <- 5:30
        kmeans_totss <- rep(length(numCenters_all))
        
        for (l in seq_along(numCenters_all)) {
          
          numCenters <- numCenters_all[l]
          
          kmeansfit <- kmeans(acceptedPoints2, centers = numCenters)
          
          kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints2)) 
        }
        
        best_k <- numCenters_all[which.min(kmeans_totss)]
        kmeansfit <- kmeans(acceptedPoints2, centers = best_k)
        
        mixtureMeans2 <- matrix(NA, nrow = best_k, ncol = 2)
        mixtureSd2 <- matrix(NA, nrow = best_k, ncol = 2)
        normConstNormals2 <- rep(NA, best_k)
        
        for (i in 1:best_k) {
          mixtureMeans2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, mean)
          mixtureSd2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, sd)
          
          normConstNormals2[i] <- mean(sapply(1:M, function(j){
            x <- c(rnorm(1, mixtureMeans2[i,1], mixtureSd2[i,1]), 
                   rnorm(1, mixtureMeans2[i,2], mixtureSd2[i,2]))
            checkPointIsInRegionTraps(x, allTraps, R_traps)
          }))
        }
        
        normConst2 <- c(normConstUnif, normConstNormals2)
        
        mixture2Updated <- T
      }    
    }
    
    # intermediate plots
    {
      if(iter %% 250 == 0){
        # plotCurrent1 <- qplot(1:iter, params_iter[1:iter,1], geom =  "line")
        # plotCurrent2 <- qplot(1:iter, log(params_iter[1:iter,3]), geom =  "line")
        # plotCurrent3 <- qplot(1:iter, log(params_iter[1:iter,5]), geom =  "line")
        # ggsave(filename = paste0("plot_param1_iter",iter,".jpg"), plotCurrent1)
        # ggsave(filename = paste0("plot_param2_iter",iter,".jpg"), plotCurrent2)
        # ggsave(filename = paste0("plot_param3_iter",iter,".jpg"), plotCurrent3)
      }  
    }
    
  }
  
  # Write results in MCMC output
  {
    sigma_output[chain,,] <- sigma_iter
    N_output[chain,,] <- N_iter
    p0_output[chain,,] <- p0_iter
    params_output[chain,,] <- params_iter
  }
  
}

qplot(1:niter, N_iter[,1]) + geom_hline(aes(yintercept = N1_0))
qplot(1:niter, N_iter[,2]) + geom_hline(aes(yintercept = N2_0))

# CHAIN DIAGNOSTICS -------------------------------------------------------

# plot chains
{
  plotVar <- function(variable_output){
    
    nchain <- nrow(variable_output)
    niter <- ncol(variable_output)
    
    # beta_psi_output <- matrix(beta_psi_output[,,1], nrow = nchain, ncol = niter)
    
    variable_output_long <- melt(variable_output)
    
    ggplot() + geom_line(data = variable_output_long, aes(x = Var2, y = value, group = Var1, 
                                                          color = factor(Var1))) + 
      xlab("Iterations") + ylab("Value") +
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black"),
            legend.position = "none") 
    
  }
  
  
}

# diagnostics p
{
  plotVar(p0_output[,,1,drop = F])
  plotVar(p0_output[,,2])
  plotVar(sigma_output[,,1])
  plotVar(sigma_output[,,2])
}

# diagnostics N
{
  plotVar(N_output[,,1]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
  plotVar(N_output[,,2]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
}

# diagnostics theta
{
  plotVar(params_output[,,1]) 
  plotVar(params_output[,,2]) 
  plotVar(log(params_output[,,3]))
  plotVar(log(params_output[,,4]))
  plotVar(log(params_output[,,5]))
}
# DIAGNOSTICS -------------------------------------------------------------

# diagnostics p
{
  plot(1:niter, p0_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = p0_1_true))
  plot(1:niter, p0_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = p0_2_true))
  plot(1:niter, sigma_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = sigma_1_true))
  plot(1:niter, sigma_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = sigma_2_true))
}

# diagnostics N
{
  plot(1:niter, N_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
  plot(1:niter, N_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = N2_0))
}

# diagnostics params
{
  plot(1:niter, params_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = beta1_true))
  plot(1:niter, params_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = beta2_true))
  plot(1:niter, params_iter[,3]) #+ geom_hline(data = NULL, aes(yintercept = theta1_true))
  plot(1:niter, params_iter[,4]) #+ geom_hline(data = NULL, aes(yintercept = theta2_true))
  plot(1:niter, params_iter[,5]) #+ geom_hline(data = NULL, aes(yintercept = theta12_true))
}

params_iter_long <- melt(params_iter)

CI_params <- apply(params_iter, 2, function(x){
  quantile(x, probs = c(.05,.5,.95), na.rm = T)
})
c(D1,D2)
intParamsPlot <- ggplot(data = NULL, aes(x = c("theta1","theta2","theta12"), y = CI_params[2,3:5], 
                                         ymin = CI_params[1,3:5], 
                                         ymax = CI_params[3,3:5])) + 
  geom_point(stat = "identity", size = 4, shape = 1) +
  geom_errorbar() + 
  # geom_point(aes(x = c("theta1","theta2","theta12"), y = c(theta1_true, theta2_true, theta12_true)), shape = 4, size = 4) + 
  scale_x_discrete(name = "Parameter", labels = c(expression(theta[1]),expression(theta[12]),expression(theta[2]))) +
  scale_y_continuous(name = "", limits = c(0, 0.15)) +
  geom_text(data = NULL, aes(x = 2, y = .145), 
            label = expression(paste(D[1]," = ",32, " - ", D[2]," = ",37)),
            size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black")) 
intParamsPlot
setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Simulation")

ggsave(filename = paste0("beta = ",beta1_true, "- p0 = ",p0_true,".jpeg"), intParamsPlot,
       width = 6.36, height = 5.25, units = "in")


qplot(1:niter, params_iter[,1])
qplot(1:niter, params_iter[,2])
qplot(1:niter, params_iter[,5])

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

qplot(1:niter, p0_iter)

qplot(1:niter, N_iter[,1])
qplot(1:niter, N_iter[,2])

qplot(1:niter, s_iter[,1,1,1])
qplot(1:niter, s_iter[,1,21,2])
qplot(1:niter, s_iter[,2,11,2])

ggplot(data = params_iter, aes(x = c))

# PLOTS -------------------------------------------------------------------

# plots for gamma
{
  CI_params <- apply(params_iter, 2, function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  plotgamma <- ggplot(data = NULL, aes(x = reorder(factor(c("theta1","theta2","theta12")),c(1,3,2)), y = CI_params[2,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[3,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c(expression(theta[1]),expression(theta[1][2]),expression(theta[2])),
                     name = "Interaction parameters") +
    scale_y_continuous(name = "") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "gamma.jpeg",  plotgamma,
         width = 6.36, height = 5.25, units = "in")
}

# plots for N
{
  plotN <- ggplot(data = NULL, aes(x = N_iter[,1])) + 
    geom_bar(aes(y = (..count..)/sum(..count..)), 
             fill = "cornsilk", color = "black", stat = "count", width = .2, size = .5) +
    scale_x_continuous(breaks = sort(unique(N_iter[,1])), 
                       labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       minor_breaks = sort(unique(N_iter[,1])),   
                       name = expression(paste("Number of individuals per ",km^2))) +
    scale_y_continuous(name = "Density") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13, face = "bold"))
  plotN
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "Popsize1.jpeg",  plotN,
         width = 6.36, height = 5.25, units = "in")
  
  plotN2 <- ggplot(data = NULL, aes(x = N_iter[,2])) + 
    geom_bar(aes(y = (..count..)/sum(..count..)), 
             fill = "cornsilk", color = "black", stat = "count", width = .2, size = .5) +
    scale_x_continuous(breaks = sort(unique(N_iter[,2])), 
                       labels =  round(sort(unique(N_iter[,2])) /  925,3),
                       minor_breaks = sort(unique(N_iter[,2])),   
                       name = expression(paste("Number of individuals per ",km^2))) +
    scale_y_continuous(name = "Density") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13, face = "bold"))
  plotN2
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "Popsize2.jpeg",  plotN2,
         width = 6.36, height = 5.25, units = "in")
  
  
}


# CREATE BOUNDARIES -------------------------------------------------------

pointsBoundaries <- matrix(NA, nrow = 10, ncol = 2)
pointsBoundaries[1,] <- c(-2.4, 0.5)
pointsBoundaries[2,] <- c(-1.3, 2.65)
pointsBoundaries[3,] <- c(.85, -2.65)
pointsBoundaries[4,] <- c(2.2, -1.35)
pointsBoundaries[5,] <- c(1.8, .6)
pointsBoundaries <- pointsBoundaries[complete.cases(pointsBoundaries),]

coefsLine <- function(x1, y1, x2, y2){
  
  m <- (y2 - y1) / (x2 - x1)
  
  q <- y1 - m * x1
  
  list("m" = m, "q" = q)
}

line1 <- function(x){
  
  m <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[2,1], pointsBoundaries[2,2])$m
  
  q <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[2,1], pointsBoundaries[2,2])$q
  
  return(m * x + q)
}

line2 <- function(x){
  
  m <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$m
  
  q <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$q
  
  return(m * x + q)
}

line3 <- function(x){
  
  m <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$m
  
  q <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$q
  
  return(m * x + q)
}

line4 <- function(x){
  
  m <- coefsLine(pointsBoundaries[2,1], pointsBoundaries[2,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$m
  
  q <- coefsLine(pointsBoundaries[2,1], pointsBoundaries[2,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$q
  
  return(m * x + q)
}

line5 <- function(x){
  
  m <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$m
  
  q <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$q
  
  return(m * x + q)
}

{
  ggplot() +
    geom_point(data = NULL, aes(x = traps$x, y = traps$y)) + 
    scale_x_continuous(name = "Longitude") + 
    scale_y_continuous(name = "Latitude") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                  fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                  fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
    stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                  fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                  fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                  fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  
  
}




# HOME RANGES -------------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

sd_s1 <- apply(s_iter[,1,1:D1,], c(2,3), sd)
sd_s2 <- apply(s_iter[,2,1:D2,], c(2,3), sd)

sd_s1_mean <- apply(sd_s1, 1, mean)
sd_s2_mean <- apply(sd_s2, 1, mean)

ggplot() + geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2], size = as.numeric(4 * sd_s1_mean)), color = "red", alpha = .3) + 
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2], size = 4 * sd_s2_mean), color = "black", alpha = .3) + 
  scale_size_continuous() +
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  



# DENSITY PLOT ------------------------------------------------------------

N1 <- N_iter[niter,1]
N2 <- N_iter[niter,2]
s1_mean <- apply(s_iter[,1,1:N1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:N2,], c(2,3), mean)

ggplot() + geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]), color = "black") 

gridLength_x <- seq(min(as.vector(s_iter[,,,1]), na.rm = T)  - .05, 
                    max(as.vector(s_iter[,,,1]), na.rm = T) + .05, length.out = 25)
gridLength_y <- seq(min(as.vector(s_iter[,,,2]), na.rm = T) - .05, 
                    max(as.vector(s_iter[,,,2]), na.rm = T) + .05, length.out = 25)
pointsInGrid <- array(0, dim = c(niter, 2, length(gridLength_x) - 1, length(gridLength_y) - 1))

for (iter in 1:niter) {
  print(iter)
  N1 <- N_iter[niter,1]
  for (i in 1:N1) {
    cellX <- findInterval(s_iter[iter,1,i,1], gridLength_x)
    cellY <- findInterval(s_iter[iter,1,i,2], gridLength_y)
    pointsInGrid[iter, 1, cellX, cellY] <- pointsInGrid[iter, 1, cellX, cellY] + 1
  }
  N2 <- N_iter[niter,2]
  for (i in 1:N2) {
    cellX <- findInterval(s_iter[iter,2,i,1], gridLength_x)
    cellY <- findInterval(s_iter[iter,2,i,2], gridLength_y)
    pointsInGrid[iter, 2, cellX, cellY] <- pointsInGrid[iter, 2, cellX, cellY] + 1
  }
}

binwidth <- gridLength_x[2] - gridLength_x[1]
midPoints_x <- (gridLength_x[-1] + gridLength_x[-length(gridLength_x)]) / 2 
midPoints_y <- (gridLength_y[-1] + gridLength_y[-length(gridLength_y)]) / 2 
allGrid <- expand.grid(midPoints_x, midPoints_y)
meanDensityGrid <- apply(pointsInGrid, c(2,3,4), mean)
meanDensityGrid_vec <- as.vector(meanDensityGrid[1,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensityGrid_vec), 
             size = 100 * binwidth, shape = 15) + 
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  


# PAPANGELOU DENSITY ------------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

papangelou_density_mean <- apply(papangelou_density_iter, c(2,3,4), mean)

allGrid <- expand.grid(gridLength_x, gridLength_y)
binwidth <- gridLength_x[2] - gridLength_x[1]

meanDensity1Grid_vec <- as.vector(papangelou_density_mean[1,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensity1Grid_vec), 
             size = 100 * binwidth, shape = 15) + 
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 50 * binwidth, color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 50 * binwidth, color = "red") +
  scale_x_continuous(name = "Longitude") + 
  scale_color_continuous(low = "black", high = "white") +
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "blue") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "blue")  

meanDensity2Grid_vec <- as.vector(papangelou_density_mean[2,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensity2Grid_vec), 
             size = 100 * binwidth, shape = 15) + 
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 50 * binwidth, color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 50 * binwidth, color = "red") +
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  scale_color_continuous(low = "black", high = "white") +
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  

# HOME CENTERS  -----------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

s1_sd <- apply(s_iter[,1,1:D1,], c(2,3), sd)
s2_sd <- apply(s_iter[,2,1:D2,], c(2,3), sd)

binwidth <- gridLength_x[2] - gridLength_x[1]

meanTrapsX <- 0
sdTrapsX <-  0.1035411
meanTrapsY <- 0
sdTrapsY <- 0.07364998

plotHomeRange <- ggplot() +
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 35 * binwidth,  color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 35 * binwidth,  color = "red") +
  geom_point(data = NULL, aes(x = trapsX, y = trapsY), size = 1, shape = 3) +
  # geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2], size = 4 * s1_sd[,1]),
  # shape = 1, color = "blue") +
  # geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2], size = 4 * s2_sd[,1]),
  # shape = 1, color = "red") +
  scale_x_continuous(name = "Longitude", breaks = -2:2, labels = round((-2:2)*sdTrapsX + meanTrapsX,2)) + 
  scale_y_continuous(name = "Latitude", breaks = -2:2, labels = round((-2:2)*sdTrapsY + meanTrapsY,2)) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "black") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black")) 

plotHomeRange
setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")

ggsave(filename = "homerangesplot.jpeg", plotHomeRange,
       width = 5, height = 5, units = "in")

