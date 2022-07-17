library(MASS); library(reshape2); library(ggplot2)
library(foreach); library(doParallel); library(coda)
# library(SpatialFunctionsCR);
library(Rcpp); library(RcppArmadillo)

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
# stopCluster(cl)

setwd("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simplified model")
sourceCpp("code.cpp")
source('functions.r', echo=TRUE)

# CREATE NEW TRAPS -------

a1 <- 0
b1 <- 1
a2 <- 0
b2 <- 1

# equispaced
{
  numTraps <- 7
  trapsX_grid <- seq(a1, b1, length.out = numTraps + 2)[-(numTraps+2)][-1]
  trapsY_grid <- seq(a2, b2, length.out = numTraps + 2)[-(numTraps+2)][-1]
  traps <- expand.grid(trapsX_grid, trapsY_grid)
}

# simulation soft
{
  list_traps <- simulate_softcore_cpp(theta = .05, beta = 4000, niter = 1000, 
                                      Sigma_prop = diag(.05, nrow = 2), Nmax = 200,
                                      lambda = 100, a1, b1, a2, b2)  
  traps <-  list_traps$data
}

# simulation soft
{
  set.seed(4)#
  list_traps <- simulate_cond_softcore_cpp(theta = exp(-100), beta = 100, niter = 1000, 
                                           Sigma_prop = diag(.05, nrow = 2), N = 200,
                                           lambda = 100, a1, b1, a2, b2)  
  traps <-  list_traps$data
}

ggplot(data = NULL, aes(x = traps[,1], y = traps[,2])) + geom_point()

trapsX <- traps[,1]
trapsY <- traps[,2]

trapsX_mean <- mean(trapsX)
trapsY_mean <- mean(trapsY)
trapsX_sd <- sd(trapsX)
trapsY_sd <- sd(trapsY)


# traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2
# 
# trapsX <- (trapsX - mean(trapsX)) / traps_meansd
# trapsY <- (trapsY - mean(trapsY)) / traps_meansd

allTraps <- cbind(trapsX, trapsY)
# R_traps <- .75
# 
# a1 <- min(trapsX) - R_traps
# b1 <- max(trapsX) + R_traps
# a2 <- min(trapsY) - R_traps
# b2 <- max(trapsY) + R_traps

K <- length(trapsX)

# SIMULATED DATA ----------------------------------------------------------

usingTrueData <- F

# starting params
{
  beta <- 300
  theta <- exp(-1000)
  
  p0 <- .02
  sigma <- .05
  
  S_usage <- rep(40, K)
  
  beta_true <- beta
  p0_true <- p0
  sigma_true <- sigma
  theta_true <- theta
}

# data simulation
{
  # N <- rpois(1, beta)
  # s <- t(sapply(1:N, function(i){
  #   c(runif(1, a1, b1), runif(1, a2, b2))
  # }))
  
  
  list_s <- simulate_softcore_cpp(theta, beta, niter = 1000, 
                                  Sigma_prop = diag(.05, nrow = 2), Nmax = 1000,
                                  lambda = beta, a1, b1, a2, b2)  
  s <-  list_s$data
  N <-  list_s$N
  print(N)
  s <- s[1:N,]
  
  ggplot(data = NULL, aes(x = s[,1], y = s[,2])) + geom_point()
  
  CH <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N) {
    for (k in 1:K) {
      p <- p0 * exp(-(1 / (2 * sigma^2)) * ((s[i,1] - trapsX[k])^2 + (s[i,2] - trapsY[k])^2))
      CH[i,k] <- rbinom(1, S_usage[k], p)
    }
  }
  indexesCaughtIndividuals <- which(apply(CH,1,sum) != 0)
  indexesUncaughtIndividuals <- which(apply(CH,1,sum) == 0)
  D <- length(indexesCaughtIndividuals)
  CH <- CH[indexesCaughtIndividuals,]
  
  N_0 <- N
  
  s_0 <- matrix(NA, nrow = N, ncol = 2)
  
  s_0[seq_len(D),] <- s[indexesCaughtIndividuals,]
  if(N > D){
    s_0[(D + 1):N,] <- s[indexesUncaughtIndividuals,]
  }
  
  print(D)
}

# plot
{
  ggplot() + geom_point(data = NULL, aes(x = s[indexesCaughtIndividuals,1],
                                         y = s[indexesCaughtIndividuals,2]), color = "red") +  
    geom_point(data = NULL, aes(x = s[indexesUncaughtIndividuals,1],
                                y = s[indexesUncaughtIndividuals,2]), color = "blue") + 
    geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "black", size = 2)
}

# PRIOR -------------------------------------------------------------------

Nmax <- 1200

# intensity parameter
{
  meanBeta <- N
  varBeta <- 10000
  b_beta <- meanBeta / varBeta
  a_beta <- meanBeta * b_beta  
  
  # ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,400))
}

# proposal value
sigma_prop <- .01
# lambda_newPoints1 <- 3
# lambda_newPoints2 <- 3

# capture probability params
{
  a_p <- 1
  b_p <- 1
  
  sigma_0 <- 1.5 #/ traps_meansd
  sd_sigma <- 2 #/ traps_meansd
  
  sigma_p0_prop <- .0002
  sigma_sigma_prop <- .005
}

# points proposal parameters
{
  probMixture <- c(.95,.95)
}

# MCMC --------------------------------------------------------------------

niter <- 20000
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
      N <- 2 * D
    } else {
      N <- N_0
    }
    
    CH <- rbind(CH[1:D,],
                matrix(0, nrow = N - D, ncol = K))
    
    if(usingTrueData){
      
      s <- matrix(NA, nrow = Nmax, ncol = 2)
      for (i in 1:N) {
        if(i <= D){ # assign the home center to the average location of the traps where the indiv. was captured
          s[i,1] <- mean(trapsX[CH[i,] != 1])
          s[i,2] <- mean(trapsY[CH[i,] != 1])#apply(trapsX[CH1[i,] == 1, 2:3], 2, mean)
          s[i,] <- mvrnorm(1, mu = s[i,], diag(.005, nrow = 2))
        } else { # otherwise a random point
          s[i,] <- proposeNewPointPP(a1, b1, a2, b2)
        }
      }
      
    } else {
      
      s <- matrix(NA, nrow = Nmax, ncol = 2)
      
      s[1:N,] <- s_0
      
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
      p0 <- p0_true
      sigma <- sigma_true
    }
    
  }
  
  # point process parameters
  {
    if(usingTrueData){
      beta <- beta_true
    } else {
      beta <- beta_true
      theta <- theta_true
    }
    
  }
  
  # new proposal parameters
  {
    lambda_newPoints <- 5
    
    numAcceptedPoints <- c(0,0)
    acceptedPoints <- matrix(NA, pointsToAccept, 2)
    N_burnin <- rep(NA, pointsToAccept * 10)
    
    # compute normalizing constant for uniform
    {
      M <- 10000
      
      normConstUnif <- 1
      
      S_area <- (b1 - a1) * (b2 - a2) 
      
      normConst <- c(normConstUnif)
      
      mixtureMeans <- matrix(0, nrow = 0, ncol = 2)
      mixtureSd <- matrix(0, nrow = 0, ncol = 2)
      
    }
    
  }
  
  # output variables
  {
    sigma_iter <- rep(NA, niter)
    # s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
    N_iter <- rep(NA, niter)
    p0_iter <- rep(NA, niter)
  }
  
  iter <- 0
  burnin_it <- 1
  outOfBurnIn <- T
  
  et1 <- Sys.time()
  
  while(iter < niter){
    
    if(!outOfBurnIn){
      print(paste0("Chain = ",chain," / AcceptedPoints = ",numAcceptedPoints[1])) 
    } else if(((iter - nburn)/nthin) %% 1 == 0){
      print(paste0("Chain = ",chain," / Iteration = ",(iter - nburn)/nthin))
    }
    
    print(paste0("N = ",N))
    
    # UPDATE P ----------------------------------------------------------------
    
    # list_p <- update_p(p0, sigma, s, CH, trapsX, trapsY, 
    #                    S_usage, N, sigma_0, sd_sigma_1,
    #                    sigma_p0_prop, sigma_sigma_prop, a_p, b_p)
    # p0 <- list_p$p0
    # sigma <- list_p$sigma
    
    # UPDATE HOME CENTERS -----------------------------------------------
    
    if(!outOfBurnIn){
      s <- update_s(s, beta, theta, CH, N, D, p0, sigma, N, Sigma_prop = diag(.1, nrow = 2),
                    Sigma_newpoint = diag(.1, nrow = 2), allTraps, 
                    S_usage, trapsX, trapsY, a1, b1, a2, b2)  
    } else {
      s <- update_s(s, beta, theta, CH, N, D, p0, sigma, N, Sigma_prop = diag(.1, nrow = 2),
                    Sigma_newpoint = diag(.1, nrow = 2), allTraps, 
                    S_usage, trapsX, trapsY, a1, b1, a2, b2)  
    }
    
    # UPDATE N -----------------------------------------------
    
    list_s <- update_N_withacceptance(s, beta, theta, CH, N, D, p0, sigma, S_usage, trapsX,
                                      trapsY, Sigma_prop = diag(.1, nrow = 2),
                                      Sigma_newpoint = diag(.1, nrow = 2), a1, b1, a2, b2)
    
    s <- list_s$s
    N <- list_s$N
    
    CH <- rbind(CH[1:D,],
                matrix(0, nrow = N - D, ncol = K))
    
   
    # UPDATE BETA ---------
    
    beta <- rgamma(1, a_beta + N, b_beta + 1)
    
    # WRITE RESULTS ---------------
    
    if(outOfBurnIn){
      
      iter <- iter + 1
      
      # trueIter <- (iter - nburn)/nthin
      trueIter <- iter
      
      p0_iter[trueIter] <- p0
      sigma_iter[trueIter] <- sigma
      N_iter[trueIter] <- N  
      
      
    }
    
    
    
  }
  
  # Write results in MCMC output
  {
    sigma_output[chain,,] <- sigma_iter
    N_output[chain,,] <- N_iter
    p0_output[chain,,] <- p0_iter
    # params_output[chain,,] <- params_iter
  }
  
}

et1 <- Sys.time()

diffSecs <- as.numeric(difftime(et1,st1), units = "secs")

effectiveSize(mcmc(N_iter)) / diffSecs # 0.05197856 
effectiveSize(mcmc(N_iter[(1:(niter/2))*2])) / diffSecs

qplot(1:niter, N_iter, geom = "line") + geom_hline(aes(yintercept = N_0))
qplot(1:100, N_iter[1:100], geom = "line") + geom_hline(aes(yintercept = N_0))

ggplot(data = NULL, aes(x = mixtureMeans[,1],
                        y = mixtureMeans[,2])) + geom_point(color = "red") +
  geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "black", size = 2)

# DIAGNOSTICS -------------------------------------------------------------

# burn_in_iter <- sum(!is.na(N_burnin))
# 
# qplot(1:burn_in_iter, N_burnin[1:burn_in_iter])
# 
# var(N_burnin[1:burn_in_iter])

#



allSamples <- 600
thins <- 3
idx <- (1:(allSamples / thins))*thins
effectiveSize(mcmc(N_iter[idx]))
qplot(idx, N_iter[idx], geom = "line") + geom_hline(aes(yintercept = N_0))
#

ggplot(data = NULL, aes(x = mixtureMeans[,1],
                        y = mixtureMeans[,2])) + geom_point(color = "red", size = 2) +
  geom_point(data = NULL, aes(x = trapsX, y = trapsY))

#

library(microbenchmark)

microbenchmark({
  list_s <- update_N_all(s, beta, CH, N, D, p0, sigma, 
                         Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint = diag(.1, nrow = 2), 
                         lambda_newPoints, probMixture[1], mixtureMeans, mixtureSd, normConst, 
                         S_area, S_usage, trapsX, trapsY, a1, b1, a2, b2)
}, times = 1000)

microbenchmark({
  s <- update_s(s, beta, CH, N, D, p0, sigma, Sigma_prop = diag(.1, nrow = 2),
                Sigma_newpoint = diag(.1, nrow = 2), allTraps, R_traps, S_area, 
                S_usage, trapsX, trapsY, a1, b1, a2, b2)
},{
  list_s <- update_N(s, beta, CH, N, D, p0, sigma, 
                     Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint = diag(.1, nrow = 2), 
                     lambda_newPoints, probMixture[1], mixtureMeans, mixtureSd, normConst, 
                     S_area, S_usage, trapsX, trapsY, a1, b1, a2, b2)
})
