library(Rcpp); library; library(MASS); library(ggplot2); 
library(reshape2)
library(foreach); library(doParallel)
# library(SpatialFunctionsCR);
library(Rcpp); library(RcppArmadillo)
library(here)

#See how many cores we have
# ncl<- detectCores()
# cl <- makeCluster(ncl)
# registerDoParallel(cl)
# stopCluster(cl)

sourceCpp("CR model/code.cpp")
source('CR model/functions.r', echo=TRUE)

# REAL DATA ---------------------------------------------------------------

usingTrueData <- T
load("/cluster/home/osr/ad625/Spatial CR/CR model/JimCorbettData.rda")
load("/cluster/home/osr/ad625/Spatial CR/CR model/leop_poly.rda")
load("/cluster/home/osr/ad625/Spatial CR/CR model/tiger_poly.rda")

leop_polycoord <- leop_polycoord / 1000
tiger_polycoord <- tiger_polycoord / 1000
leop_polyboundaries <- leop_polyboundaries / 1000
tiger_polyboundaries <- tiger_polyboundaries / 1000

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

trapsX <- (trapsX - trapsX_mean) / traps_meansd
trapsY <- (trapsY - trapsY_mean) / traps_meansd

tiger_polycoord[,1] <- (tiger_polycoord[,1] - trapsX_mean) / traps_meansd
tiger_polycoord[,2] <- (tiger_polycoord[,2] - trapsY_mean) / traps_meansd
leop_polycoord[,1] <- (leop_polycoord[,1] - trapsX_mean) / traps_meansd
leop_polycoord[,2] <- (leop_polycoord[,2] - trapsY_mean) / traps_meansd

leop_polyboundaries[c(1,3)] <- (leop_polyboundaries[c(1,3)] - trapsX_mean) / traps_meansd
leop_polyboundaries[c(2,4)] <- (leop_polyboundaries[c(2,4)] - trapsY_mean) / traps_meansd
tiger_polyboundaries[c(1,3)] <- (tiger_polyboundaries[c(1,3)] - trapsX_mean) / traps_meansd
tiger_polyboundaries[c(2,4)] <- (tiger_polyboundaries[c(2,4)] - trapsY_mean) / traps_meansd

tiger_polyboundaries[c(2,3)] <- tiger_polyboundaries[c(3,2)]
leop_polyboundaries[c(2,3)] <- leop_polyboundaries[c(3,2)]

a1_leop <- leop_polyboundaries[1]
b1_leop <- leop_polyboundaries[2]
a2_leop <- leop_polyboundaries[3]
b2_leop <- leop_polyboundaries[4]

a1_tiger <- tiger_polyboundaries[1]
b1_tiger <- tiger_polyboundaries[2]
a2_tiger <- tiger_polyboundaries[3]
b2_tiger <- tiger_polyboundaries[4]

allTraps <- cbind(trapsX, trapsY)
R_traps <- 8 / traps_meansd

# a1 <- min(trapsX) - R_traps
# b1 <- max(trapsX) + R_traps
# a2 <- min(trapsY) - R_traps
# b2 <- max(trapsY) + R_traps

# REAL DATA 2 ---------------------------------------------------------------

usingTrueData <- T
load("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/CR Model/Data/Valeria/valeria_data_magmedio.rda")

# traps$x <- traps$x / 1000
# traps$y <- traps$y / 1000

trapsX <- traps$X.coordinate
trapsY <- traps$Y.coordinate

trapsX_mean <- mean(trapsX)
trapsY_mean <- mean(trapsY)
trapsX_sd <- sd(trapsX)
trapsY_sd <- sd(trapsY)

traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2

trapsX <- (trapsX - mean(trapsX)) / traps_meansd
trapsY <- (trapsY - mean(trapsY)) / traps_meansd

a1 <- min(trapsX) - .05
b1 <- max(trapsX) + .05
a2 <- min(trapsY) - .05
b2 <- max(trapsY) + .05

D1 <- nrow(CH_jaguars)
D2 <- nrow(CH_ocelots)

K <- ncol(CH_jaguars) # number of traps

CH1 <- CH_jaguars
CH2 <- CH_ocelots

# CREATE NEW TRAPS -------------------

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

usingTrueData <- T

# starting params
{
  # theta1 <- exp(-5)
  # theta2 <- exp(-5)
  # theta12 <- exp(-5)
  # beta1 <- 400
  # beta2 <- 400
  theta1 <- exp(-10)
  theta2 <- exp(-10)
  theta12 <- exp(-10)
  beta1 <- 250
  beta2 <- 250
  
  beta1_true <- beta1
  beta2_true <- beta2
  theta1_true <- theta1
  theta2_true <- theta2
  theta12_true <- theta12
  
  p0_1 <- .005#.005
  p0_2 <- .02#.005
  sigma_1 <- 2.5 / traps_meansd
  sigma_2 <- 2.5 / trapsX_sd
  
  p0_1_true <- p0_1
  p0_2_true <- p0_2
  
  sigma_1_true <- sigma_1
  sigma_2_true <- sigma_2
  
}

# data simulation
{
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                        Sigma_newpoint = diag(1, nrow = 2), 
                                        leop_polycoord, leop_polyhole, leop_polystart, 
                                        leop_polyboundaries, tiger_polycoord, 
                                        tiger_polyhole, tiger_polystart, 
                                        tiger_polyboundaries, allTraps, R_traps)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  
  s1 <- s1[1:N1,]
  s2 <- s2[1:N2,]
  
  print(N1)
  print(N2)
  ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") +
    geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black") +
    geom_point(data = NULL, aes(x = trapsX,
                                y = trapsY), color = "green")
    
  
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
  
  captured <- factor(c(rep(1, D1),rep(0, N1 - D1)))
  ggplot() +  geom_point(data = NULL, aes(x = s1_0[,1], y = s1_0[,2], 
                                          color = captured)) + 
    geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "red", size = 2)
  
  captured <- rep(0,N1)
  captured[indexesCaughtIndividuals1] <- 1
  # captured <- factor(c(rep(1, D1),rep(0, N1 - D1)))
  ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2], 
                                          color = factor(captured))) + 
    geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "red", size = 2)
}

qplot(apply(CH_tigers, 1, sum))
qplot(as.vector(CH_tigers))

# CLEAN DATA --------------------------------------------------------------

D1 <- nrow(CH_leopards)
D2 <- nrow(CH_tigers)

K <- ncol(CH_leopards) # number of traps

CH1 <- CH_leopards
CH2 <- CH_tigers

# PLOT THE DATA -----------------------------------------------------------

# 
# {
#   ggplot() +
#     geom_point(data = NULL, aes(x = trapsX, y = trapsY)) + 
#     scale_x_continuous(name = "Kilometers") + 
#     scale_y_continuous(name = "Kilometers") + 
#     stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
#                   fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
#     stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
#                   fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
#     stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
#                   fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
#     stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
#                   fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
#     stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
#                   fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  
#   
# }
# 
# i <- 1
# 
# indexesTraps <- which(CH1[i,] == 1)
# 
# trapsActive <- rep(F, nrow(traps))
# trapsActive[indexesTraps] <- T
# 
# ggplot() +
#   geom_point(data = NULL, aes(x = traps$x, y = traps$y)) + 
#   scale_x_continuous(name = "Kilometers") +  scale_y_continuous(name = "Kilometers") + 
#   scale_color_manual(values=c("black", "red")) +
#   stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
#                 fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
#   stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
#                 fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
#   stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
#                 fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
#   stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
#                 fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
#   stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
#                 fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  
# i <- i+1


grid_x <- seq(a1, b1, length.out = 100)
grid_y <- seq(a2, b2, length.out = 100)
grid_all <- expand.grid(grid_x, grid_y)
pointPresent <- apply(grid_all, 1, function(x){
  checkPointIsInRegionTraps(x, allTraps, R_traps)
})

ggplot() + 
  geom_point(data = NULL, aes(x = grid_all[,1],
                              y = grid_all[,2],
                              color = pointPresent)) +
  geom_point(data = NULL, aes(x = trapsX,
                                       y = trapsY)) 

# PRIOR -------------------------------------------------------------------

Nmax <- 2000

# interaction parameter
{
  # a_theta <- 1
  # b_theta <- 20
  # ggplot(data = NULL) + xlim(c(0,exp(-1))) +
  #   stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))
  
  mu_logtheta <- -30
  sd_logtheta <- 20
  ggplot(data = NULL) + xlim(c(-100,-1)) +
    stat_function(fun = dnorm, args = list(mean = mu_logtheta, sd = sd_logtheta))
}

# proposal interaction parameters
{
  # epsilon_beta <- 70^2
  # epsilon_logtheta1 <- .75
  # epsilon_logtheta2 <- .75
  # epsilon_logtheta12 <- .75
  
  epsilon_beta <- 1000
  epsilon_theta1 <- .0005^2
  epsilon_theta2 <- .0005^2
  epsilon_theta12 <- .0005^2 
}

# adaptation parameters for proposal of interaction parameters
{
  beta_proposal <- .1
  iterAfterAdapting <- 1000
  iterToDiscard <- 250
}

# intensity parameter
{
  meanBeta1 <- 250
  varBeta1 <- 10000
  b_beta1 <- meanBeta1 / varBeta1
  a_beta1 <- meanBeta1 * b_beta1  
  
  ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta1, rate = b_beta1)) + 
    xlim(c(0,400))
  
  meanBeta2 <- 250
  varBeta2 <- 10000
  b_beta2 <- meanBeta2 / varBeta2
  a_beta2 <- meanBeta2 * b_beta2  
  
  ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta2, rate = b_beta2)) + 
    xlim(c(0,400))
}

# proposal value
{
  sigma_prop <- .01
  Sigma_newpoint <- cov(cbind(trapsX, trapsY)) 
}

# capture probability params
{
  a_p <- 1
  b_p <- 1
  
  # leopards (.5 to 2.5)
  # sd_logsigma_0_1 <- 1.5 / trapsX_sd
  sigma_0_1 <- 1.5 / traps_meansd
  sd_sigma_1 <- 2 / traps_meansd
  
  # tigers ( 1.5 to 3.5)
  sigma_0_2 <- 2.5 / traps_meansd
  sd_sigma_2 <- 2 / traps_meansd
  
  # ggplot(data = NULL) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1)) + xlim(c(0.01, .5))
  # ggplot(data = NULL) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1)) + xlim(c(0.01, .5))
  
  # sigma_p0_1_prop <- .004
  # sigma_p0_2_prop <- .004
  sigma_sigma1_prop <- .002
  sigma_sigma2_prop <- .002
}

# points proposal parameters
{
  probMixture <- c(.95,.95)
}

# proposal adaptation
{
  w <- 5 # how much to change lambda
  iterForLambda <- 50 # every how many iter to update lambda
}

# MCMC --------------------------------------------------------------------

pointsToAccept <- 1000
niter <- 10000
nchain <- 2
nthin <- 1

# variables for output
{
  sigma_output <- array(NA, dim = c(nchain, niter, 2))
  s_output <- array(NA, dim = c(nchain, niter, 2, Nmax, 2))
  N_output <- array(NA, dim = c(nchain, niter, 2))
  p0_output <- array(NA, dim = c(nchain, niter, 2))
  params_output <- array(NA, dim = c(nchain, niter, 5))
  
}

# params to update
{
  updateP <- T
  updateS <- T
  updateN <- T
  updateTheta <- T
}

for(chain in 1:nchain) {
  
  # starting values
  {
    
    # home centers
    {
      if(usingTrueData){
        N1 <- 1.25 * D1
        N2 <- 1.25 * D2
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
            s1[i,] <- proposeNewPoint(leop_polyboundaries[1], 
                                      leop_polyboundaries[2], 
                                      leop_polyboundaries[3], 
                                      leop_polyboundaries[4], 
                                      leop_polycoord, 
                                      leop_polyhole,
                                      leop_polystart,
                                      allTraps, R_traps)
          }
        }
        
        s2 <- matrix(NA, nrow = Nmax, ncol = 2)
        for (i in 1:N2) {
          if(i <= D2){ # assign the home center to the average location of the traps where the indiv. was captured
            s2[i,1] <- mean(trapsX[CH2[i,] != 1])
            s2[i,2] <- mean(trapsY[CH2[i,] != 1])
            s2[i,] <- mvrnorm(1, mu = s2[i,], diag(.005, nrow = 2))
          } else { # otherwise a random point
            s2[i,] <- proposeNewPoint(tiger_polyboundaries[1], 
                                      tiger_polyboundaries[2], 
                                      tiger_polyboundaries[3], 
                                      tiger_polyboundaries[4], 
                                      tiger_polycoord, 
                                      tiger_polyhole,
                                      tiger_polystart,
                                      allTraps, R_traps)
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
        
        {
          # c_i <- sapply(1:N1, function(i){
          #   sapply(1:K, function(k){
          #     exp(- (1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
          #   })
          # })
          # 
          # p0_1 <- rbeta(1, sum(CH1) + 1, sum(mean(S_usage) - CH1) + 1) / mean(c_i)
          
          
          c_bar <- sum(CH1)
          
          k_i <- matrix(NA, N1, K)
          
          for (i in 1:N1) {
            for (j in 1:K) {
              k_i[i, j] <- S_usage[j] - CH1[i,j]
            }
          }
          
          c_i <- matrix(NA, N1, K)
          
          for (i in 1:N1) {
            for (j in 1:K) {
              c_i[i, j] <- exp(- (1 / (2 * sigma_1^2)) * ( (s1[i,1] - trapsX[j])^2 + (s1[i,2] - trapsY[j])^2))
            }
          }
          
          c_i <- as.vector(c_i)
          k_i <- as.vector(k_i)
          
          a_fp <- c_bar
          b_fp <-  sum(- c_i * k_i)
          c_fp <- sum(k_i * c_i^2)
          
          p0_1 <- (- b_fp - sqrt(b_fp^2 - 4 * a_fp * c_fp)) / (2 * c_fp)
        }
        
        sigma_2 <- sigma_0_2
        
        {
          # c_i <- sapply(1:N2, function(i){
          #   sapply(1:K, function(k){
          #     exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
          #   })
          # })
          # 
          # p0_2 <- rbeta(1, sum(CH2) + 1, sum(mean(S_usage) - CH2) + 1) / mean(c_i)
          
          c_bar <- sum(CH2)
          
          k_i <- matrix(NA, N2, K)
          
          for (i in 1:N2) {
            for (j in 1:K) {
              k_i[i, j] <- S_usage[j] - CH2[i,j]
            }
          }
          
          c_i <- matrix(NA, N2, K)
          
          for (i in 1:N2) {
            for (j in 1:K) {
              c_i[i, j] <- exp(- (1 / (2 * sigma_2^2)) * ( (s2[i,1] - trapsX[j])^2 + (s2[i,2] - trapsY[j])^2))
            }
          }
          
          c_i <- as.vector(c_i)
          k_i <- as.vector(k_i)
          
          a_fp <- c_bar
          b_fp <-  sum(- c_i * k_i)
          c_fp <- sum(k_i * c_i^2)
          
          p0_2 <- (- b_fp - sqrt(b_fp^2 - 4 * a_fp * c_fp)) / (2 * c_fp)
        }
        
        
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
        # theta1 <- .01
        # theta2 <- .01
        # theta12 <- .01  
        beta1 <- beta1_true
        beta2 <- beta2_true
      }
      
      M <- 50
      
      x_all <- array(NA, dim = c(M, 2, Nmax, 2))
      N_all <- matrix(NA, nrow = M, ncol = 2)
      
      for (m in 1:M) {
        
        print(m)
        
        list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                              niter = 2000,
                                              Sigma_prop = diag(.01, nrow = 2), 
                                              Nmax = Nmax, lambda = beta1,
                                              Sigma_newpoint = diag(1, nrow = 2),
                                              leop_polycoord, 
                                              leop_polyhole,
                                              leop_polystart,
                                              leop_polyboundaries, 
                                              tiger_polycoord, 
                                              tiger_polyhole,
                                              tiger_polystart,
                                              tiger_polyboundaries,
                                              allTraps, R_traps)
        N1_sim <- list_sims$N1
        N2_sim <- list_sims$N2
        x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
        x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
        N_all[m,] <- c(N1_sim, N2_sim)
        
      }
      
      params_values <- matrix(NA, nrow = pointsToAccept * 10 + niter, ncol = 5)
      
    }
    
    # new proposal parameters
    {
      
      numAcceptedPoints <- c(0,0)
      acceptedPoints1 <- matrix(NA, pointsToAccept, 2)
      acceptedPoints2 <- matrix(NA, pointsToAccept, 2)
      N_burnin1 <- rep(NA, pointsToAccept * 10)
      N_burnin2 <- rep(NA, pointsToAccept * 10)
      
      # compute normalizing constant for uniform
      {
        M_proposal <- 10000
        
        normConstUnif_leop <- mean(sapply(1:M_proposal, function(i){
          x <- c(runif(1, leop_polyboundaries[1], leop_polyboundaries[2]), 
                 runif(1, leop_polyboundaries[3], leop_polyboundaries[4]))
          checkPointIsInRegionPolygons(x, leop_polycoord, leop_polyhole, leop_polystart)
        }))
        
        S_area_leop <- (b1_leop - a1_leop) * (b2_leop - a2_leop) * normConstUnif_leop
        
        normConst1 <- c(normConstUnif_leop)
        
        normConstUnif_tiger <- mean(sapply(1:M_proposal, function(i){
          x <- c(runif(1, tiger_polyboundaries[1], tiger_polyboundaries[2]), 
                 runif(1, tiger_polyboundaries[3], tiger_polyboundaries[4]))
          checkPointIsInRegionPolygons(x, tiger_polycoord, tiger_polyhole, tiger_polystart)
        }))
        
        S_area_tiger <- (b1_tiger - a1_tiger) * (b2_tiger - a2_tiger) * normConstUnif_tiger
        
        normConst2 <- c(normConstUnif_tiger)
        
        mixtureMeans1 <- matrix(0, nrow = 0, ncol = 2)
        mixtureSd1 <- matrix(0, nrow = 0, ncol = 2)
        mixtureMeans2 <- matrix(0, nrow = 0, ncol = 2)
        mixtureSd2 <- matrix(0, nrow = 0, ncol = 2)
        
      }
      
    }
    
  }
  
  # output variables
  {
    sigma_iter <- array(NA, dim = c(niter, 2))
    s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
    N_iter <- matrix(NA, nrow = niter, ncol = 2)
    p0_iter <- array(NA, dim = c(niter, 2))
    params_iter <- matrix(NA, nrow = niter, ncol = 5)
    
    # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
    # params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
    
    # acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
    # theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
  }
  
  iter <- 0
  iterInteractionParams <- 1
  burnin_it <- 1
  if(updateN){
    outOfBurnIn <- F
    mixture1Updated <- F
    mixture2Updated <- F
  } else {
    outOfBurnIn <- T
    mixture1Updated <- T
    mixture2Updated <- T
  }
  
  lambda_newPoints1 <- 0
  lambda_newPoints2 <- 0
  mh_ratios1 <- rep(NA, iterForLambda)
  mh_ratios2 <- rep(NA, iterForLambda)
  
  while(iter < niter){
    
    if(!outOfBurnIn){
      print(paste0("Chain = ",chain," / AcceptedPoints1 = ",numAcceptedPoints[1]," / AcceptedPoints2 = ",numAcceptedPoints[2])) 
    # } else if(((iter - nburn)/nthin) %% 5 == 0){
    } else if(iter %% 5 == 0){
      # print(paste0("Chain = ",chain," / Iteration = ",(iter - nburn)/nthin))
      print(paste0("Chain = ",chain," / Iteration = ",iter))
    }
    
    print(paste0("N1 = ",N1," - N2 = ",N2))
    print(paste0("p0_1 = ",p0_1," - p0_2 = ",p0_2))
    # print(paste0("sigma_1 = ",sigma_1," - sigma_2 = ",sigma_2))
    
    # UPDATE P ----------------------------------------------------------------
    
    if(updateP){
      
      # list_p <- update_p(p0_1, sigma_1, p0_2, sigma_2,
      #                    s1, s2, CH1, CH2, trapsX, trapsY,
      #                    S_usage, N1, N2,
      #                    sigma_0_1, sd_sigma_1,
      #                    sigma_0_2, sd_sigma_2,
      #                    sigma_p0_1_prop, sigma_p0_2_prop,
      #                    sigma_sigma1_prop, sigma_sigma2_prop,
      #                    a_p, b_p)
      # p0_1 <- list_p$p0_1
      # p0_2 <- list_p$p0_2
      # sigma_1 <- list_p$sigma_1
      # sigma_2 <- list_p$sigma_2
      
      list_p <- update_psigma_cpp(p0_1, sigma_1, CH1, s1, N1, 
                                  K, S_usage, trapsX, trapsY,
                                  sigma_0_1, sd_sigma_1, sigma_sigma1_prop,
                                  a_p, b_p)
      p0_1 <- list_p[1]
      sigma_1 <- list_p[2]
      
      list_p2 <- update_psigma_cpp(p0_2, sigma_2, CH2, s2, N2, 
                                   K, S_usage, trapsX, trapsY,
                                   sigma_0_2, sd_sigma_2, sigma_sigma2_prop,
                                   a_p, b_p)
      p0_2 <- list_p2[1]
      sigma_2 <- list_p2[2]
      
    }
    
    # UPDATE HOME CENTERS -----------------------------------------------
    
    # print("Sample home centers")
    
    if(updateS){
      
      list_s <- update_s_cpp(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
                             CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                             Sigma_newpoint, 
                             lambda_movedPoints1 = 20, 
                             lambda_movedPoints2 = 20,
                             lambda_newPoints1 = 0, 
                             lambda_newPoints2 = 0, 
                             probMixture[1],
                             probMixture[2],
                             mixtureMeans1, mixtureMeans2, 
                             mixtureSd1, mixtureSd2,
                             normConst1, normConst2, 
                             leop_polycoord, 
                             leop_polyhole,
                             leop_polystart,
                             leop_polyboundaries, 
                             tiger_polycoord, 
                             tiger_polyhole,
                             tiger_polystart,
                             tiger_polyboundaries,
                             allTraps, R_traps,
                             S_area_leop, S_area_tiger, 
                             S_usage, trapsX, trapsY)
      
      s1 <- list_s$s1
      s2 <- list_s$s2
      N1 <- list_s$N1
      N2 <- list_s$N2
      
    }
    
    # UPDATE N ------------------
    
    # print("Sample N")
    
    if(updateN){
      
      if(outOfBurnIn){
        
        list_s <- update_N_cpp(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
                                   CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                                   Sigma_newpoint, lambda_movedPoints1 = 20, lambda_movedPoints2 = 20,
                                   lambda_newPoints1, lambda_newPoints2, probMixture[1],
                                   probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
                                   normConst1, normConst2, 
                                   leop_polycoord, 
                                   leop_polyhole,
                                   leop_polystart,
                                   leop_polyboundaries, 
                                   tiger_polycoord, 
                                   tiger_polyhole,
                                   tiger_polystart,
                                   tiger_polyboundaries,
                               allTraps, R_traps,
                                   S_area_leop, S_area_tiger, 
                                   S_usage, trapsX,
                                   trapsY)
        
      } else {
        list_s <-  update_s_cpp_withacceptance(s1, s2,  theta1, theta2, theta12, beta1, beta2, CH1,
                                               CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                               S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
                                               leop_polycoord, 
                                               leop_polyhole,
                                               leop_polystart,
                                               leop_polyboundaries, 
                                               tiger_polycoord, 
                                               tiger_polyhole,
                                               tiger_polystart,
                                               tiger_polyboundaries,
                                               allTraps, R_traps)
      }
      
      s1 <- list_s$s1
      s2 <- list_s$s2
      N1 <- list_s$N1
      N2 <- list_s$N2
      mh_ratio1 <- list_s$mh_ratio1
      mh_ratio2 <- list_s$mh_ratio2
      
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
      
      if(outOfBurnIn){
        mh_ratios1[((iter - 1) %% iterForLambda) + 1] <- mh_ratio1
        mh_ratios2[((iter - 1) %% iterForLambda) + 1] <- mh_ratio2
      }
      
    }
    
    # UPDATE LAMBDA ----------
    
    # print("Sample Lambda")
    
    if(outOfBurnIn){
      if(iter %% iterForLambda == 0){
        mh_ratios1[mh_ratios1 > 1] <- 1
        meanMHratios <- mean(mh_ratios1[mh_ratios1 > 0])
        lambda_newPoints1 <- lambda_newPoints1 + (meanMHratios - .2) * w * (1 / floor(iter / iterForLambda))
        mh_ratios1 <- rep(NA, iterForLambda)
        
        mh_ratios2[mh_ratios2 > 1] <- 1
        meanMHratios <- mean(mh_ratios2[mh_ratios2 > 0])
        lambda_newPoints2 <- lambda_newPoints2 + (meanMHratios - .2) * w * (1 / floor(iter / iterForLambda))
        mh_ratios2 <- rep(NA, iterForLambda)
      }  
    }
    
    # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------

    if(updateTheta){
      
      if(iterInteractionParams > iterAfterAdapting){
        Sigma_n <- cov(params_values[iterToDiscard:(iterInteractionParams-1),c(1,3,5)])
        Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
          beta_proposal * diag(c(epsilon_beta, epsilon_theta1, epsilon_theta12), nrow = 3) 
      } else {
        Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta1, epsilon_theta12), nrow = 3) 
      }
      
      proposed_theta <- mvrnorm(1, c(beta1, theta1, theta12),
                                Sigma = Sigma_proposal)
      beta1_star <- proposed_theta[1]
      theta1_star <- proposed_theta[2]
      theta12_star <- proposed_theta[3]
      
      if(theta1_star > 0 & theta12_star > 0){
        
        ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                    x_all[,1,,], 
                                    x_all[,2,,],
                                    c(beta1, beta2,
                                      theta1, theta2, theta12),
                                    c(beta1_star, beta2,
                                      theta1_star, theta2, theta12_star),
                                    N_all,
                                    mu_logtheta, sd_logtheta,
                                    a_beta1, b_beta1)
        
        if(!is.na(ratio) & theta1_star > 0){
          
          if(runif(1) < ratio){
            beta1 <- beta1_star
            theta1 <- theta1_star
            theta12 <- theta12_star
            
            # list_xall <- update_x_all_foreach(x_all, N_all, 100,
            #                                   theta1, theta2, theta12,
            #                                   beta1, beta2, diag(0.05, nrow = 2),
            #                                   Sigma_newpoint = diag(1, nrow = 2),
            #                                   allTraps, R_traps,
            #                                   a1, b1, a2, b2)
            # x_all <- list_xall$x_all
            # N_all <- list_xall$N_all
            
            x_all1 <- x_all[,1,,]
            x_all2 <- x_all[,2,,]
            
            list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2),
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          leop_polycoord, 
                                          leop_polyhole,
                                          leop_polystart,
                                          leop_polyboundaries, 
                                          tiger_polycoord, 
                                          tiger_polyhole,
                                          tiger_polystart,
                                          tiger_polyboundaries,
                                          allTraps, R_traps)
            x_all[,1,,] <-  list_xall$x_all1
            x_all[,2,,] <-  list_xall$x_all2
            N_all <- list_xall$N_all
            
          }
          
        }
        
      }
      
      if(iterInteractionParams <= nrow(params_values)){
        params_values[iterInteractionParams,c(1,3,5)] <- c(beta1, theta1, theta12)
      } else {
        params_values <- rbind(params_values,
                               matrix(NA, nrow = pointsToAccept * 10 + niter, ncol = 5))
        params_values[iterInteractionParams,c(1,3,5)] <- c(beta1, theta1, theta12)
      }
      
      
      
    } else {
      
      beta1 <- rgamma(1, a_beta1 + N1, b_beta1 + 1)
      
    }
    
    # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
    
    if(updateTheta){
      
      if(iter > iterAfterAdapting){
        Sigma_n <- cov(params_values[iterToDiscard:(iterInteractionParams-1),c(2,4,5)])
        Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
          beta_proposal * diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
      } else {
        Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
      }
      
      proposed_theta <- mvrnorm(1, c(beta2, theta2, theta12),
                                Sigma = Sigma_proposal)
      beta2_star <- proposed_theta[1]
      theta2_star <- proposed_theta[2]
      theta12_star <- proposed_theta[3]
      
      if(theta2_star > 0 & theta12_star > 0){
        
        # ratio <- hastings_ratio_2(s1, s2, N1, N2,
        #                           x_all[,1,,], 
        #                           x_all[,2,,],
        #                           N_all,
        #                         beta1, beta2,
        #                         beta1, beta2_star,
        #                         theta1, theta2, theta12,
        #                         theta1, theta2_star, theta12_star,
        #                         mu_logtheta, sd_logtheta,
        #                         # a_theta, b_theta,
        #                         a_beta, b_beta)
        ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                    x_all[,1,,], 
                                    x_all[,2,,],
                                    c(beta1, beta2,
                                      theta1, theta2, theta12),
                                    c(beta1, beta2_star,
                                      theta1, theta2_star, theta12_star),
                                    N_all,
                                    mu_logtheta, sd_logtheta,
                                    # a_theta, b_theta,
                                    a_beta2, b_beta2)
        
        # print(paste0("mh_ratio_beta2 = ",ratio))
        if(!is.na(ratio) & theta2_star > 0){
          
          if(runif(1) < ratio){
            beta2 <- beta2_star
            theta2 <- theta2_star
            theta12 <- theta12_star
            
            x_all1 <- x_all[,1,,]
            x_all2 <- x_all[,2,,]
            
            list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2),
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          leop_polycoord, 
                                          leop_polyhole,
                                          leop_polystart,
                                          leop_polyboundaries, 
                                          tiger_polycoord, 
                                          tiger_polyhole,
                                          tiger_polystart,
                                          tiger_polyboundaries,
                                          allTraps, R_traps)
            x_all[,1,,] <-  list_xall$x_all1
            x_all[,2,,] <-  list_xall$x_all2
            N_all <- list_xall$N_all
            
            # list_xall <- update_x_all_foreach(x_all, N_all, 100,
            #                                   theta1, theta2, theta12,
            #                                   beta1, beta2, diag(0.05, nrow = 2),
            #                                   Sigma_newpoint = diag(1, nrow = 2),
            #                                   allTraps, R_traps,
            #                                   a1, b1, a2, b2)
            # x_all <- list_xall$x_all
            # N_all <- list_xall$N_all
            
            # print("accepted 2")
            # acceptances_theta[iter,1] <- 1
          }
          
        }
        
      }
      
      params_values[iterInteractionParams,c(2,4,5)] <- c(beta2, theta2, theta12)  
      
      iterInteractionParams <- iterInteractionParams + 1
      
    } else {
      
      beta2 <- rgamma(1, a_beta2 + N2, b_beta2 + 1)
      
    }
    
    # WRITE RESULTS ---------------
    
    if(updateN) {
      outOfBurnIn <- numAcceptedPoints[1] >= pointsToAccept & numAcceptedPoints[2] >= pointsToAccept
      
    }
    
    if(!outOfBurnIn){
      N_burnin1[burnin_it] <- N1
      N_burnin2[burnin_it] <- N2
      burnin_it <- burnin_it + 1  
    }
    
    # update proposals
    {
      if(outOfBurnIn & !mixture1Updated){
        
        # update lambda 
        
        burn_in_iter <- sum(!is.na(N_burnin1))
        
        lambda_newPoints1 <- (sd(N_burnin1[iterAfterAdapting:burn_in_iter])) / 4
        
        # update q_b
        
        numCenters_all <- 5:20
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
          
          normConstNormals1[i] <- mean(sapply(1:M_proposal, function(j){
            x <- c(rnorm(1, mixtureMeans1[i,1], mixtureSd1[i,1]), 
                   rnorm(1, mixtureMeans1[i,2], mixtureSd1[i,2]))
            checkPointIsInRegionPolygons(x, leop_polycoord, leop_polyhole, leop_polystart)
          }))
        }
        
        normConst1 <- c(normConstUnif_leop, normConstNormals1)
        
        mixture1Updated <- T
      }
      
      if(outOfBurnIn & !mixture2Updated){
        
        # update lambda 
        
        burn_in_iter <- sum(!is.na(N_burnin2))
        
        lambda_newPoints2 <- (sd(N_burnin2[iterAfterAdapting:burn_in_iter])) / 4
        
        # update q_b
        
        numCenters_all <- 5:20
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
          
          normConstNormals2[i] <- mean(sapply(1:M_proposal, function(j){
            x <- c(rnorm(1, mixtureMeans2[i,1], mixtureSd2[i,1]), 
                   rnorm(1, mixtureMeans2[i,2], mixtureSd2[i,2]))
            checkPointIsInRegionPolygons(x, tiger_polycoord, tiger_polyhole, tiger_polystart)
          }))
        }
        
        normConst2 <- c(normConstUnif_tiger, normConstNormals2)
        
        mixture2Updated <- T
      }    
    }
    
    if(outOfBurnIn){
      
      iter <- iter + 1
      
      # trueIter <- (iter - nburn)/nthin
      trueIter <- iter
      
      s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
      s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
      p0_iter[trueIter,] <- c(p0_1, p0_2)
      sigma_iter[trueIter,] <- c(sigma_1, sigma_2)
      N_iter[trueIter,] <- c(N1,N2)  
      params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
      # papangelou_density_iter[trueIter,,,] <- computeNewPointsDensity(gridLength_x, gridLength_y,
      # s1, s2, N1, N2,
      # theta1,  theta2,  theta12,
      # beta1,  beta2)
      
    }
    
    if(iter %% 250 == 0){
      # plotCurrent1 <- qplot(1:iter, params_iter[1:iter,1], geom =  "line")
      # plotCurrent2 <- qplot(1:iter, log(params_iter[1:iter,3]), geom =  "line")
      # plotCurrent3 <- qplot(1:iter, log(params_iter[1:iter,5]), geom =  "line")
      # ggsave(filename = paste0("plot_param1_iter",iter,".jpg"), plotCurrent1)
      # ggsave(filename = paste0("plot_param2_iter",iter,".jpg"), plotCurrent2)
      # ggsave(filename = paste0("plot_param3_iter",iter,".jpg"), plotCurrent3)
    }
    
  }
  
  # Write results in MCMC output
  {
    sigma_output[chain,,] <- sigma_iter
    s_output[chain,,,,] <- s_iter#array(NA, dim = c(nchain, niter, 2, Nmax, 2))
    N_output[chain,,] <- N_iter#array(NA, dim = c(nchain, niter, 2))
    p0_output[chain,,] <- p0_iter#array(NA, dim = c(nchain, niter, 2))
    params_output[chain,,] <- params_iter#array(NA, dim = c(nchain, niter, 5))
  }
  
}

# DIAGNOSTICS MIXTURE PROPOSAL --------------

ggplot(data = NULL, aes(x = trapsX, y = trapsY)) + geom_point() +
  geom_point(data = NULL, aes(x = mixtureMeans1[,1], y = mixtureMeans1[,2]), size = 2, color = "red")

# SIMULATED DIAGNOSTICS ---------------------------------------------------

qplot(1:niter, p0_iter[,1]) + geom_hline(data = NULL, aes(yintercept = p0_1_true))
qplot(1:niter, p0_iter[,2]) + geom_hline(data = NULL, aes(yintercept = p0_2_true))
qplot(1:niter, sigma_iter[,1]) + geom_hline(data = NULL, aes(yintercept = sigma_1_true))
qplot(1:niter, sigma_iter[,2]) + geom_hline(data = NULL, aes(yintercept = sigma_2_true))

qplot(1:niter, N_iter[,1]) + geom_hline(data = NULL, aes(yintercept = N1_0)) +
  geom_hline(data = NULL, aes(yintercept = D1), color = "red")
qplot(1:niter, N_iter[,2]) + geom_hline(data = NULL, aes(yintercept = N2_0)) +
  geom_hline(data = NULL, aes(yintercept = D2), color = "red")

qplot(1:niter, params_iter[,1]) + geom_hline(data = NULL, aes(yintercept = beta1_true))
qplot(1:niter, params_iter[,2]) + geom_hline(data = NULL, aes(yintercept = beta2_true))
qplot(1:niter, log(params_iter[,3])) + geom_hline(data = NULL, aes(yintercept = log(theta1_true)))
qplot(1:niter, log(params_iter[,4])) + geom_hline(data = NULL, aes(yintercept = log(theta2_true)))
qplot(1:niter, log(params_iter[,5])) + geom_hline(data = NULL, aes(yintercept = log(theta12_true)))

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

p0_output <- p0_output[,-(1:5000),]
sigma_output <- sigma_output[,-(1:5000),]
N_output <- N_output[,-(1:5000),]
params_output <- params_output[,-(1:5000),]

# diagnostics p
{
  plotVar(p0_output[,,1])
  plotVar(p0_output[,,2])
  plotVar(N_output[,,1])
  plotVar(N_output[,,2])
  plotVar(sigma_output[,,1])
  plotVar(sigma_output[,,2])
  plotVar(log(sigma_output[,,2]))
  plotVar(sigma_output[,,1] * log(p0_output[,,1]))
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
  scale_y_continuous(name = "", limits = c(0, 0.02)) +
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

qplot(1:niter, log(params_iter[,3] / params_iter[,5])) #+ ylim(c(0, 2))


# DETERMINE AREA SIZE --------

M_proposal <- 300000
ratio_area_leop <- mean(sapply(1:M_proposal, function(j){
  x <- c(runif(1, a1_leop, b1_leop), 
         runif(1, a2_leop, b2_leop))
  checkPointIsInRegionPolygonsAndTraps(x, leop_polycoord, leop_polyhole, 
                                       leop_polystart, allTraps, R_traps)
}))

totalArea_leop <- ratio_area_leop * (b1_leop - a1_leop) * traps_meansd * 
  (b2_leop - a2_leop) * traps_meansd 

ratio_area_tiger <- mean(sapply(1:M_proposal, function(j){
  x <- c(runif(1, a1_tiger, b1_tiger), 
         runif(1, a2_tiger, b2_tiger))
  checkPointIsInRegionPolygonsAndTraps(x, tiger_polycoord, tiger_polyhole, tiger_polystart,
                               allTraps, R_traps)
}))

totalArea_tiger <- ratio_area_tiger * (b1_tiger - a1_tiger) * traps_meansd * 
  (b2_tiger - a2_tiger) * traps_meansd 

# PLOTS -------------------------------------------------------------------

# plots for gamma
{
  params_output_all <- apply(params_output, 3, c)
  
  CI_params <- apply(params_output_all[,3:5], 2, function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  # CI_params <- apply(params_iter[,3:5], 2, function(x){
  # quantile(x, probs = c(.05,.5,.95))
  # })
  
  plotgamma <- ggplot(data = NULL, aes(x = reorder(factor(c("theta1","theta2","theta12")),c(1,3,2)), 
                                       y = CI_params[2,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[3,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c(expression(theta[1]),expression(theta[1][2]),expression(theta[2])),
                     name = "Interaction parameters") +
    scale_y_continuous(name = "") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), 
          panel.grid = element_line(colour = "grey92"), 
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

# plots for gamma log ratios
{
  params_output_all <- apply(params_output, 3, c)
  
  CI_params_theta1_12 <- quantile(log(params_output_all[,3] / params_output_all[,5]), probs = c(.025,.5,.975))
  CI_params_theta2_12 <- quantile(log(params_output_all[,4] / params_output_all[,5]), probs = c(.025,.5,.975))
  # CI_params_theta1_12 <- quantile(log(params_output_all[,5]) - log(params_output_all[,3]), probs = c(.025,.5,.975))
  # CI_params_theta2_12 <- quantile(log(params_output_all[,5]) - log(params_output_all[,4]), probs = c(.025,.5,.975))
  
  CI_params <- cbind(CI_params_theta1_12, CI_params_theta2_12)
  # CI_params <- apply(params_iter[,3:5], 2, function(x){
  # quantile(x, probs = c(.05,.5,.95))
  # })
  
  plotgamma <- ggplot(data = NULL, aes(x = reorder(factor(c("theta12 - theta1","theta12 - theta2")),c(1,2)),
                                       y = CI_params[2,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[3,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    geom_hline(aes(yintercept = 0), color = "red",size = 1) +
    scale_x_discrete(labels = c(expression(log(gamma[1] / gamma[1][2])),
                                expression(log(gamma[2] / gamma[1][2]))),
                     name = "") +
    scale_y_continuous(name = "Value") +
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
    # theme(panel.background = element_rect(fill = "white", 
    #                                       colour = NA), 
    #       panel.border = element_rect(fill = NA, 
    #                                   colour = "grey20"), 
    #       panel.grid = element_line(colour = "grey92"), 
    #       panel.grid.minor = element_line(size = rel(0.5)), 
    #       strip.background = element_rect(fill = "grey85", 
    #                                       colour = "grey20"), 
    #       legend.key = element_rect(fill = "white", 
    #                                 colour = NA),
    #       plot.title = element_text(hjust = 0.5, size = 20),
    #       axis.title = element_text(size = 13, face = "bold"),
    #       axis.text = element_text(size = 13, face = "bold"))
  plotgamma
  setwd("/cluster/home/osr/ad625/Spatial CR/CR model/ImgToExport")
  
  ggsave(filename = "gamma_ratios.jpeg",  plotgamma,
         width = 6.36, height = 5.25, units = "in")
}

# plots for N
{
  N_output_all <- apply(N_output, 3, c)
  
  plotN <- ggplot(data = NULL, aes(x = N_output_all[,1])) + 
    geom_bar(aes(y = (..count..)/sum(..count..)), 
             fill = "cornsilk", color = "black", stat = "count", width = .3, size = .5) +
    scale_x_continuous(breaks = sort(unique(N_iter[,1])), 
                       labels =  round(sort(unique(N_iter[,1])) /  totalArea_leop,3),
                       minor_breaks = sort(unique(N_iter[,1])),   
                       name = expression(paste("Number of leopards per ",km^2))) +
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
                       labels =  round(sort(unique(N_iter[,2])) /  totalArea,3),
                       minor_breaks = sort(unique(N_iter[,2])),   
                       name = expression(paste("Number of tigers per ",km^2))) +
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
  setwd("/cluster/home/osr/ad625/Spatial CR/CR model/ImgToExport")
  
  ggsave(filename = "Popsize2.jpeg",  plotN2,
         width = 6.36, height = 5.25, units = "in")
  
  
}

# boxplots for N
{
  N_output_all <- apply(N_output / 1, 3, c)
  
  i <- 0
  CI_params <- apply(N_output_all, 2, function(x){
    
    i <<- i + 1 
    if(i == 1){
      c(quantile(x / totalArea_leop, probs = c(.025,.975)), mean(x / totalArea_leop))  
    } else {
      c(quantile(x / totalArea_tiger, probs = c(.025,.975)), mean(x / totalArea_tiger))  
    }
    
  })
  
  # CI_params <- apply(params_iter[,3:5], 2, function(x){
  # quantile(x, probs = c(.05,.5,.95))
  # })
  
  plot_N <- ggplot(data = NULL, aes(x = reorder(factor(c("Leopards","Tigers")),c(1,2)), y = CI_params[3,], 
                                    ymin = CI_params[1,], 
                                    ymax = CI_params[2,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c("Leopards","Tigers"),
                     name = "Species") +
    scale_y_continuous(name = expression(paste("Number of individuals per ",km^2))) + 
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
  plot_N
  setwd("/cluster/home/osr/ad625/Spatial CR/CR model/ImgToExport")
  
  ggsave(filename = "boxplot_N.jpeg",  plot_N,
         width = 6.36, height = 5.25, units = "in")
}

# histograms for N
{
  plotN <- ggplot(data = NULL, aes(x = N_output_all[,1])) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", alpha = .6, binwidth = 3) +
    geom_vline(aes(xintercept = mean(N_output_all[,1])), color = "red", size = 2) +
    scale_x_continuous(breaks = seq(120,200,by = 10), 
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Number of leopards") +
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
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotN
  
  ggsave(filename = "N_leopards.jpeg",  plotN,
         width = 6.36, height = 5.25, units = "in")
  
  plotN2 <- ggplot(data = NULL, aes(x = N_output_all[,2])) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", binwidth = 3) +
    geom_vline(aes(xintercept = mean(N_output_all[,2])), color = "red", size = 2) +
    scale_x_continuous(breaks = seq(230,320,by = 10), 
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Number of tigers") +
    scale_y_continuous(name = "Density") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), 
          panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotN2
  
  ggsave(filename = "N_tigers.jpeg",  plotN2,
         width = 6.36, height = 5.25, units = "in")
}

# histograms for p
{
  p_output_all <- apply(p0_output, 3, c)
  
  plotp <- ggplot(data = NULL, aes(x = p_output_all[,1])) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", binwidth = .00015) +
    geom_vline(aes(xintercept = mean(p_output_all[,1])), color = "red", size = 2) +
    coord_cartesian(xlim = c(.010,.018)) + 
    scale_x_continuous(breaks = seq(.010,.018,by = .002),
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Baseline capture probability of leopards") +
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
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotp
  
  ggsave(filename = "p_leopards.jpeg",  plotp,
         width = 6.36, height = 5.25, units = "in")
  
  plotp2 <- ggplot(data = NULL, aes(x = p_output_all[,2])) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", binwidth = .0001) +
    geom_vline(aes(xintercept = mean(p_output_all[,2])), color = "red", size = 2) +
    coord_cartesian(xlim = c(.025,.029)) + 
    scale_x_continuous(breaks = seq(.02,.029,by = .001),
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Baseline capture probability of tigers") +
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
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotp2
  
  ggsave(filename = "p_tigers.jpeg",  plotp2,
         width = 6.36, height = 5.25, units = "in")
}

# boxplots for p
{
  p_output_all <- apply(p0_output, 3, c)
  
  CI_params <- apply(p_output_all, 2, function(x){
    c(quantile(x, probs = c(.025,.975)), mean(x))
  })
  
  # CI_params <- apply(params_iter[,3:5], 2, function(x){
  # quantile(x, probs = c(.05,.5,.95))
  # })
  
  plot_p0 <- ggplot(data = NULL, aes(x = reorder(factor(c("Leopards","Tigers")),c(1,2)), y = CI_params[3,], 
                                     ymin = CI_params[1,], 
                                     ymax = CI_params[2,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c("Leopards","Tigers"),
                     name = "Species") +
    scale_y_continuous(name = expression(p[0])) + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), 
          panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plot_p0
  setwd("/cluster/home/osr/ad625/Spatial CR/CR model/ImgToExport")
  
  ggsave(filename = "boxplot_p0.jpeg",  plot_p0,
         width = 6.36, height = 5.25, units = "in")
}

# boxplots for sigma
{
  sigma_output_all <- apply(sigma_output, 3, c)
  
  CI_params <- apply(sigma_output_all, 2, function(x){
    c(quantile(x, probs = c(.025,.975)), mean(x))
  })
  
  CI_params[,1] <- CI_params[,1] * traps_meansd
  CI_params[,2] <- CI_params[,2] * traps_meansd
  
  # CI_params <- apply(params_iter[,3:5], 2, function(x){
  # quantile(x, probs = c(.05,.5,.95))
  # })
  
  plot_sigma <- ggplot(data = NULL, aes(x = reorder(factor(c("Leopards","Tigers")),c(1,2)), y = CI_params[3,], 
                                        ymin = CI_params[1,], 
                                        ymax = CI_params[2,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c("Leopards","Tigers"),
                     name = "Species") +
    scale_y_continuous(name = expression(sigma)) + 
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
  plot_sigma
  setwd("/cluster/home/osr/ad625/Spatial CR/CR model/ImgToExport")
  
  ggsave(filename = "boxplot_sigma.jpeg",  plot_sigma,
         width = 6.36, height = 5.25, units = "in")
}

# histograms for sigma
{
  sigma_output_all <- apply(sigma_output, 3, c)
  
  plotsigma <- ggplot(data = NULL, aes(x = sigma_output_all[,1] * trapsX_sd)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", binwidth = .02) +
    geom_vline(aes(xintercept = mean(sigma_output_all[,1] * trapsX_sd)), color = "red", size = 2) +
    coord_cartesian(xlim = c(1.5, 2.1)) +
    scale_x_continuous(breaks = seq(1.5,2.1,by = .1), 
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Sigma estimate of leopards") +
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
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotsigma
  
  ggsave(filename = "sigma_leopards.jpeg",  plotsigma,
         width = 6.36, height = 5.25, units = "in")
  
  plotsigma2 <- ggplot(data = NULL, aes(x = sigma_output_all[,2] * trapsY_sd)) + 
    geom_histogram(aes(y = (..count..)/sum(..count..)), 
                   fill = "cornsilk", color = "black", binwidth = .004) +
    geom_vline(aes(xintercept = mean(sigma_output_all[,2] * trapsY_sd)), color = "red", size = 2) +
    scale_x_continuous(breaks = seq(1.6,2,by = .025), 
                       # labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       # labels =  sort(unique(N_iter[,1])),
                       # minor_breaks = sort(unique(N_iter[,1])),   
                       name = "Sigma estimate of tigers") +
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
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  plotsigma2
  
  ggsave(filename = "sigma_tigers.jpeg",  plotsigma2,
         width = 6.36, height = 5.25, units = "in")
  
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


# HOME CENTERS2 -----

grid_x <- seq(a1, b1, length.out = 100)
grid_y <- seq(a2, b2, length.out = 100)
grid_all <- expand.grid(grid_x, grid_y)
pointPresent <- apply(grid_all, 1, function(x){
  checkPointIsInRegionTraps(x, allTraps, R_traps)
})

ggplot() + 
  geom_point(data = NULL, aes(x = grid_all[,1],
                              y = grid_all[,2],
                              color = pointPresent)) +
  geom_point(data = NULL, aes(x = trapsX,
                              y = trapsY)) + 
  geom_point(data = NULL, aes(x = s2[1:D2,1],
                              y = s2[1:D2,2]), color = "green") + 
  geom_point(data = NULL, aes(x = s2[(D2 + 1):N2,1],
                              y = s2[(D2 + 1):N2,2]), color = "red") 

# HOME CENTERS ---

grid_x <- seq(a1_tiger, b1_tiger, length.out = 200)
grid_y <- seq(a2_tiger, b2_tiger, length.out = 200)

grid_all <- expand.grid(grid_x, grid_y)
pointInside <- rep(NA, nrow(grid_all))
for (i in 1:nrow(grid_all)) {
  point <- as.numeric(grid_all[i,])
  pointInside[i] <- checkPointIsInRegionPolygons(point, tiger_polycoord, 
                                         tiger_polyhole, tiger_polystart)
}

library(ggplot2)

ggplot() + 
  geom_point(data = NULL, aes(x = grid_all$Var1,
                              y = grid_all$Var2,
                              color = pointInside)) +
  scale_color_manual(values = c("grey","white")) +
  theme_bw() + 
  geom_point(data = NULL, aes(x = trapsX,
                              y = trapsY), color = "black") + 
  geom_point(data = NULL, aes(x = s2[1:D2,1],
                              y = s2[1:D2,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2[(D2 + 1):N2,1],
                              y = s2[(D2 + 1):N2,2]), color = "green")
