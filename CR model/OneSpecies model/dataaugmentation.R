library(MASS); library(reshape2); library(ggplot2)
library(foreach); library(doParallel)
# library(SpatialFunctionsCR);
library(Rcpp); library(RcppArmadillo)

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
# stopCluster(cl)

setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/CR Model/Simplified model")
sourceCpp("code.cpp")
source('functions.r', echo=TRUE)

# CREATE NEW TRAPS -------

a1 <- 0
b1 <- 1
a2 <- 0
b2 <- 1

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


traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2

trapsX <- (trapsX - mean(trapsX)) / traps_meansd
trapsY <- (trapsY - mean(trapsY)) / traps_meansd

allTraps <- cbind(trapsX, trapsY)
R_traps <- .75

a1 <- min(trapsX) - R_traps
b1 <- max(trapsX) + R_traps
a2 <- min(trapsY) - R_traps
b2 <- max(trapsY) + R_traps

K <- length(trapsX)
S_usage <- rep(1, K)

# SIMULATED DATA ----------------------------------------------------------

usingTrueData <- F

# starting params
{
  beta <- 200
  
  p0 <- .005
  sigma <- 3 / trapsX_sd
  
  beta_true <- beta
  p0_true <- p0
  sigma_true <- sigma
}

# data simulation
{
  N <- rpois(1, beta)
  s <- t(sapply(1:N, function(i){
    c(runif(1, a1, b1), runif(1, a2, b2))
  }))
  
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

# PRIOR -------------------------------------------------------------------

# intensity parameter
{
  M <- 5 * N_0
  
  phi <- beta_true / M

}

# proposal value
sigma_prop <- .01

# capture probability params
{
  a_p <- 1
  b_p <- 1
  
  sigma_0 <- sigma_true
  sd_sigma <- 1
  
  sigma_p0_prop <- .0002
  sigma_sigma_prop <- .005
}

# area
{
  S_area <- (b2 - a2) * (b1 - a1)
}


# MCMC --------------------------------------------------------------------

nburn <- 500
niter <- 1000
nchain <- 1
nthin <- 1

# home centers
{
  if(usingTrueData){
    
  } else {
    N <- N_0
    z <- rep(0, M)
    z[1:N] <- 1
  }
  
  CH <- rbind(CH[1:D,],
              matrix(0, nrow = N - D, ncol = K))
  
  if(usingTrueData){
    
    s <- matrix(NA, nrow = M, ncol = 2)
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
    
    s <- matrix(NA, nrow = M, ncol = 2)
    
    s[1:N,] <- s_0
    s[(N+1):M,] <- cbind(runif(M - N, a1, b1),
                         runif(M - N, a2, b2))
    
  }
  
}

# capture probability parameter
{
  if(usingTrueData){
    
    
  } else {
    p0 <- p0_true
    sigma <- sigma_true
  }
  
}

# point process parameters
{
  if(usingTrueData){
    beta <- N
  } else {
    # beta <- beta_true
    theta <- theta_true
  }
  
}

# output variables
{
  sigma_iter <- rep(NA, niter)
  # s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
  N_iter <- rep(NA, niter)
  p0_iter <- rep(NA, niter)
}

for(iter in seq_len(niter + nburn)){
  
  if(iter == nburn + 1){
    st1 <- Sys.time()
  }
  
  if(iter < nburn){
    print(paste0("Burn In Iteration = ",iter)) 
  } else if(((iter - nburn)/nthin) %% 5 == 0){
    print(paste0("Iteration = ",(iter - nburn)/nthin))
  }
  
  print(paste0("N = ",N))
  
  # UPDATE P ----------------------------------------------------------------
  
  # list_p <- update_p(p0, sigma, s, CH, trapsX, trapsY, 
  #                    S_usage, N, sigma_0, sd_sigma_1,
  #                    sigma_p0_prop, sigma_sigma_prop, a_p, b_p)
  # p0 <- list_p$p0
  # sigma <- list_p$sigma
  
  # UPDATE HOME CENTERS -----------------------------------------------
  
  if(iter < nburn){
    s <- update_s_da(s, beta, theta, z, CH, N, D, p0, sigma, Sigma_prop = diag(.1, nrow = 2),
                     Sigma_newpoint = diag(.1, nrow = 2), allTraps, R_traps, S_area, 
                     S_usage, trapsX, trapsY, a1, b1, a2, b2)  
  } else {
    s <- update_s_da(s, beta, theta, z, CH, N, D, p0, sigma, Sigma_prop = diag(.1, nrow = 2),
                     Sigma_newpoint = diag(.1, nrow = 2), allTraps, R_traps, S_area, 
                     S_usage, trapsX, trapsY, a1, b1, a2, b2)  

  }
  
  # UPDATE N -----------------------------------------------
  
  list_z <- update_z_cpp(z, s, N, D, phi, theta, K, S_usage, trapsX, trapsY, p0, sigma)
  # list_z <- update_z(z, s, beta, N, phi)
  z <- list_z$z
  N <- list_z$N
  
  CH <- rbind(CH[1:D,],
              matrix(0, nrow = N - D, ncol = K))
  
  phi <- rbeta(1,sum(z)+ 1,sum(1-z) + 1)
  
  # WRITE RESULTS ---------------
  
  if(iter > nburn){
    
    trueIter <- (iter - nburn)/nthin
    # trueIter <- iter
    
    p0_iter[trueIter] <- p0
    sigma_iter[trueIter] <- sigma
    N_iter[trueIter] <- N  
    
  }
  
}

et1 <- Sys.time()
diffSecs <- as.numeric(difftime(et1,st1), units = "secs")

effectiveSize(mcmc(N_iter)) / diffSecs
qplot(1:niter, N_iter, geom = "line") + geom_hline(aes(yintercept = N_0))

qplot(N_iter)

# DIAGNOSTICS ----------

library(coda)

ess <- function(mcmc_output){
  
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i){
    mcmc(mcmc_output[i,])
  })
  mcmc_output_list_2 <- as.mcmc.list(mcmc_output_list)
  
  effectiveSize(mcmc_output_list_2)
  
}

qplot(1:niter, N_iter, geom = "line") + geom_hline(aes(yintercept = N_0))

str(mcmc(N_iter))

effectiveSize(mcmc(N_iter))

# SPEED --------

library(microbenchmark)

microbenchmark({
  s <- update_s_da(s, beta, z, CH, N, D, p0, sigma, Sigma_prop = diag(.1, nrow = 2),
                   Sigma_newpoint = diag(.1, nrow = 2), allTraps, R_traps, S_area, 
                   S_usage, trapsX, trapsY, a1, b1, a2, b2)
},{
  list_z <- update_z_cpp(z, s, N, D, phi, K, S_usage, trapsX, trapsY, p0, sigma)
  # list_z <- update_z(z, s, beta, N, phi)
  z <- list_z$z
  N <- list_z$N
  
  CH <- rbind(CH[1:D,],
              matrix(0, nrow = N - D, ncol = K))
})

####

pointsToPropose <- lambda_newPoints
x_new <- matrix(NA, pointsToPropose, 2)
for (l in 1:pointsToPropose) {
  x_new[l,] <- proposeNewPointNewMixture(probMixture[1], mixtureMeans, mixtureSd,
                                         a1, b1, a2, b2)
}
  

microbenchmark({
  # computePointsAdditionRatio(x_new, pointsToPropose, N,
  #                            s, beta, probMixture[1], mixtureMeans,
  #                            mixtureSd, normConst,
  #                            a1, b1, a2, b2, S_area,
  #                            D, K, p0, sigma,
  #                            S_usage, trapsX, trapsY)
  
  # log(probIndividualNotCaptured(s[100,], p0, sigma, D, N + l,
                                # K, S_usage, trapsX, trapsY))
  # logDensityBirth(x_new[1,], probMixture[1],
  #                 a1, b1, a2, b2,
  #                 mixtureMeans, mixtureSd, normConst)
},{
  exp(loglikelihood_xi_uncaptured(s[100,], K, S_usage, trapsX, trapsY, p0, sigma))  
})

logRatio_f = loglikelihood_addition(ptsToPropose, beta / S_area);

double logbirthDensity = 0;
for(int l = 0; l < pointsToPropose; l++){
  arma::vec x1_new = arma::conv_to<arma::vec>::from(x_new.row(l));
  logbirthDensity += logDensityBirth(x1_new, probMixture,
                                     a1, b1, a2, b2,
                                     mixtureMeans, mixtureSd, normConst);
}

double logDeathDensity = 0;
for(int l = 0; l < ptsToPropose; l++){
  logDeathDensity -= log(N + 1 + l);
}

double loglikelihood_xi_star = 0;
for(int l = 0; l < ptsToPropose; l++){
  arma::vec x_current = arma::conv_to<arma::vec>::from(x_new.row(l));
  loglikelihood_xi_star += log(probIndividualNotCaptured(x_current, p0, sigma, D, N + l,
                                                         K, S_usage, trapsX, trapsY));
}



