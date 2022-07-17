library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2)
setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
sourceCpp("rcpp/code.cpp")

# FUNCTIONS ---------------------------------------------------------------

# loglikelihood_p <- function(p0, alpha, CH1, CH2, s_1, s_2, trapsX, trapsY){
#   
#   likelihood <- 0
#   
#   for (i in 1:D1) {
#     for (j in 1:K) {
#       p <- p0 * exp(-alpha * ((s_1[i,1] - trapsX[j])^2 + (s_1[i,2] - trapsY[j])^2))
#       likelihood <- likelihood + dbinom(CH1[i,j], S, p, log = T)
#     }
#   }
#   
#   for (i in 1:D2) {
#     for (j in 1:K) {
#       p <- p0 * exp(-alpha * ((s_2[i,1] - trapsX[j])^2 + (s_2[i,2] - trapsY[j])^2))
#       likelihood <- likelihood + dbinom(CH2[i,j], S, p, log = T)
#     }
#   }
#   
#   likelihood
# }

update_p <- function(p0, alpha, s1, s2, CH1, CH2, trapsX, trapsY, 
                     S, D1, D2, N1, N2, alpha_0, sd_alpha,
                     sigma_p0_prop, sigma_alpha_prop){
  
  # update p0
  p0_star <- rnorm(1, p0, sd = sigma_p0_prop)
  if(p0_star > 0 & p0_star < 1){
    
    loglikelihood_p0_star <- loglikelihood_p_cpp(p0_star, alpha, CH1, CH2, s1, s2,
                                                 trapsX, trapsY, S, D1, D2, N1, N2)
    
    loglikelihood_p0 <- loglikelihood_p_cpp(p0, alpha, CH1, CH2, s1, s2,
                                            trapsX, trapsY, S, D1, D2, N1, N2)
    
    if(runif(1) < exp(loglikelihood_p0_star - loglikelihood_p0)){
      p0 <- p0_star
    }
    
  }
  
  # update alpha
  alpha_star <- rnorm(1, alpha, sd = sigma_alpha_prop)
  if(alpha_star > 0){
    
    loglikelihood_alpha_star <- loglikelihood_p_cpp(p0, alpha_star, CH1, CH2, s1, s2,
                                                    trapsX, trapsY, S, D1, D2, N1, N2)
    
    loglikelihood_alpha <- loglikelihood_p_cpp(p0, alpha, CH1, CH2, s1, s2,
                                               trapsX, trapsY, S, D1, D2, N1, N2)
    
    logprior_alpha_star <- dnorm(alpha_star, alpha_0, sd_alpha, log = T)
    
    logprior_alpha <- dnorm(alpha, alpha_0, sd_alpha, log = T)
    
    logposterior_alpha_star <- loglikelihood_alpha_star + logprior_alpha_star
    
    logposterior_alpha <- loglikelihood_alpha + logprior_alpha
    
    if(runif(1) < exp(loglikelihood_alpha_star - loglikelihood_alpha)){
      alpha <- alpha_star
    }
    
  }
  
  list("p0" = p0, "alpha" = alpha)
}

loglikelihood_xi_captured <- function(si, i, K, S, trapsX, trapsY, CH, p0, alpha){
  
  loglikelihood <- 0
  
  for (j in 1:K) {
    p <- p0 * exp(-alpha * ((si[1] - trapsX[j])^2 + (si[2] - trapsY[j])^2))
    loglikelihood <- loglikelihood + dbinom(CH[i,j], S, p, log = T)
  }  
  
  loglikelihood
  
}

loglikelihood_xi_uncaptured <- function(si, K, S, trapsX, trapsY, p0, alpha){
  
  loglikelihood <- 0
  
  for (j in 1:K) {
    p <- p0 * exp(-alpha * ((si[1] - trapsX[j])^2 + (si[2] - trapsY[j])^2))
    loglikelihood <- loglikelihood + dbinom(0, S, p, log = T)
  }  
  
  loglikelihood
  
}

update_s <- function(s1, s2, theta1, theta2, theta3, beta1, beta2, sigma_prop, N1, N2){
  
  log_hastings_ratio_den <- log_f_bivsoftcore_cpp(s1, s2, N1, N2,
                                                  beta1, beta2,
                                                  theta1, theta2, theta3)
  
  # move
  
  for (i in 1:N1) {
    
    old_xi <- s1[i,]
    xi <- mvrnorm(1, mu = old_xi, Sigma = diag(sigma_prop, nrow = 2))
    
    if(checkPointIsInRegion(xi)){
      s1[i,] <- xi
      
      log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_move(i - 1, 1, log_hastings_ratio_den,
                                                                 s1, s2,
                                                                 N1, N2,
                                                                 old_xi,
                                                                 theta1, theta2, theta3)
      
      if(i <= D1){
        loglikelihood_xi_star <- loglikelihood_xi_captured(xi, i, K, S, trapsX, trapsY, CH1, p0, alpha)
        loglikelihood_xi_old <- loglikelihood_xi_captured(old_xi, i, K, S, trapsX, trapsY, CH1, p0, alpha)
      } else {
        loglikelihood_xi_star <- loglikelihood_xi_uncaptured(xi, K, S, trapsX, trapsY, p0, alpha)
        loglikelihood_xi_old <- loglikelihood_xi_uncaptured(old_xi, K, S, trapsX, trapsY, p0, alpha)
      }
      
      logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
      logposterior_xi <- log_hastings_ratio_den + loglikelihood_xi_old
      
      if(runif(1) < exp(logposterior_xistar - logposterior_xi)){
        log_hastings_ratio_den <- log_hastings_ratio_num
      } else {
        s1[i,] <- old_xi
      }
    }
    
  }
  
  for (i in 1:N2) {
    
    old_xi <- s2[i,]
    xi <- mvrnorm(1, mu = old_xi, Sigma = diag(sigma_prop, nrow = 2))
    
    if(checkPointIsInRegion(xi)){
      s2[i,] <- xi
      
      log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_move(i - 1, 2, log_hastings_ratio_den,
                                                                 s1, s2,
                                                                 N1, N2,
                                                                 old_xi,
                                                                 theta1, theta2, theta3)
      
      if(i <= D2){
        loglikelihood_xi_star <- loglikelihood_xi_captured(xi, i, K, S, trapsX, trapsY, CH2, p0, alpha)
        loglikelihood_xi_old <- loglikelihood_xi_captured(old_xi, i, K, S, trapsX, trapsY, CH2, p0, alpha)  
      } else {
        loglikelihood_xi_star <- loglikelihood_xi_uncaptured(xi, K, S, trapsX, trapsY, p0, alpha)
        loglikelihood_xi_old <- loglikelihood_xi_uncaptured(old_xi, K, S, trapsX, trapsY, p0, alpha)
      }
      
      logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
      logposterior_xi <- log_hastings_ratio_den + loglikelihood_xi_old
      
      if(runif(1) < exp(logposterior_xistar - logposterior_xi)){
        log_hastings_ratio_den <- log_hastings_ratio_num
      } else {
        s2[i,] <- old_xi
      }
    }
    
  }
  
  # birth
  x1_new <- proposeNewPoint(sigma_prop)
  
  if(checkPointIsInRegion(x1_new)){
    
    log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_birth(x1_new, 1, log_hastings_ratio_den,
                                                                s1, s2,
                                                                N1, N2,
                                                                beta1, beta2,
                                                                theta1, theta2, theta3)
    
    loglikelihood_xi_star <- loglikelihood_xi_uncaptured(x1_new, K, S, trapsX, trapsY, p0, alpha)
    
    logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
    logposterior_xi <- log_hastings_ratio_den
    
    if(runif(1) < exp(logposterior_xistar - logposterior_xi)){
      log_hastings_ratio_den <- log_hastings_ratio_num
      N1 <- N1 + 1
      s1[N1,] <- x1_new
    } 
    
  }
  
  x2_new <- proposeNewPoint(sigma_prop)
  
  if(checkPointIsInRegion(x2_new)){
    
    log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_birth(x2_new, 2, log_hastings_ratio_den,
                                                                s1, s2,
                                                                N1, N2,
                                                                beta1, beta2,
                                                                theta1, theta2, theta3)
    
    loglikelihood_xi_star <- loglikelihood_xi_uncaptured(x2_new, K, S, trapsX, trapsY, p0, alpha)
    
    logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
    logposterior_xi <- log_hastings_ratio_den
    
    if(runif(1) < exp(logposterior_xistar - logposterior_xi)){
      log_hastings_ratio_den <- log_hastings_ratio_num
      N2 <- N2 + 1
      s2[N2,] <- x2_new
    } 
    
  }
  
  # death
  
  if(N1 > D1){
    itemToRemove <- sample((D1 + 1):N1, 1)
    xi <- s1[itemToRemove,]
    
    log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_death(itemToRemove - 1,
                                                                1, log_hastings_ratio_den,
                                                                s1, s2,
                                                                N1, N2,
                                                                beta1, beta2,
                                                                theta1, theta2, theta3)
    
    loglikelihood_xi_star <- - loglikelihood_xi_uncaptured(xi, K, S, trapsX, trapsY, p0, alpha)
    
    logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
    logposterior_xi <- log_hastings_ratio_den
    
    if(runif(1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
      s1[itemToRemove,] <- s1[N1,]
      log_hastings_ratio_den <- log_hastings_ratio_num
      N1 <- N1 - 1
    }
  }
  
  if(N2 > D2){
    itemToRemove <- sample((D2 + 1):N2, 1)
    xi <- s2[itemToRemove,]
    
    log_hastings_ratio_num = log_f_bivsoftcore_cpp_quick_death(itemToRemove - 1,
                                                               2, log_hastings_ratio_den,
                                                               s1, s2,
                                                               N1, N2,
                                                               beta1, beta2,
                                                               theta1, theta2, theta3)
    
    loglikelihood_xi_star <- - loglikelihood_xi_uncaptured(xi, K, S, trapsX, trapsY, p0, alpha)
    
    logposterior_xistar <- log_hastings_ratio_num + loglikelihood_xi_star
    logposterior_xi <- log_hastings_ratio_den
    
    if(runif(1) < exp(log_hastings_ratio_num - log_hastings_ratio_den)){
      s2[itemToRemove,] <- s2[N2,]
      log_hastings_ratio_den <- log_hastings_ratio_num
      N2 <- N2 - 1
    }
  }
  
  list("s1" = s1, "s2" = s2, 
       "N1" = N1, "N2" = N2)
}

log_prior <- function(theta1, theta2, theta3,
                      a_theta, b_theta){
  
  dgamma(theta1, a_theta, b_theta, log = T) + dgamma(theta2, a_theta, b_theta, log = T) + 
    dgamma(theta3, a_theta, b_theta, log = T)
  
}

hastings_ratio <- function(data1, data2, N1, N2, 
                           x_all, N_all,
                           beta1, beta2, 
                           theta1, theta2, theta3,
                           theta1_star, theta2_star, theta3_star,
                           a_theta, b_theta){
  
  ratio_constants <- importance_sampling_estimate(x_all, N_all,
                                                  beta1, beta2, 
                                                  theta1, theta2, theta3,
                                                  theta1_star, theta2_star, theta3_star)
  
  den <- log_prior(theta1, theta2, theta3,
                   a_theta, b_theta) + log_f_bivsoftcore_cpp(data1, data2, N1, N2, 
                                                             beta1, beta2, theta1, theta2, theta3)
  
  num <- log_prior(theta1_star, theta2_star, theta3_star,
                   a_theta, b_theta) + log_f_bivsoftcore_cpp(data1, data2, N1, N2, 
                                                             beta1, beta2, theta1_star, theta2_star, theta3_star)
  
  exp(num - den) * ratio_constants
}

importance_sampling_estimate <- function(x_all, N_all,
                                         beta1, beta2, 
                                         theta1, theta2, theta3,
                                         theta1_star, theta2_star, theta3_star){
  
  constant <- 0
  M <- dim(x_all)[1]
  
  for (m in 1:M) {
    
    N1 <- N_all[m,1]
    N2 <- N_all[m,2]
    x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
    x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
    
    num <- log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                 beta1, beta2, 
                                 theta1, theta2, theta3)
    
    den <- log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                 beta1, beta2, 
                                 theta1_star, theta2_star, theta3_star)
    
    constant <- constant + exp(num - den)
    
  }
  
  constant / M
}

# REAL DATA ---------------------------------------------------------------

usingTrueData <- T
load("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/Data/Shameer Data/data.rda")
S <- 30

# SIMULATED DATA ----------------------------------------------------------

usingTrueData <- F

theta1 <- .4
theta2 <- .4
theta3 <- .01
beta1 <- 100
beta2 <- 350

theta1_true <- theta1
theta2_true <- theta2
theta3_true <- theta3

p0 <- 1
alpha <- 1

p0_true <- p0
alpha_true <- alpha

list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta3, beta1, beta2, 2000, 
                                      Sigma_prop = diag(0.01, nrow = 2), 500, 100, 
                                      Sigma_newpoint = diag(1, nrow = 2))
s1 <- list_s1s2$data1
s2 <- list_s1s2$data2
N1 <- list_s1s2$N1
N2 <- list_s1s2$N2

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")

K <- nrow(traps)

# normalize traps locations
traps$X <- (traps$X - mean(traps$X)) / sd(traps$X)
traps$Y <- (traps$Y - mean(traps$Y)) / sd(traps$Y)

trapsX <- traps$X
trapsY <- traps$Y

CH_1 <- matrix(NA, nrow = N1, ncol = K)
for (i in 1:N1) {
  for (k in 1:K) {
    p <- p0 * exp(-alpha * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
    CH_1[i,k] <- rbinom(1, S, p)
  }
}
indexesCaughtIndividuals1 <- which(apply(CH_1,1,sum) != 0)
indexesUncaughtIndividuals1 <- which(apply(CH_1,1,sum) == 0)
D1 <- length(indexesCaughtIndividuals1)
CH_1 <- CH_1[indexesCaughtIndividuals1,]

CH_2 <- matrix(NA, nrow = N2, ncol = K)
for (i in 1:N2) {
  for (k in 1:K) {
    p <- p0 * exp(-alpha * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
    CH_2[i,k] <- rbinom(1, S, p)
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

# CLEAN DATA --------------------------------------------------------------

D1 <- nrow(CH_leopards)
D2 <- nrow(CH_tigers)

K <- ncol(CH_leopards) # number of traps

CH1 <- CH_leopards
CH2 <- CH_tigers

# normalize traps locations
traps$X <- (traps$X - mean(traps$X)) / sd(traps$X)
traps$Y <- (traps$Y - mean(traps$Y)) / sd(traps$Y)

trapsX <- traps$X
trapsY <- traps$Y

# PLOT THE DATA -----------------------------------------------------------

# leopards
ggplot(data = traps, aes(x = X, y = Y, color = ACTIVE_LEOPARD)) + geom_point(size = 2) + 
  theme(legend.title = element_blank()) + scale_color_discrete() + theme_bw() + theme(legend.title = element_blank())
# tigers
ggplot(data = traps, aes(x = X, y = Y)) + geom_point(aes(color = as.factor(ACTIVE_TIGER)), size = 2) + 
  theme(legend.title = element_blank()) + theme_bw()  + theme(legend.title = element_blank())

# PRIOR -------------------------------------------------------------------

Nmax <- 500

beta1 <- 200
beta2 <- 200

a_theta <- 1.5
b_theta <- 20

sigma_prop <- .01
Sigma_newpoint <- cov(cbind(trapsX, trapsY))

ggplot(data = NULL, aes(x = c(0,1))) + stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))

alpha_0 <- 1
sd_alpha <- .05

ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = alpha_0, sd = sd_alpha))

# MCMC --------------------------------------------------------------------

nburn <- 1000
niter <- 1000
nchain <- 1
nthin <- 2

# grid fornew points density
{
  gridLength_x <- seq(min(trapsX) - .05, 
                      max(trapsY) + .05, length.out = 100)
  gridLength_y <- seq(min(trapsY) - .05, 
                      max(trapsY) + .05, length.out = 100)
}

# variables for output
{
  
}

for(chain in 1:nchain) {
  
  # starting values
  {
    # capture probability parameter
    {
      if(usingTrueData){
        p0 <- .001
        alpha <- 1  
      } else {
        p0 <- p0_true
        alpha <- alpha_true
      }
      
      sigma_p0_prop <- .001
      sigma_alpha_prop <- .01
    }
    
    # population sizes
    {
      if(usingTrueData){
        N1 <- 2 * D1 - 10
        N2 <- 2 * D2 - 10  
      } else {
        N1 <- N1_0
        N2 <- N2_0
      }
    }
    
    # home centers
    {
      if(usingTrueData){
        
        s1 <- matrix(NA, nrow = Nmax, ncol = 2)
        for (i in 1:N1) {
          if(i <= D1){ # assign the home center to the average location of the traps where the indiv. was captured
            s1[i,] <- apply(traps[CH1[i,] == 1, 1:2], 2, mean)
            s1[i,] <- mvrnorm(1, mu = s1[i,], diag(.005, nrow = 2))
          } else { # otherwise a random point
            s1[i,] <- proposeNewPoint(Sigma_newpoint)
          }
        }
        s2 <- matrix(NA, nrow = Nmax, ncol = 2)
        for (i in 1:N2) {
          if(i <= D2){ # assign the home center to the average location of the traps where the indiv. was captured
            s2[i,] <- apply(traps[CH2[i,] == 1, 1:2], 2, mean)
            s2[i,] <- mvrnorm(1, mu = s2[i,], diag(.005, nrow = 2))
          } else { # otherwise a random point
            s2[i,] <- proposeNewPoint(Sigma_newpoint)
          }
        }
        
      } else {
        
        s1 <- matrix(NA, nrow = Nmax, ncol = 2)
        s2 <- matrix(NA, nrow = Nmax, ncol = 2)
        
        s1[1:N1,] <- s1_0
        s2[1:N2,] <- s2_0
        
      }
      
    }
    
    # point process parameters
    {
      if(usingTrueData){
        theta1 <- .2
        theta2 <- .2
        theta3 <- .2  
      } else {
        # theta1 <- theta1_true
        # theta2 <- theta2_true
        # theta3 <- theta3_true
        theta1 <- .2
        theta2 <- .2
        theta3 <- .2  
      }
      
    }
    
    # additional random variables
    {
      M <- 10
      
      x_all <- array(NA, dim = c(M, 2, Nmax, 2))
      N_all <- matrix(NA, nrow = M, ncol = 2)
      
      for (m in 1:M) {
        
        print(m)
        
        list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta3, beta1, beta2,
                                              niter = 500, 
                                              Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta1, 
                                              Sigma_newpoint = diag(1, nrow = 2))
        N1_sim <- list_sims$N1
        N2_sim <- list_sims$N2
        x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
        x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
        N_all[m,] <- c(N1_sim, N2_sim)
        
      }
      
    }
    
  }
  
  # output variables
  {
    alpha_iter <- rep(NA, niter)
    s_iter <- array(NA, dim = c(niter, 2, 3 * D1, 2))
    N_iter <- matrix(NA, nrow = niter, ncol = 2)
    p0_iter <- rep(NA, niter)
    params_iter <- matrix(NA, nrow = niter, ncol = 3)
    
    papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter <= nburn){
      print(iter) 
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print((iter - nburn)/nthin)  
    } 
    
    # UPDATE P ----------------------------------------------------------------
    
    list_p <- update_p(p0, alpha, s1, s2, CH1, CH2, trapsX, trapsY,
                       S, D1, D2, N1, N2, alpha_0, sd_alpha,
                       sigma_p0_prop, sigma_alpha_prop)
    p0 <- list_p$p0
    alpha <- list_p$alpha
    
    # print(p0)
    # print(alpha)
    
    # UPDATE HOME CENTERS -----------------------------------------------
    
    list_s <-  update_s_cpp(s1, s2, theta1, theta2, theta3, beta1, beta2, CH1, 
                            CH2, N1, N2, D1, D2, p0, alpha, S, trapsX, trapsY, 
                            Sigma_prop= diag(.01, nrow = 2), Sigma_newpoint)
    s1 <- list_s$s1
    s2 <- list_s$s2
    N1 <- list_s$N1
    N2 <- list_s$N2
    
    # UPDATE BETA, THETA ---------------------------------------------------------
    
    theta1_star <- rnorm(1, theta1, sd = .01)
    theta2_star <- rnorm(1, theta2, sd = .01)
    theta3_star <- rnorm(1, theta3, sd = .01)
    
    ratio <- hastings_ratio(s1, s2, N1, N2,
                            x_all, N_all,
                            beta1, beta2,
                            theta1, theta2, theta3,
                            theta1_star, theta2_star, theta3_star,
                            a_theta, b_theta)
    
    # print(ratio)
    
    if(theta1_star > 0 & theta2_star > 0 & theta3_star > 0){
      
      if(runif(1) < ratio){
        theta1 <- theta1_star
        theta2 <- theta2_star
        theta3 <- theta3_star
        
        x_all1 <- x_all[,1,,]
        x_all2 <- x_all[,2,,]
        
        list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 10,
                                      theta1, theta2, theta3, beta1, beta2, diag(0.01, nrow = 2), Sigma_newpoint)
        x_all[,1,,] <- list_xall$x_all1
        x_all[,2,,] <- list_xall$x_all2
        N_all <- list_xall$N_all
      }
      
    }
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn)/nthin
      s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
      s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
      p0_iter[trueIter] <- p0
      alpha_iter[trueIter] <- alpha
      N_iter[trueIter,] <- c(N1,N2)  
      params_iter[trueIter,] <- c(theta1, theta2, theta3)
      papangelou_density_iter[trueIter,,,] <- computeNewPointsDensity(gridLength_x, gridLength_y,
                                                                    s1, s2, N1, N2,
                                                                    theta1,  theta2,  theta3,
                                                                    beta1,  beta2)
    }
    
  }
  
  # Write results in MCMC output
  if(iter > nburn & (iter - nburn) %% nthin == 0) {
    
  }
  
}

# DIAGNOSTICS -------------------------------------------------------------

qplot(1:niter, params_iter[,1])
qplot(1:niter, params_iter[,2])
qplot(1:niter, params_iter[,3])

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

qplot(1:niter, p0_iter)

qplot(1:niter, N_iter[,1])
qplot(1:niter, N_iter[,2])

qplot(1:niter, s_iter[,1,1,1])
qplot(1:niter, s_iter[,1,21,2])
qplot(1:niter, s_iter[,2,11,2])

# PLOTS -------------------------------------------------------------------

# plots for gamma
{
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  p1 <- ggplot(data = NULL, aes(x = params_iter[,1], y = ..density..)) + 
    geom_histogram(binwidth = .015, fill = "cornsilk", color = "black") +
    geom_vline(aes(xintercept = mean(params_iter[,1])), color = "red", size = 1.5) +
    theme_bw() + scale_x_continuous(name = expression(gamma[1]), limits = c(0,.7)) +
    scale_y_continuous(name = "")
  
  p2 <- ggplot(data = NULL, aes(x = params_iter[,2], y = ..density..)) + 
    geom_histogram(binwidth = .015, fill = "cornsilk", color = "black") +
    geom_vline(aes(xintercept = mean(params_iter[,2])), color = "red", size = 1.5) +
    theme_bw() + scale_x_continuous(name = expression(gamma[2]), limits = c(0,.7)) +
    scale_y_continuous(name = "")
  
  p3 <- ggplot(data = NULL, aes(x = params_iter[,3], y = ..density..)) + 
    geom_histogram(binwidth = .015, fill = "cornsilk", color = "black") +
    geom_vline(aes(xintercept = mean(params_iter[,3])), color = "red", size = 1.5) +
    theme_bw() + scale_x_continuous(name = expression(gamma[3]), limits = c(0,.7)) +
    scale_y_continuous(name = "")
  
  multiplot(p1, p2, p3, cols = 1)
  
  p1 <- qplot(1:niter, params_iter[,1])
  p2 <- qplot(1:niter, params_iter[,2])
  p3 <- qplot(1:niter, params_iter[,3])
  
  multiplot(p1, p2, p3, p4, cols = 1)
  
  p4 <- qplot(1:niter, p0_iter)
}

# plots for N
{
  p1 <- ggplot(data = NULL, aes(x = N_iter[,1], y = ..density..)) + 
    geom_histogram(binwidth = 1, fill = "cornsilk", color = "black") +
    geom_vline(aes(xintercept = mean(N_iter[,1])), color = "red", size = 1.5) +
    theme_bw() + scale_x_continuous(name = expression(N[1]), limits = c(20,120)) +
    scale_y_continuous(name = "")
  
  p2 <- ggplot(data = NULL, aes(x = N_iter[,2], y = ..density..)) + 
    geom_histogram(binwidth = 1, fill = "cornsilk", color = "black") +
    geom_vline(aes(xintercept = mean(N_iter[,2])), color = "red", size = 1.5) +
    theme_bw() + scale_x_continuous(name = expression(N[2]), limits = c(20,120)) +
    scale_y_continuous(name = "")
  
  multiplot(p1, p2)
  
  p1 <- qplot(1:niter, N_iter[,1])
  p2 <- qplot(1:niter, N_iter[,2])
  
  multiplot(p1, p2, p4)
}


# CREATE BOUNDARIES -------------------------------------------------------

pointsBoundaries <- matrix(NA, nrow = 10, ncol = 2)
pointsBoundaries[1,] <- c(-2.6, -0.3)
pointsBoundaries[2,] <- c(-.39, 2.05)
pointsBoundaries[3,] <- c(.85, -2.65)
pointsBoundaries[4,] <- c(1.9, .45)
pointsBoundaries[5,] <- c(1.5, 1.9)
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

# DENSITY PLOT ------------------------------------------------------------

N1 <- N_iter[niter,1]
N2 <- N_iter[niter,2]
s1 <- s_iter[niter,1,1:N1,]
s2 <- s_iter[niter,2,1:N2,]

ggplot() + geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black") 

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
  scale_color_continuous(low = "white", high = "blue") +
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
  scale_color_continuous(low = "white", high = "blue") +
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

ggplot() +
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2], color = "blue"),
             size = 50 * binwidth) +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2], color = "red"),
             size = 50 * binwidth) +
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "black") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "black")  


