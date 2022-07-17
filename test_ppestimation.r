library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2); library(reshape2)
library(foreach); library(doParallel)
# library(SpatialFunctionsCR); 

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
# stopCluster(cl)

setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
sourceCpp("CR Model/code.cpp")
source('CR Model/functions.r', echo=TRUE)


# SIMULATION --------------------------------------------------------------

usingTrueData <- T

# starting params
{
  theta1 <- exp(-100)
  theta2 <- exp(-100)
  theta12 <- exp(-100)
  beta1 <- 200
  beta2 <- 300
  
  beta1_true <- beta1
  beta2_true <- beta2
  theta1_true <- theta1
  theta2_true <- theta2
  theta12_true <- theta12
  
}

# data simulation
{
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 1000, beta1, 
                                        Sigma_newpoint = diag(1, nrow = 2),
                                        a1, b1, a2, b2)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  
  s1 <- s1[1:N1,]
  s2 <- s2[1:N2,]
  
  print(N1)
  print(N2)
  ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
    geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
  
}

# PRIOR -------------------------------------------------------------------

Nmax <- 1000

# interaction parameter
a_theta <- .01
b_theta <- .01
ggplot(data = NULL, aes(x = c(0,.3))) + stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))

beta_proposal <- .1

# epsilon_theta1 <- .00005^2
# epsilon_theta12 <- .00005^2
epsilon_logtheta1 <- .00005^2
epsilon_logtheta12 <- .00005^2
epsilon_beta <- 2^2

# intensity parameter
meanBeta <- 300
varBeta <- 20000
b_beta <- meanBeta / varBeta
a_beta <- meanBeta * b_beta

ggplot(data = NULL, aes(x = c(0,300))) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta))

# proposal value
sigma_prop <- .01
Sigma_newpoint <- cov(cbind(trapsX, trapsY))

# capture probability params

{
  a_p <- 1
  b_p <- 1
  
  sigma_0_1 <- 1.5 / trapsX_sd
  sd_sigma_1 <- 2
  
  sigma_0_2 <- 2.5 / trapsX_sd
  sd_sigma_2 <- 2
  
  ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1))
  ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_2, sd = sd_sigma_2))
  
  sigma_p0_prop <- .0001
  sigma_sigma_prop <- .005
}

# MCMC --------------------------------------------------------------------

nburn <- 3000
niter <- 3000
nchain <- 1
nthin <- 1
# nburn <- 7500
# niter <- 7500
# nchain <- 1
# nthin <- 3

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
    
    # point process parameters
    {
      if(usingTrueData){
        theta1 <- .002
        theta2 <- .002
        theta12 <- .002  
        beta1 <- 200
        beta2 <- 300
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
      
      params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
      
      # M <- 50
      # 
      # x_all <- array(NA, dim = c(M, 2, Nmax, 2))
      # N_all <- matrix(NA, nrow = M, ncol = 2)
      # 
      # for (m in 1:M) {
      #   
      #   # print(m)
      #   
      #   list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
      #                                         niter = 2000, 
      #                                         Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta1, 
      #                                         Sigma_newpoint = diag(1, nrow = 2),
      #                                         a1, b1, a2, b2)
      #   N1_sim <- list_sims$N1
      #   N2_sim <- list_sims$N2
      #   x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
      #   x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
      #   N_all[m,] <- c(N1_sim, N2_sim)
      #   
      # }
      
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
    
    acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
    theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter <= nburn){
      print(paste0("Burn In Iteration = ",iter)) 
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print(paste0("Iteration = ",(iter - nburn)/nthin))
    } 
    
    # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
    
    if(iter > 400){
      Sigma_n <- cov(params_values[1:(iter-1),c(1,3,5)])
      Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
        beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
    } else {
      Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta1, epsilon_logtheta12), nrow = 3) / 3
    }

    proposed_theta <- mvrnorm(1, c(beta1, log(theta1), log(theta12)),
                              Sigma = Sigma_proposal)
    beta1_star <- proposed_theta[1]
    theta1_star <- exp(proposed_theta[2])
    theta12_star <- exp(proposed_theta[3])
    # Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta1, epsilon_theta12),
    #                        nrow = 3) / 3
    #
    # proposed_theta <- mvrnorm(1, c(beta1, theta1, theta12),
    #                           Sigma = Sigma_proposal)
    # beta1_star <- proposed_theta[1]
    # theta1_star <- proposed_theta[2]
    # theta12_star <- proposed_theta[3]

    ratio <- hastings_ratio(s1, s2, N1, N2,
                            x_all, N_all,
                            beta1, beta2,
                            beta1_star, beta2,
                            theta1, theta2, theta12,
                            theta1_star, theta2, theta12_star,
                            a_theta, b_theta,
                            a_beta, b_beta)
    # print(ratio)
    if(beta1_star > 0 & theta1_star > 0.000 & theta12_star > 0.0000){

      if(runif(1) < ratio){
        beta1 <- beta1_star
        theta1 <- theta1_star
        theta12 <- theta12_star

        list_xall <- update_x_all_foreach(x_all, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2),
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          a1, b1, a2, b2)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all

      }

    }

    params_values[iter,c(1,3,5)] <- c(beta1, log(theta1), log(theta12))
    
    # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
    
    if(iter > 400){
      Sigma_n <- cov(params_values[1:(iter-1),c(2,4,5)])
      Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
        beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
    } else {
      Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta1, epsilon_logtheta12),
                             nrow = 3) / 3
    }

    proposed_theta <- mvrnorm(1, c(beta2, log(theta2), log(theta12)),
                              Sigma = Sigma_proposal)
    beta2_star <- proposed_theta[1]
    theta2_star <- exp(proposed_theta[2])
    theta12_star <- exp(proposed_theta[3])
    # Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta1, epsilon_theta12),
    #                        nrow = 3) / 3
    #
    # proposed_theta <- mvrnorm(1, c(beta2, theta2, theta12),
    #                           Sigma = Sigma_proposal)
    # beta2_star <- proposed_theta[1]
    # theta2_star <- proposed_theta[2]
    # theta12_star <- proposed_theta[3]

    ratio <- hastings_ratio(s1, s2, N1, N2,
                            x_all, N_all,
                            beta1, beta2,
                            beta1, beta2_star,
                            theta1, theta2, theta12,
                            theta1, theta2_star, theta12_star,
                            a_theta, b_theta,
                            a_beta, b_beta)

    if(beta2_star > 0 & theta2_star > 0 & theta12_star > 0){

      if(runif(1) < ratio){
        beta2 <- beta2_star
        theta2 <- theta2_star
        theta12 <- theta12_star

        list_xall <- update_x_all_foreach(x_all, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2),
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          a1, b1, a2, b2)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all

        # acceptances_theta[iter,1] <- 1
      }

    }

    params_values[iter,c(2,4,5)] <- c(beta2, log(theta2), log(theta12))
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      
      trueIter <- (iter - nburn)/nthin
      
      params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
    }
    
  }
  
  # Write results in MCMC output
  if(iter > nburn & (iter - nburn) %% nthin == 0) {
    
  }
  
}

