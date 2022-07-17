library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2); library(reshape2)
library(SpatialFunctions)
library(foreach); library(doParallel)

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)

setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
# setwd("C:/Users/ad603/Dropbox/R Folder/PhD/Spatial")
# sourceCpp("rcpp/code_basic.cpp")

# FUNCTIONS ---------------------------------------------------------------

log_prior <- function(theta1, theta2, theta12, 
                      a_theta, b_theta){

  dgamma(theta1, a_theta, b_theta, log = T) + 
    dgamma(theta2, a_theta, b_theta, log = T) +
    dgamma(theta12, a_theta, b_theta, log = T) 

}

hastings_ratio <- function(data1, data2, 
                           N1, N2, 
                           x_all, N_all,
                           beta1, beta2,
                           beta1_star, beta2_star, 
                           theta1, theta2, theta12,
                           theta1_star, theta2_star, theta12_star,
                           a_theta, b_theta){
  
  ratio_constants <- importance_sampling_estimate_foreach(x_all, N_all,
                                                  beta1, beta2,
                                                  beta1_star, beta2_star,
                                                  theta1, theta2, theta12,
                                                  theta1_star, theta2_star, theta12_star)
  
  den <- log_prior(theta1, theta2, theta12,
                   a_theta, b_theta) + 
    log_f_bivsoftcore_cpp(data1, data2, N1, N2,
                           beta1, beta2, 
                           theta1, theta2, theta12)
  
  num <- log_prior(theta1_star, theta2_star, theta12_star,
                   a_theta, b_theta) + 
    log_f_bivsoftcore_cpp(data1, data2, N1, N2, 
                           beta1_star, beta2_star, 
                           theta1_star, theta2_star, theta12_star)
  
  exp(num - den) * ratio_constants
}

importance_sampling_estimate_foreach <- function(x_all, N_all,
                                         beta1, beta2,
                                         beta1_star, beta2_star,
                                         theta1, theta2, theta12,
                                         theta1_star, theta2_star, theta12_star){
  
  r <- foreach(m = 1:M, .combine=c, 
               .packages = "SpatialFunctions") %dopar% {
                 
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
                 x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
                 
                 num <- log_f_bivsoftcore_cpp(x1, x2, 
                                              N1, N2,  
                                              beta1, beta2, 
                                              theta1, theta2, theta12) 
                 
                 den <- log_f_bivsoftcore_cpp(x1, x2,
                                              N1, N2,  
                                              beta1_star, beta2_star, 
                                              theta1_star, theta2_star, theta12_star)
                 
                 exp(num - den)         
               }
  
  mean(r)
}

update_x_all_foreach <- function(x_all, N_all, niter,
                         theta1, theta2, theta12,
                         beta1, beta2, Sigma_prop, a, b){
  
  M <- dim(x_all)[1]
  
  r <- foreach(m = 1:M, .combine='c', .multicombine=F,
               .packages = "SpatialFunctions") %dopar% {
                 
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 x1 <- x_all[m,1,,]
                 x2 <- x_all[m,2,,]
                 
                 list_sims <- simulate_bivsoftcore_from_startingpoint(x1, x2, N1, N2, 
                                                                      theta1, theta2, theta12,
                                                                      beta1, beta2, 
                                                                      niter, Sigma_prop, 
                                                                      a, b)
                 
                 list(list("N1" = list_sims$N1,
                           "N2" = list_sims$N2,
                           "data1" = list_sims$data1,
                           "data2" = list_sims$data2))
                 
               }
  
  for (m in 1:M) {
    
    list_r <- r[[m]]
    
    N_all[m,1] <- list_r$N1
    N_all[m,2] <- list_r$N2
    x_all[m,1,,] <- list_r$data1
    x_all[m,2,,] <- list_r$data2
    
  }
  
  list("x_all" = x_all, "N_all" = N_all)
}

# SIMULATING DATA ---------------------------------------------------------

usingTrueData <- F

theta1 <- .001
theta2 <- .001
theta12 <- .002
beta1 <- 200
beta2 <- 200

beta1_true <- beta1
beta2_true <- beta2
theta1_true <- theta1
theta2_true <- theta2
theta12_true <- theta12

a <- 0
b <- 10

list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, 
                                      beta1, beta2, 5000, 
                                      Sigma_prop = diag(0.01, nrow = 2), 500, 100, 
                                      a, b)
s1 <- list_s1s2$data1
s2 <- list_s1s2$data2
N1 <- list_s1s2$N1
N2 <- list_s1s2$N2

c(N1,N2)

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black") 

# PRIOR -------------------------------------------------------------------

Nmax <- 700

beta1 <- beta1_true
beta2 <- beta2_true

a_theta <- 0.001
b_theta <- 0.001

ggplot(data = NULL, aes(x = c(0,1))) + stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))

sigma_prop <- .05

beta_proposal <- .2

epsilon_theta <- .001^2
epsilon_beta <- 1^2

# MCMC --------------------------------------------------------------------

nburn <- 5000
niter <- 5000
nchain <- 1
nthin <- 2

# variables for output
{
  
}

for(chain in 1:nchain) {
  
  # starting values
  {
    
    # point process parameters
    {
      beta1 <- beta1_true
      beta2 <- beta2_true
      
      theta1 <- theta1_true
      theta2 <- theta2_true
      theta12 <- theta12_true
    }
    
    # additional random variables
    {
      M <- 50
      
      x_all <- array(NA, dim = c(M, 2, Nmax, 2))
      N_all <- matrix(NA, nrow = M, ncol = 2)
      
      for (m in 1:M) {
        
        print(m)
        
        list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, 
                                               beta1, beta2, 
                                              niter = 10000, 
                                              Sigma_prop = diag(.05, nrow = 2), Nmax = Nmax, 
                                              lambda = beta1, a, b)
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
    params_iter <- matrix(NA, nrow = niter, ncol = 5)
    acceptances_iter <- matrix(0, nrow = nburn, ncol = 5)
    
    theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter < nburn){
      print(iter)
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print((iter - nburn)/nthin) 
      print(paste0("Theta 1 = ",theta1," / Theta 2 = ",theta2," / Theta 3 = ",theta3))
    }
    
    # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
    
    if(iter > 200){
      Sigma_theta <- cov(theta_all[1:(iter - 1),c(1,3,5)])
      Sigma_proposal <- (1 - beta_proposal) * (2.38)^2 * Sigma_theta / 3 +
        beta_proposal * diag(c(epsilon_beta, epsilon_theta, epsilon_theta), 
                             nrow = 3) / 3
    } else {
      Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta, epsilon_theta), 
                             nrow = 3) / 3 
    }
    
    proposed_theta <- mvrnorm(1, c(beta1, theta1, theta12),
                              Sigma = Sigma_proposal)
    beta1_star <- proposed_theta[1]
    theta1_star <- proposed_theta[2]
    theta12_star <- proposed_theta[3]
    beta1_star <- beta1
    
    ratio <- hastings_ratio(s1, s2, N1, N2,
                            x_all, N_all,
                            beta1, beta2,
                            beta1_star, beta2,
                            theta1, theta2, theta12,
                            theta1_star, theta2, theta12_star,
                            a_theta, b_theta)
    
    if(beta1_star > 0 & theta1_star > 0 & theta12_star > 0){
      
      if(runif(1) < ratio){
        beta1 <- beta1_star
        theta1 <- theta1_star
        theta12 <- theta12_star
        
        list_xall <- update_x_all_foreach(x_all, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2), 
                                          a, b)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
        
      }
      
    }
    
    # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
    
    if(iter > 200){
      Sigma_theta <- cov(theta_all[1:(iter - 1),c(2,4,5)])
      Sigma_proposal <- (1 - beta_proposal) * (2.38)^2 * Sigma_theta / 3 +
        beta_proposal * diag(c(epsilon_beta, epsilon_theta, epsilon_theta), 
                             nrow = 3) / 3
    } else {
      Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta, epsilon_theta), 
                             nrow = 3) / 3 
    }
    
    proposed_theta <- mvrnorm(1, c(beta2, theta2, theta12),
                              Sigma = Sigma_proposal)
    beta2_star <- proposed_theta[1]
    theta2_star <- proposed_theta[2]
    theta12_star <- proposed_theta[3]
    beta2_star <- beta2
    
    ratio <- hastings_ratio(s1, s2, N1, N2,
                            x_all, N_all,
                            beta1, beta2,
                            beta1, beta2_star,
                            theta1, theta2, theta12,
                            theta1, theta2_star, theta12_star,
                            a_theta, b_theta)
    
    if(beta2_star > 0 & theta2_star > 0 & theta12_star > 0){
      
      if(runif(1) < ratio){
        beta2 <- beta2_star
        theta2 <- theta2_star
        theta12 <- theta12_star
        
        list_xall <- update_x_all_foreach(x_all, N_all, 100,
                                          theta1, theta2, theta12,
                                          beta1, beta2, diag(0.05, nrow = 2),
                                          a, b)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
        
        # acceptances_theta[iter,1] <- 1
      }
      
    }
    
    # WRITE RESULTS ---------------
    
    theta_all[iter,] <- c(beta1, beta2, theta1, theta2, theta12)
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn)/nthin
      params_iter[trueIter,] <- c(theta1, theta2, theta12, beta1, beta2)
    }
    
  }
  
  # Write results in MCMC output
  if(iter > nburn & (iter - nburn) %% nthin == 0) {
    
    
    
  }
  
}

# DIAGNOSTICS -------------------------------------------------------------

trueparams <- c(theta1_true,
                theta2_true,
                theta3_true,
                theta12_true,
                theta13_true,
                theta23_true)

qplot(1:niter, params_iter[,1])
qplot(1:niter, params_iter[,2])
qplot(1:niter, params_iter[,3])
qplot(1:niter, params_iter[,4])
qplot(1:niter, params_iter[,5])
qplot(1:niter, params_iter[,6])

theta_pci <- t(apply(params_iter, 2, function(x){
  quantile(x, probs = c(0.025,.5,0.975))
}))[1:6,]

# rownames(theta_pci) <- 
colnames(theta_pci) <- c("Min","Mean","Max")#c("theta1","theta2","theta3","theta12","theta13","theta23")

theta_pci_long <- melt(theta_pci)

ggplot() + geom_errorbar(data = as.data.frame(theta_pci), aes(x = 1:6, ymin = Min, ymax = Max)) +
  geom_point(data = NULL,aes(x = 1:6, y = trueparams), size = 2) + xlab("Parameters") + ylab("") + theme_bw()
