library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2); library(reshape2)
library(SpatialFunctions)
setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
# setwd("C:/Users/ad603/Dropbox/R Folder/PhD/Spatial")
# sourceCpp("rcpp/code_basic.cpp")

# FUNCTIONS ---------------------------------------------------------------

log_prior <- function(theta1, theta2, theta3,
                      theta12, theta13, theta23,
                      a_theta, b_theta){

  dgamma(theta1, a_theta, b_theta, log = T) + 
    dgamma(theta2, a_theta, b_theta, log = T) +
    dgamma(theta3, a_theta, b_theta, log = T) + 
    dgamma(theta12, a_theta, b_theta, log = T) + 
    dgamma(theta13, a_theta, b_theta, log = T) +
    dgamma(theta23, a_theta, b_theta, log = T)

}

hastings_ratio <- function(data1, data2, data3, 
                           N1, N2, N3, 
                           x_all, N_all,
                           beta1, beta2, beta3, 
                           beta1_star, beta2_star, beta3_star, 
                           theta1, theta2, theta3,
                           theta1_star, theta2_star, theta3_star,
                           theta12, theta13, theta23,
                           theta12_star, theta13_star, theta23_star,
                           a_theta, b_theta){
  
  # ratio_constants <- importance_sampling_estimate(x_all, N_all,
  #                                                 beta1, beta2, beta3,
  #                                                 theta1, theta2, theta3,
  #                                                 theta1_star, theta2_star, theta3_star,
  #                                                 theta12, theta13, theta23,
  #                                                 theta12_star, theta13_star, theta23_star)
  
  ratio_constants <- importance_sampling_estimate_foreach(x_all, N_all,
                                                  beta1, beta2, beta3,
                                                  theta1, theta2, theta3,
                                                  theta1_star, theta2_star, theta3_star,
                                                  theta12, theta13, theta23,
                                                  theta12_star, theta13_star, theta23_star)
  
  den <- log_prior(theta1, theta2, theta3,
                   theta12, theta13, theta23,
                   a_theta, b_theta) + 
    log_f_trivsoftcore_cpp(data1, data2, data3, N1, N2, N3, 
                           beta1, beta2, beta3, 
                           theta1, theta2, theta3, 
                           theta12, theta13, theta23)
  
  num <- log_prior(theta1_star, theta2_star, theta3_star,
                   theta12_star, theta13_star, theta23_star,
                   a_theta, b_theta) + 
    log_f_trivsoftcore_cpp(data1, data2, data3, N1, N2, N3, 
                           beta1, beta2, beta3, 
                           theta1_star, theta2_star, theta3_star, 
                           theta12_star, theta13_star, theta23_star)
  
  exp(num - den) * ratio_constants
}

importance_sampling_estimate <- function(x_all, N_all,
                                         beta1, beta2, beta3,
                                         theta1, theta2, theta3,
                                         theta1_star, theta2_star, theta3_star,
                                         theta12, theta13, theta23,
                                         theta12_star, theta13_star, theta23_star){
  
  constant <- 0
  M <- dim(x_all)[1]
  
  for (m in 1:M) {
    
    N1 <- N_all[m,1]
    N2 <- N_all[m,2]
    N3 <- N_all[m,3]
    x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
    x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
    x3 <- matrix(x_all[m,3,seq_len(N3),], nrow = N3)
    
    num <- log_f_trivsoftcore_cpp(x1, x2, x3, 
                                  N1, N2, N3, 
                                  beta1, beta2, beta3, 
                                  theta1, theta2, theta3, 
                                  theta12, theta13, theta23) 
    
    den <- log_f_trivsoftcore_cpp(x1, x2, x3, 
                                  N1, N2, N3, 
                                  beta1, beta2, beta3, 
                                  theta1_star, theta2_star, theta3_star, 
                                  theta12_star, theta13_star, theta23_star)
    
    constant <- constant + exp(num - den)
    
  }
  
  constant / M
}

importance_sampling_estimate_foreach <- function(x_all, N_all,
                                         beta1, beta2, beta3,
                                         theta1, theta2, theta3,
                                         theta1_star, theta2_star, theta3_star,
                                         theta12, theta13, theta23,
                                         theta12_star, theta13_star, theta23_star){
  
  r <- foreach(m = 1:M, .combine=c, 
               .packages = "SpatialFunctions") %dopar% {
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 N3 <- N_all[m,3]
                 x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
                 x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
                 x3 <- matrix(x_all[m,3,seq_len(N3),], nrow = N3)
                 
                 num <- log_f_trivsoftcore_cpp(x1, x2, x3, 
                                               N1, N2, N3, 
                                               beta1, beta2, beta3, 
                                               theta1, theta2, theta3, 
                                               theta12, theta13, theta23) 
                 
                 den <- log_f_trivsoftcore_cpp(x1, x2, x3,
                                               N1, N2, N3,
                                               beta1, beta2, beta3,
                                               theta1_star, theta2_star, theta3_star,
                                               theta12_star, theta13_star, theta23_star)
                 
                 exp(num - den)         
               }
  
  mean(r)
}

update_x_all <- function(x_all1, x_all2,  x_all3, N_all, niter,
                         theta1, theta2, theta3, theta12, theta13, theta23,
                         beta1, beta2, beta3, Sigma_prop, a, b){
  
  M <- dim(x_all)[1]
  
  for (m in 1:M) {
    
    N1 <- N_all[m,1]
    N2 <- N_all[m,2]
    N3 <- N_all[m,3]
    x1 <- x_all[m,1,,]
    x2 <- x_all[m,2,,]
    x3 <- x_all[m,3,,]
    
    list_sims <- simulate_trivsoftcore_from_startingpoint(x1, x2, x3, N1, N2, N3, 
                                                         theta1, theta2, theta3,
                                                         theta12, theta13, theta23,
                                                         beta1, beta2, beta3,
                                                         niter, Sigma_prop, 
                                                         a, b)
    
    
    N1 <- list_sims$N1
    N2 <- list_sims$N2
    N3 <- list_sims$N3
    x_all[m,1,,] <- list_sims$data1
    x_all[m,2,,] <- list_sims$data2
    x_all[m,3,,] <- list_sims$data3
    
  }
  
  
  list("x_all" = x_all, "N_all" = N_all)
}

update_x_all_foreach <- function(x_all, N_all, niter,
                         theta1, theta2, theta3, theta12, theta13, theta23,
                         beta1, beta2, beta3, Sigma_prop, a, b){
  
  M <- dim(x_all)[1]
  
  r <- foreach(m = 1:M, .combine='c', .multicombine=F,
               .packages = "SpatialFunctions") %dopar% {
                 
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 N3 <- N_all[m,3]
                 x1 <- x_all[m,1,,]
                 x2 <- x_all[m,2,,]
                 x3 <- x_all[m,3,,]
                 
                 list_sims <- simulate_trivsoftcore_from_startingpoint(x1, x2, x3, N1, N2, N3, 
                                                                       theta1, theta2, theta3,
                                                                       theta12, theta13, theta23,
                                                                       beta1, beta2, beta3,
                                                                       niter, Sigma_prop, 
                                                                       a, b)
                 
                 list(list("N1" = list_sims$N1,
                      "N2" = list_sims$N2,
                      "N3" = list_sims$N3,
                      "data1" = list_sims$data1,
                      "data2" = list_sims$data2,
                      "data3" = list_sims$data3))
                 
               }
  
  for (m in 1:M) {
    
    list_r <- r[[m]]
    
    N_all[m,1] <- list_r$N1
    N_all[m,2] <- list_r$N2
    N_all[m,3] <- list_r$N3
    x_all[m,1,,] <- list_r$data1
    x_all[m,2,,] <- list_r$data2
    x_all[m,3,,] <- list_r$data3
    
  }
  
  list("x_all" = x_all, "N_all" = N_all)
}

# SIMULATING DATA ---------------------------------------------------------

usingTrueData <- F

theta1 <- .01
theta2 <- .1
theta3 <- .2
theta12 <- .1
theta23 <- .01
theta13 <- .2
beta1 <- 25
beta2 <- 25
beta3 <- 25

beta1_true <- beta1
beta2_true <- beta2
beta3_true <- beta3
theta1_true <- theta1
theta2_true <- theta2
theta3_true <- theta3
theta12_true <- theta12
theta23_true <- theta23
theta13_true <- theta13

a <- 0
b <- 10

list_s1s2s3 <- simulate_trivsoftcore_cpp(theta1, theta2, theta3, theta12, theta13, theta23, 
                                         beta1, beta2, beta3, 2000, 
                                      Sigma_prop = diag(0.01, nrow = 2), 500, 100, 
                                      a, b)
s1 <- list_s1s2s3$data1
s2 <- list_s1s2s3$data2
s3 <- list_s1s2s3$data3
N1 <- list_s1s2s3$N1
N2 <- list_s1s2s3$N2
N3 <- list_s1s2s3$N3

c(N1,N2,N3)
sum(N1, N2, N3)

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]
s3 <- s3[1:N3,]

ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black") +
  geom_point(data = NULL, aes(x = s3[,1], y = s3[,2]), color = "blue")

# PRIOR -------------------------------------------------------------------

Nmax <- 700

beta1 <- beta1_true
beta2 <- beta2_true
beta3 <- beta3_true

a_theta <- 0.001
b_theta <- 0.001

ggplot(data = NULL, aes(x = c(0,1))) + stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))

sigma_prop <- .05

# MCMC --------------------------------------------------------------------

nburn <- 5000
niter <- 10
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
      beta3 <- beta3_true
      
      theta1 <- theta1_true
      theta2 <- theta2_true
      theta3 <- theta3_true
      theta13 <- theta13_true
      theta12 <- theta12_true
      theta23 <- theta23_true
    }
    
    # additional random variables
    {
      M <- 100
      
      x_all <- array(NA, dim = c(M, 3, Nmax, 2))
      N_all <- matrix(NA, nrow = M, ncol = 3)
      
      for (m in 1:M) {
        
        print(m)
        
        list_sims <- simulate_trivsoftcore_cpp(theta1, theta2, theta3, 
                                               theta12, theta13, theta23, 
                                               beta1, beta2, beta3,
                                              niter = 1000, 
                                              Sigma_prop = diag(.05, nrow = 2), Nmax = Nmax, 
                                              lambda = beta1, a, b)
        N1_sim <- list_sims$N1
        N2_sim <- list_sims$N2
        N3_sim <- list_sims$N3
        x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
        x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
        x_all[m,3,1:N3_sim,] <- list_sims$data3[1:N3_sim,]
        N_all[m,] <- c(N1_sim, N2_sim, N3_sim)
        
      }
      
    }
    
  }
  
  # output variables
  {
    params_iter <- matrix(NA, nrow = niter, ncol = 9)
    acceptances_iter <- matrix(0, nrow = nburn, ncol = 9)
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter < nburn){
      print(iter)
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print((iter - nburn)/nthin) 
      print(paste0("Theta 1 = ",theta1," / Theta 2 = ",theta2," / Theta 3 = ",theta3))
    }
    
    # UPDATE BETA ---------------------------------------------------------
    
    beta1_star <- rnorm(1, beta1, sd = .01)
    beta2_star <- rnorm(1, beta2, sd = .01)
    beta3_star <- rnorm(1, beta3, sd = .01)

    ratio <- hastings_ratio(s1, s2, s3, 
                            N1, N2, N3, 
                            x_all, N_all,
                            beta1, beta2, beta3, 
                            beta1_star, beta2_star, beta3_star, 
                            theta1, theta2, theta3,
                            theta1, theta2, theta3,
                            theta12, theta13, theta23,
                            theta12, theta13, theta23,
                            a_theta, b_theta)

    if(beta1_star > 0 & beta2_star > 0 & beta3_star > 0){

      if(runif(1) < ratio){
        beta1 <- beta1_star
        beta2 <- beta2_star
        beta3 <- beta3_star

        # x_all1 <- x_all[,1,,]
        # x_all2 <- x_all[,2,,]
        # x_all3 <- x_all[,3,,]
        # 
        # list_xall <- update_x_all_tri_cpp(x_all1, x_all2,  x_all3, N_all, 10,
        #                                   theta1, theta2, theta3, theta12, theta13, theta23,
        #                                   beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        # x_all[,1,,] <- list_xall$x_all1
        # x_all[,2,,] <- list_xall$x_all2
        # x_all[,3,,] <- list_xall$x_all3
        # N_all <- list_xall$N_all
        
        list_xall <- update_x_all_foreach(x_all, N_all, 10,
                                          theta1, theta2, theta3, theta12, theta13, theta23,
                                          beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
        
        # if(iter > nburn & ((iter - nburn) %% nthin == 0)){
          acceptances_iter[iter,c(1,2,3)] <- 1  
        # }
        
      }

    }

    # UPDATE THETA ---------------------------------------------------------
    
    theta1_star <- rnorm(1, theta1, sd = .0005)
    theta2_star <- rnorm(1, theta2, sd = .0005)
    theta3_star <- rnorm(1, theta3, sd = .0005)
    
    ratio <- hastings_ratio(s1, s2, s3, 
                            N1, N2, N3, 
                            x_all, N_all,
                            beta1, beta2, beta3, 
                            beta1, beta2, beta3, 
                            theta1, theta2, theta3,
                            theta1_star, theta2_star, theta3_star,
                            theta12, theta13, theta23,
                            theta12, theta13, theta23,
                            a_theta, b_theta)
    
    if(theta1_star > 0 & theta2_star > 0 & theta3_star > 0){
      
      if(runif(1) < ratio){
        
        theta1 <- theta1_star
        theta2 <- theta2_star
        theta3 <- theta3_star
        
        # x_all1 <- x_all[,1,,]
        # x_all2 <- x_all[,2,,]
        # x_all3 <- x_all[,3,,]
        # 
        # list_xall <- update_x_all_tri_cpp(x_all1, x_all2,  x_all3, N_all, 10,
        #                               theta1, theta2, theta3, theta12, theta13, theta23,
        #                               beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        # x_all[,1,,] <- list_xall$x_all1
        # x_all[,2,,] <- list_xall$x_all2
        # x_all[,3,,] <- list_xall$x_all3
        # N_all <- list_xall$N_all
        
        list_xall <- update_x_all_foreach(x_all, N_all, 10,
                                          theta1, theta2, theta3, theta12, theta13, theta23,
                                          beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
        
        acceptances_iter[iter,c(4,5,6)] <- 1
      }
      
    }
    
    # UPDATE THETA ---------------------------------------------------------
    
    theta12_star <- rnorm(1, theta12, sd = .0005)
    theta13_star <- rnorm(1, theta13, sd = .0005)
    theta23_star <- rnorm(1, theta23, sd = .0005)
    
    ratio <- hastings_ratio(s1, s2, s3, 
                            N1, N2, N3, 
                            x_all, N_all,
                            beta1, beta2, beta3, 
                            beta1, beta2, beta3, 
                            theta1, theta2, theta3,
                            theta1, theta2, theta3,
                            theta12, theta13, theta23,
                            theta12_star, theta13_star, theta23_star,
                            a_theta, b_theta)
    
    if(theta12_star > 0 & theta13_star > 0 & theta23_star > 0){
      
      if(runif(1) < ratio){
        
        theta12 <- theta12_star
        theta13 <- theta13_star
        theta23 <- theta23_star
        
        # x_all1 <- x_all[,1,,]
        # x_all2 <- x_all[,2,,]
        # x_all3 <- x_all[,3,,]
        # 
        # list_xall <- update_x_all_tri_cpp(x_all1, x_all2,  x_all3, N_all, 10,
        #                               theta1, theta2, theta3, theta12, theta13, theta23,
        #                               beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        # x_all[,1,,] <- list_xall$x_all1
        # x_all[,2,,] <- list_xall$x_all2
        # x_all[,3,,] <- list_xall$x_all3
        # N_all <- list_xall$N_all
        
        list_xall <- update_x_all_foreach(x_all, N_all, 10,
                                          theta1, theta2, theta3, theta12, theta13, theta23,
                                          beta1, beta2, beta3, diag(0.05, nrow = 2), a, b)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
        
        acceptances_iter[iter,c(7,8,9)] <- 1
      }
      
    }
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn)/nthin
      params_iter[trueIter,] <- c(theta1, theta2, theta3, theta12, theta13, theta23, beta1, beta2, beta3)
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
