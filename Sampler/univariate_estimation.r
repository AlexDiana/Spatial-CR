library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2)
# setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
setwd("C:/Users/ad603/Dropbox/R Folder/PhD/Spatial")
sourceCpp("rcpp/code_basic.cpp")

# FUNCTIONS ---------------------------------------------------------------

log_prior <- function(beta, theta, a_theta, b_theta, a_beta, b_beta){
  
  dgamma(theta, a_theta, b_theta, log = T) + dgamma(beta, a_beta, b_beta, log = T)
  
}

hastings_ratio <- function(data, N,  
                           x_all, N_all,
                           beta, beta_star, 
                           theta, theta_star,
                           a_theta, b_theta,
                           a_beta, b_beta){
  
  ratio_constants <- importance_sampling_estimate(x_all, N_all,
                                                  beta, beta_star, 
                                                  theta, theta_star)
  
  den <- log_prior(beta, theta, a_theta, b_theta, a_beta, b_beta) + log_f_softcore_cpp(data, N, beta, theta)
  
  num <- log_prior(beta_star, theta_star, a_theta, b_theta, a_beta, b_beta) + log_f_softcore_cpp(data, N, beta_star, theta_star)
  
  exp(num - den) * ratio_constants
}

hastings_ratio_PS <- function(data, N,  
                              grid_theta, grid_beta,
                              grid_pathsampling,
                              beta, beta_star, 
                              theta, theta_star,
                              a_theta, b_theta,
                              a_beta, b_beta){
  
  ratio_constants <- path_sampling_estimate(beta, beta_star, 
                                            theta, theta_star,
                                            grid_theta, grid_beta,
                                            grid_pathsampling)
  
  den <- log_prior(beta, theta, a_theta, b_theta, a_beta, b_beta) + log_f_softcore_cpp(data, N, beta, theta)
  
  num <- log_prior(beta_star, theta_star, a_theta, b_theta, a_beta, b_beta) + log_f_softcore_cpp(data, N, beta_star, theta_star)
  
  exp(num - den) * ratio_constants
}

# compute c_theta / c_theta0
importance_sampling_estimate <- function(x_all, N_all,
                                         beta, beta_star, 
                                         theta, theta_star){
  
  constant <- 0
  M <- dim(x_all)[1]
  output <- rep(NA, M)
  for (m in 1:M) {
    
    N <- N_all[m]
    x <- matrix(x_all[m,seq_len(N),], nrow = N)
    
    num <- log_f_softcore_cpp(x, N, beta, theta)
    
    den <- log_f_softcore_cpp(x, N, beta_star, theta_star)
    
    constant <- constant + exp(num - den)
    output[m] <- exp(num - den)
  }
  
  constant / M
}

# compute c_theta / c_theta0
importance_sampling_estimate_new <- function(x, N,
                                         beta, beta_star, 
                                         theta, theta_star, niter){
  
  list_simsip <- simulate_softcore_importancesampling(x, N, theta, theta_star, beta, niter = 200000, 
                                                      Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta,
                                                      Sigma_newpoint = diag(1, nrow = 2))
  likelihood_output1 <- list_simsip$likelihood1
  likelihood_output2 <- list_simsip$likelihood2
  
  # thin
  # cumsumThinned <- 
  # qplot(1:200000, cumsum(exp(likelihood_output1 - likelihood_output2))/(1:200000))
  
  mean(exp(likelihood_output1 - likelihood_output2))
}

path_sampling_estimate <- function(beta, beta_star, 
                                   theta, theta_star,
                                   grid_theta, grid_beta,
                                   grid_pathsampling){
  
  # approximate theta and theta_star to closest in grid
  index_theta <- which.min(abs(grid_theta - theta))
  index_theta_star <- which.min(abs(grid_theta - theta_star))
  
  theta_grid <- grid_theta[index_theta]
  theta_star_grid <- grid_theta[index_theta_star]
  
  # approximate beta and beta_star to closest in grid
  index_beta <- which.min(abs(grid_beta - beta))
  index_beta_star <- which.min(abs(grid_beta - beta_star))
  
  beta_grid <- grid_beta[index_beta]
  beta_star_grid <- grid_beta[index_beta_star]
  
  current_grid_theta <- grid_pathsampling[index_beta,seq(index_theta, index_theta_star),1]
  nsteps_theta <- length(current_grid_theta) - 1
  
  current_grid_beta <- grid_pathsampling[seq(index_beta, index_beta_star),index_theta_star,1]
  nsteps_beta <- length(current_grid_beta) - 1
  
  if(nsteps_theta > 0){
    logctheta <- (current_grid_theta[1] / 2) + sum(current_grid_theta[2:(nsteps_theta - 1)]) + 
      (current_grid_theta[nsteps_theta] / 2)
    logctheta <- logctheta * (theta - theta_star)
    
    return(exp(logctheta / (nsteps_theta - 1)))
  } else {
    return(1)
  }

}

generateEVTheta <- function(beta_min, beta_max, step_beta, 
                theta_min, theta_max, step_theta, M){
  
  grid_beta <- seq(beta_min, beta_max, by = step_beta)
  grid_theta <- seq(theta_min, theta_max, by = step_theta)
  
  EVTheta <- array(NA, dim = c(length(grid_beta), length(grid_theta), 2))
  
  for (i in seq_along(grid_beta)) {
    
    for (j in seq_along(grid_theta)) {
      
      print(paste0("i = ",i," / j = ",j))
      
      EVTheta[i,j,] <- estimateEVTheta_cpp(grid_beta[i], grid_theta[j], 2000, 
                                           Sigma_prop = diag(0.01, nrow = 2), 
                                           Sigma_newpoint = diag(1, nrow = 2), 500, 300)
      
    } 
     
  }
  
  EVTheta
}

estimateEVTheta <- function(beta, theta, M){
  
  VTheta <- matrix(NA, nrow = M, ncol = 2)
  
  for (m in 1:M) {
    
    # print(m)
    
    list_s <- simulate_softcore_cpp(theta, beta, 2000, 
                                    Sigma_prop = diag(0.01, nrow = 2), 500, 100, 
                                    Sigma_newpoint = diag(1, nrow = 2))
    s <- list_s$data
    N <- list_s$N
    
    s <- s[1:N,,drop = F]
    
    VTheta[m,1] <- dlog_dtheta_softcore_cpp(s, N, theta)   
    VTheta[m,2] <- dlog_dbeta_softcore_cpp(N, beta)   
    
  }
  
  EVTheta <- apply(VTheta, 2, mean)
  
  EVTheta
}

# GENERATE PATH SAMPLING MATRIX -------------------------------------------

beta_min <- 150
beta_max <- 150
length_beta <- 1
step_beta <- ifelse(length_beta > 1,(beta_max - beta_min) / (length_beta - 1), 0)
grid_beta <- seq(beta_min, beta_max, by = step_beta)

theta_min <- .01
theta_max <- .0001
length_theta <- 301
# step_theta <- ifelse(length_theta > 1,(theta_max - theta_min) / (length_theta - 1),0)
grid_theta <- exp(seq(log(theta_min), log(theta_max), length.out = length_theta))

length(grid_theta) * length(grid_beta)

grid_pathsampling <- generateEVTheta_cpp(grid_beta, grid_theta, M = 2000, 
                                         Sigma_prop = diag(0.01, nrow = 2), 
                                         Sigma_newpoint = diag(1, nrow = 2),
                                         niter = 300, Nmax = 500) 

qplot(1:length_theta, grid_pathsampling[1,,1])

# COMPARE ESTIMATES -------------------------------------------------------

theta <- .0002
theta_star <- .00024
beta <- 150
beta_star <- 150

path_sampling_estimate(beta, beta, 
                       theta, theta_star,
                       grid_theta, grid_beta,
                       grid_pathsampling)

{
  M <- 100
  niter <- 10000

  x_all <- array(NA, dim = c(M, 500, 2))
  N_all <- rep(NA, M)
  
  likelihood_output <- matrix(NA, nrow = M, ncol = niter)

  for (m in 1:M) {

    print(m)

    list_sims <- simulate_softcore_cpp(theta, beta, niter = niter,
                                       Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta,
                                       Sigma_newpoint = diag(1, nrow = 2))
    N_sim <- list_sims$N
    x_all[m,1:N_sim,] <- list_sims$data[1:N_sim,]
    N_all[m] <- N_sim
    
    likelihood_output[m,] <- list_sims$likelihood

    qplot(1:niter, likelihood_output[m,])
    
  }
  
  
   likelihood_output

  importance_sampling_estimate(x_all, N_all,
                               beta, beta,
                               theta, theta_star)

  #
  
  list_sims <- simulate_softcore_cpp(theta, beta, niter = 100000,
                                     Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta,
                                     Sigma_newpoint = diag(1, nrow = 2))
  N <- list_sims$N
  x <- list_sims$data
  
  importance_sampling_estimate_new(x, N,
                                   beta, beta_star, 
                                   theta, theta_star, 50000)
}

# SIMULATING DATA ---------------------------------------------------------

usingTrueData <- F

theta <- .001
beta <- 150

beta_true <- beta
theta_true <- theta

list_s <- simulate_softcore_cpp(theta, beta, 2000, 
                                      Sigma_prop = diag(0.01, nrow = 2), 500, 100, 
                                      Sigma_newpoint = diag(1, nrow = 2))
s <- list_s$data
N <- list_s$N

s <- s[1:N,]
print(N)
ggplot() +  geom_point(data = NULL, aes(x = s[,1], y = s[,2]), color = "red") 

# PRIOR -------------------------------------------------------------------

Nmax <- 300

meanBeta <- 200
varBeta <- 4000
b_beta <- meanBeta / varBeta
a_beta <- meanBeta * b_beta

a_theta <- 1.5
b_theta <- 20

sigma_prop <- .01
Sigma_newpoint <- diag(1, nrow = 2)

ggplot(data = NULL, aes(x = c(0,1))) + stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))

# MCMC PATH SAMPLING --------------------------------------------------------------------

nburn <- 10000
niter <- 10000
nchain <- 1
nthin <- 5

# variables for output
{
  
}

for(chain in 1:nchain) {
  
  # starting values
  {
    
    # point process parameters
    {
      beta <- beta_true
      theta <- theta_true 
    }
    
  }
  
  # output variables
  {
    params_iter <- matrix(NA, nrow = niter, ncol = 2)
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter > nburn){
      print(iter)
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print((iter - nburn)/nthin) 
      print(paste0("Theta = ",theta," / Beta = ",beta))
    }
    
    # UPDATE THETA ---------------------------------------------------------
    
    theta_star <- rnorm(1, theta, sd = .0001)
    
    ratio <- hastings_ratio_PS(s, N,  
                               grid_theta, grid_beta,
                               grid_pathsampling,
                               beta, beta_star, 
                               theta, theta_star,
                               a_theta, b_theta,
                               a_beta, b_beta)
    
    if(theta_star > 0){
      
      if(runif(1) < ratio){
        
        theta <- theta_star
        
      }
      
    }
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn)/nthin
      params_iter[trueIter,] <- c(theta, beta)
    }
    
  }
  
  # Write results in MCMC output
  if(iter > nburn & (iter - nburn) %% nthin == 0) {
    
    
    
  }
  
}

summary(params_iter[,1])
qplot(params_iter[,1]) + geom_histogram()

# MCMC IMPORTANCE SAMPLING --------------------------------------------------------------------

nburn <- 5000
niter <- 5000
nchain <- 1
nthin <- 5

# variables for output
{
  
}

for(chain in 1:nchain) {
  
  # starting values
  {
    
    # point process parameters
    {
      beta <- 150
      theta <- theta_true 
    }
    
    # additional random variables
    {
      M <- 100
      
      x_all <- array(NA, dim = c(M, Nmax, 2))
      N_all <- rep(NA, M)
      
      for (m in 1:M) {
        
        print(m)
        
        list_sims <- simulate_softcore_cpp(theta, beta, niter = 10000, 
                                           Sigma_prop = diag(.01, nrow = 2), Nmax = 500, lambda = beta, 
                                           Sigma_newpoint = diag(1, nrow = 2))
        N_sim <- list_sims$N
        x_all[m,1:N_sim,] <- list_sims$data[1:N_sim,]
        N_all[m] <- N_sim
        
      }
      
    }
    
  }
  
  # output variables
  {
    params_iter <- matrix(NA, nrow = niter, ncol = 2)
  }
  
  #iterations
  for(iter in seq_len(nburn + nthin*niter)){
    
    if(iter < nburn){
      print(iter) 
      print(paste0("Theta = ",theta," / Beta = ",beta))
    } else if(((iter - nburn)/nthin) %% 10 == 0){
      print((iter - nburn)/nthin) 
      print(paste0("Theta = ",theta," / Beta = ",beta))
    }
    
    # UPDATE THETA ---------------------------------------------------------
    
    theta_star <- rnorm(1, theta, sd = .0001)
    
    ratio <- hastings_ratio(s, N,
                            x_all, N_all,
                            beta, beta,
                            theta, theta_star,
                            a_theta, b_theta,
                            a_beta, b_beta)
    
    if(theta_star > 0){
      
      if(runif(1) < ratio){
        
        theta <- theta_star
        
        list_xall <- update_x_all_uni_cpp(x_all, N_all, 200,
                                      theta, beta, diag(0.01, nrow = 2), Sigma_newpoint)
        x_all <- list_xall$x_all
        N_all <- list_xall$N_all
      }
      
    }
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn & (iter - nburn) %% nthin == 0){
      trueIter <- (iter - nburn)/nthin
      params_iter[trueIter,] <- c(theta, beta)
    }
    
  }
  
  # Write results in MCMC output
  if(iter > nburn & (iter - nburn) %% nthin == 0) {
    
    
    
  }
  
}

summary(params_iter[,1])
qplot(params_iter[,1]) + geom_histogram()

