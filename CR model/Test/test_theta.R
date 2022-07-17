theta1 <- exp(-5)
theta2 <- exp(-10)
theta12 <- exp(-5)
beta1 <- 500
beta2 <- 500

beta1_true <- beta1
beta2_true <- beta2
theta1_true <- theta1
theta2_true <- theta2
theta12_true <- theta12

niter <- 20000

list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, niter, 
                                      Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                      Sigma_newpoint = diag(1, nrow = 2), allTraps, R_traps,
                                      a1, b1, a2, b2)
s1 <- list_s1s2$data1
s2 <- list_s1s2$data2
N1 <- list_s1s2$N1
N2 <- list_s1s2$N2
likelihoods <- list_s1s2$loglikelihoods

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

# print(N1)
# print(N2)
# ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
#   geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
# 
# qplot(1:niter, likelihoods)

# PROFILE LIKELIHOOD -------

logtheta2_grid <- seq(log(theta2_true) - 3, log(theta2_true) + .5, length.out = 30)
fun_grid <- seq_len(length(theta2_grid))

M <- 150

for (i in seq_along(logtheta2_grid)) {
  
  print(i)
  
  if(i == 1){
    
    x_all <- array(NA, dim = c(M, 2, Nmax, 2))
    N_all <- matrix(NA, nrow = M, ncol = 2)
    
    for (m in 1:M) {
      
      # print(m)
      
      list_sims <- simulate_bivsoftcore_cpp(theta1, exp(logtheta2_grid[1]), theta12, beta1, beta2,
                                            niter = 5000,
                                            Sigma_prop = diag(.01, nrow = 2), Nmax = Nmax, lambda = beta1,
                                            Sigma_newpoint = diag(1, nrow = 2),
                                            allTraps, R_traps,
                                            a1, b1, a2, b2)
      N1_sim <- list_sims$N1
      N2_sim <- list_sims$N2
      x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
      x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
      N_all[m,] <- c(N1_sim, N2_sim)
      
    }
    
  } else {
    
    x_all1 <- x_all[,1,,]
    x_all2 <- x_all[,2,,]
    
    list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 1000,
                                  theta1, exp(logtheta2_grid[i]), theta12,
                                  beta1, beta2, diag(0.05, nrow = 2),
                                  Sigma_newpoint = diag(1, nrow = 2),
                                  allTraps, R_traps,
                                  a1, b1, a2, b2)
    x_all[,1,,] <-  list_xall$x_all1
    x_all[,2,,] <-  list_xall$x_all2
    N_all <- list_xall$N_all
    
  }
  
  lognormConstants <- sapply(1:M, function(m){
    
    N1_sim <- N_all[m,1]
    N2_sim <- N_all[m,2]
    
    log_f_bivsoftcore_cpp(x_all[m,1,1:N1_sim,], 
                          x_all[m,2,1:N2_sim,], 
                          N_all[m,1], 
                          N_all[m,2],
                          beta1, beta2,
                          theta1, exp(logtheta2_grid[i]), theta12)
  })
  
  maxNormConst <- max(lognormConstants)
  
  value1 <- log_f_bivsoftcore_cpp(s1, s2, N1, N2,
                                  beta1, beta2,
                                  theta1, exp(logtheta2_grid[i]), theta12)
  
  fun_grid[i] <- exp(value1 - maxNormConst)  / mean(exp(lognormConstants - maxNormConst))
  
}

qplot(logtheta2_grid, fun_grid) + 
  geom_vline(aes(xintercept = log(theta2_true)))

# ESTIMATE THETA

# prior
{
  Nmax <- 1000
  
  meanBeta <- 500
  varBeta <- 10000
  b_beta <- meanBeta / varBeta
  a_beta <- meanBeta * b_beta  
  
  ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,1000))
  
  a_theta <- 3
  b_theta <- 3
  
  epsilon_beta <- 5000
  epsilon_logtheta1 <- .75
  epsilon_logtheta2 <- .75
  epsilon_logtheta12 <- .25
  }

niter <- 1000

# starting values
{
  beta1 <- 600
  beta2 <- 600
  theta1 <- exp(-8)
  theta2 <- exp(-8)
  theta12 <- exp(-8)
  
  M <- 40
  
  x_all <- array(NA, dim = c(M, 2, Nmax, 2))
  N_all <- matrix(NA, nrow = M, ncol = 2)
  
  for (m in 1:M) {
    
    print(m)
    
    list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                          niter = 1500,
                                          Sigma_prop = diag(.01, nrow = 2), Nmax = Nmax, lambda = beta1,
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          allTraps, R_traps,
                                          a1, b1, a2, b2)
    N1_sim <- list_sims$N1
    N2_sim <- list_sims$N2
    x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
    x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
    N_all[m,] <- c(N1_sim, N2_sim)
    
  }
}

#  output
{
  params_output <- matrix(NA, nrow = niter, ncol = 5)
}

for (iter in 1:niter) {
  
  print(iter)
  
  # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
  
  # if(iter > iterAfterAdapting){
  #   Sigma_n <- cov(params_values[1:(iter-1),c(1,3,5)])
  #   Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
  #     beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
  # } else {
  Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta1, epsilon_logtheta12), nrow = 3) 
  # }
  
  proposed_theta <- mvrnorm(1, c(beta1, log(theta1), log(theta12)),
                            Sigma = Sigma_proposal)
  beta1_star <- proposed_theta[1]
  theta1_star <- exp(proposed_theta[2])
  theta12_star <- exp(proposed_theta[3])
  
  ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                              x_all[,1,,], 
                              x_all[,2,,],
                              c(beta1, beta2,
                                theta1, theta2, theta12),
                              c(beta1_star, beta2,
                                theta1_star, theta2, theta12_star),
                              N_all,
                              mu_logtheta, sd_logtheta,
                              # a_theta, b_theta,
                              a_beta, b_beta)
  
  # importanceSamplingEstimateParallel(x_all[,1,,], x_all[,2,,], 
  #                                    c(beta1, beta2,
  #                                      theta1, theta2, theta12),
  #                                    c(beta1_star, beta2,
  #                                      theta1_star, theta2, theta12_star),
  #                                    N_all)
  # 
  # x_all1 <- x_all[,1,,] 
  # x_all2 <- x_all[,2,,] 
  # 
  # sapply(1:M, function(m){
  #   
  #   N1 <- N_all[m, 1]
  #   N2 <- N_all[m, 2]
  #   
  #   x1_current = x_all1[m,1:N1,]
  #   x2_current = x_all2[m,1:N2,]
  #   
  #   num_minus_den <- log_f_bivsoftcore_cpp_ratio(x1_current, x2_current, N1, N2, 
  #                                                beta1, beta2,
  #                                                theta1, theta2, theta12,
  #                                                beta1_star, beta2_star,
  #                                                theta1_star, theta2_star, theta12_star);
  #   
  #   num_minus_den
  # })
  
  print(ratio)
  # ratio <- hastings_ratio(s1, s2, N1, N2,
  #                         x_all, N_all,
  #                         beta1, beta2,
  #                         beta1_star, beta2,
  #                         theta1, theta2, theta12,
  #                         theta1_star, theta2, theta12_star,
  #                         # mu_gamma, sigmasq_gamma, 
  #                         a_theta, b_theta,
  #                         a_beta, b_beta)
  # print(ratio)
  if(!is.na(ratio) & theta1_star > 0){
    
    if(runif(1) < ratio){
      beta1 <- beta1_star
      theta1 <- theta1_star
      theta12 <- theta12_star
      
      x_all1 <- x_all[,1,,]
      x_all2 <- x_all[,2,,]
      
      list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 100,
                                    theta1, theta2, theta12,
                                    beta1, beta2, diag(0.05, nrow = 2),
                                    Sigma_newpoint = diag(1, nrow = 2),
                                    allTraps, R_traps,
                                    a1, b1, a2, b2)
      x_all[,1,,] <-  list_xall$x_all1
      x_all[,2,,] <-  list_xall$x_all2
      N_all <- list_xall$N_all
      
      # list_xall <- update_x_all_foreach(x_all, N_all, 50,
      #                                   theta1, theta2, theta12,
      #                                   beta1, beta2, diag(0.05, nrow = 2),
      #                                   Sigma_newpoint = diag(1, nrow = 2),
      #                                   allTraps, R_traps,
      #                                   a1, b1, a2, b2)
      # x_all <- list_xall$x_all
      # N_all <- list_xall$N_all
      
      print("accepted 1")
    }
    
  }
  
  # params_values[iter,c(1,3,5)] <- c(beta1, log(theta1), log(theta12))
  
  # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
  
  # if(iter > iterAfterAdapting){
  #   Sigma_n <- cov(params_values[1:(iter-1),c(2,4,5)])
  #   Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
  #     beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
  # } else {
  Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta2, epsilon_logtheta12), nrow = 3) 
  # }
  
  proposed_theta <- mvrnorm(1, c(beta2, log(theta2), log(theta12)),
                            Sigma = Sigma_proposal)
  beta2_star <- proposed_theta[1]
  theta2_star <- exp(proposed_theta[2])
  theta12_star <- exp(proposed_theta[3])
  
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
                              a_beta, b_beta)
  
  # ratio <- hastings_ratio(s1, s2, N1, N2,
  #                         x_all, N_all,
  #                         beta1, beta2,
  #                         beta1, beta2_star,
  #                         theta1, theta2, theta12,
  #                         theta1, theta2_star, theta12_star,
  #                         # mu_gamma, sigmasq_gamma, 
  #                         a_theta, b_theta,
  #                         a_beta, b_beta)
  # print(ratio)
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
                                    allTraps, R_traps,
                                    a1, b1, a2, b2)
      x_all[,1,,] <-  list_xall$x_all1
      x_all[,2,,] <-  list_xall$x_all2
      N_all <- list_xall$N_all
      
      # list_xall <- update_x_all_foreach(x_all, N_all, 50,
      #                                   theta1, theta2, theta12,
      #                                   beta1, beta2, diag(0.05, nrow = 2),
      #                                   Sigma_newpoint = diag(1, nrow = 2),
      #                                   allTraps, R_traps,
      #                                   a1, b1, a2, b2)
      # x_all <- list_xall$x_all
      # N_all <- list_xall$N_all
      
      print("accepted 2")
      # acceptances_theta[iter,1] <- 1
    }
    
  }
  
  # params_values[iter,c(2,4,5)] <- c(beta2, log(theta2), log(theta12))
  params_output[iter,] <- c(beta1, beta2, log(theta1), log(theta2), log(theta12))
  
}

qplot(1:niter, params_output[,1]) + geom_hline(aes(yintercept = beta1_true))
qplot(1:niter, params_output[,2]) + geom_hline(aes(yintercept = beta2_true))
qplot(1:niter, params_output[,3]) + geom_hline(aes(yintercept = log(theta1_true))) + ylim(c(-20, -2))
qplot(1:niter, params_output[,4]) + geom_hline(aes(yintercept = log(theta2_true))) + ylim(c(-20, -2))
qplot(1:niter, params_output[,5]) + geom_hline(aes(yintercept = log(theta12_true)))  + ylim(c(-20, -2))


# MH  -------

M <- 150

# initial values
{
  x_all <- array(NA, dim = c(M, 2, Nmax, 2))
  N_all <- matrix(NA, nrow = M, ncol = 2)
  
  for (m in 1:M) {
    
    # print(m)
    
    list_sims <- simulate_bivsoftcore_cpp(theta1, theta2_true, theta12, beta1, beta2,
                                          niter = 5000,
                                          Sigma_prop = diag(.01, nrow = 2), Nmax = Nmax, lambda = beta1,
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          allTraps, R_traps,
                                          a1, b1, a2, b2)
    N1_sim <- list_sims$N1
    N2_sim <- list_sims$N2
    x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
    x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
    N_all[m,] <- c(N1_sim, N2_sim)
    
  }
}

logtheta2 <- log(theta2_true)

niter <- 2000

logtheta2_output <- rep(NA, niter)

for (iter in 1:niter) {
  
  print(iter)
  
  # logtheta2_star <- rnorm(1, logtheta2, sd = 1)
  
  theta2 <- exp(logtheta2)
  theta2_star <- rnorm(1, theta2, sd = .005)
  if(theta2_star > 0){
    logtheta2_star <- log(theta2_star)  
    
    
    ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                x_all[,1,,], 
                                x_all[,2,,],
                                c(beta1, beta2,
                                  theta1, exp(logtheta2), theta12),
                                c(beta1, beta2,
                                  theta1, exp(logtheta2_star), theta12),
                                N_all,
                                mu_logtheta = 1, sd_logtheta = 100,
                                # a_theta, b_theta,
                                a_beta = 1, b_beta= 1)
    
    print(paste0("mh_ratio = ",ratio))
    
    if(!is.na(ratio)){
      
      if(runif(1) < ratio){
        
        logtheta2 <- logtheta2_star
        
        x_all1 <- x_all[,1,,]
        x_all2 <- x_all[,2,,]
        
        list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 200,
                                      theta1, exp(logtheta2), theta12,
                                      beta1, beta2, diag(0.05, nrow = 2),
                                      Sigma_newpoint = diag(1, nrow = 2),
                                      allTraps, R_traps,
                                      a1, b1, a2, b2)
        x_all[,1,,] <-  list_xall$x_all1
        x_all[,2,,] <-  list_xall$x_all2
        N_all <- list_xall$N_all
        
      }
      
    }
    
  }
  
  logtheta2_output[iter] <- logtheta2
}

qplot(logtheta2_output) + 
  geom_vline(aes(xintercept = log(theta2_true)))

qplot(1:niter, logtheta2_output) + 
  geom_hline(aes(yintercept = log(theta2_true)))
# qplot(log(-logtheta2_output)) + 
# geom_vline(aes(xintercept = log(-log(theta2_true))))


# ESTIMATE THETA ---------

# prior
{
  Nmax <- 1000
  
  meanBeta <- 500
  varBeta <- 10000
  b_beta <- meanBeta / varBeta
  a_beta <- meanBeta * b_beta  
  
  ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,1000))
  
  a_theta <- 3
  b_theta <- 3
  
  # epsilon_beta <- 1000
  # epsilon_logtheta1 <- .5
  # epsilon_logtheta2 <- .5
  # epsilon_logtheta12 <- .5
  
  epsilon_beta <- 1000
  epsilon_theta1 <- .0005^2
  epsilon_theta2 <- .0005^2
  epsilon_theta12 <- .0005^2
  
  # adaptation parameters for proposal of interaction parameters
  {
    beta_proposal <- .1
    iterAfterAdapting <- 1000
  }
}

niter <- 5000

# starting values
{
  beta1 <- 600
  beta2 <- 600
  theta1 <- exp(-8)
  theta2 <- exp(-8)
  theta12 <- exp(-8)
  
  M <- 40
  
  x_all <- array(NA, dim = c(M, 2, Nmax, 2))
  N_all <- matrix(NA, nrow = M, ncol = 2)
  
  for (m in 1:M) {
    
    print(m)
    
    list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                          niter = 1500,
                                          Sigma_prop = diag(.01, nrow = 2), Nmax = Nmax, lambda = beta1,
                                          Sigma_newpoint = diag(1, nrow = 2),
                                          allTraps, R_traps,
                                          a1, b1, a2, b2)
    N1_sim <- list_sims$N1
    N2_sim <- list_sims$N2
    x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
    x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
    N_all[m,] <- c(N1_sim, N2_sim)
    
  }
  
  params_values <- matrix(NA, nrow = niter, ncol = 5)
}

# output
{
  params_output <- matrix(NA, nrow = niter, ncol = 5)
}

for (iter in 1:niter) {
  
  print(iter)
  
  # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
  
  # if(iter > iterAfterAdapting){
  #   Sigma_n <- cov(params_values[1:(iter-1),c(1,3,5)])
  #   Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
  #     beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
  # } else {
  # Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta1, epsilon_logtheta12), nrow = 3)
  # }
  # 
  # proposed_theta <- mvrnorm(1, c(beta1, log(theta1), log(theta12)),
  #                           Sigma = Sigma_proposal)
  # beta1_star <- proposed_theta[1]
  # theta1_star <- exp(proposed_theta[2])
  # theta12_star <- exp(proposed_theta[3])
  
  if(iter > iterAfterAdapting){
    Sigma_n <- cov(params_values[100:(iter-1),c(1,3,5)])
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
  
  # theta1_star > 0
  # theta12_star > 0
  
  if(theta1_star > 0 & theta12_star > 0){
    
    print("in 1")
    
    ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                x_all[,1,,], 
                                x_all[,2,,],
                                c(beta1, beta2,
                                  theta1, theta2, theta12),
                                c(beta1_star, beta2,
                                  theta1_star, theta2, theta12_star),
                                N_all,
                                mu_logtheta = 1, sd_logtheta = 100,
                                # a_theta, b_theta,
                                a_beta, b_beta)
    
    # importanceSamplingEstimateParallel(x_all[,1,,], x_all[,2,,],
    #                                    c(beta1, beta2,
    #                                      theta1, theta2, theta12),
    #                                    c(beta1_star, beta2,
    #                                      theta1_star, theta2, theta12_star),
    #                                    N_all)
    # 
    # x_all1 <- x_all[,1,,] 
    # x_all2 <- x_all[,2,,] 
    # 
    # sapply(1:M, function(m){
    #   
    #   N1 <- N_all[m, 1]
    #   N2 <- N_all[m, 2]
    #   
    #   x1_current = x_all1[m,1:N1,]
    #   x2_current = x_all2[m,1:N2,]
    #   
    #   num_minus_den <- log_f_bivsoftcore_cpp_ratio(x1_current, x2_current, N1, N2, 
    #                                                beta1, beta2,
    #                                                theta1, theta2, theta12,
    #                                                beta1_star, beta2_star,
    #                                                theta1_star, theta2_star, theta12_star);
    #   
    #   num_minus_den
    # })
    
    print(ratio)
    # ratio <- hastings_ratio(s1, s2, N1, N2,
    #                         x_all, N_all,
    #                         beta1, beta2,
    #                         beta1_star, beta2,
    #                         theta1, theta2, theta12,
    #                         theta1_star, theta2, theta12_star,
    #                         # mu_gamma, sigmasq_gamma, 
    #                         a_theta, b_theta,
    #                         a_beta, b_beta)
    # print(ratio)
    if(!is.na(ratio) & theta1_star > 0){
      
      if(runif(1) < ratio){
        beta1 <- beta1_star
        theta1 <- theta1_star
        theta12 <- theta12_star
        
        x_all1 <- x_all[,1,,]
        x_all2 <- x_all[,2,,]
        
        list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 100,
                                      theta1, theta2, theta12,
                                      beta1, beta2, diag(0.05, nrow = 2),
                                      Sigma_newpoint = diag(1, nrow = 2),
                                      allTraps, R_traps,
                                      a1, b1, a2, b2)
        x_all[,1,,] <-  list_xall$x_all1
        x_all[,2,,] <-  list_xall$x_all2
        N_all <- list_xall$N_all
        
        # list_xall <- update_x_all_foreach(x_all, N_all, 50,
        #                                   theta1, theta2, theta12,
        #                                   beta1, beta2, diag(0.05, nrow = 2),
        #                                   Sigma_newpoint = diag(1, nrow = 2),
        #                                   allTraps, R_traps,
        #                                   a1, b1, a2, b2)
        # x_all <- list_xall$x_all
        # N_all <- list_xall$N_all
        
        print("accepted 1")
      }
      
    }
    
  }
  
  # params_values[iter,c(1,3,5)] <- c(beta1, log(theta1), log(theta12))
  
  # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
  
  # # if(iter > iterAfterAdapting){
  # #   Sigma_n <- cov(params_values[1:(iter-1),c(2,4,5)])
  # #   Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
  # #     beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
  # # } else {
  # Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta2, epsilon_logtheta12), nrow = 3) 
  # # }
  # 
  # proposed_theta <- mvrnorm(1, c(beta2, log(theta2), log(theta12)),
  #                           Sigma = Sigma_proposal)
  # beta2_star <- proposed_theta[1]
  # theta2_star <- exp(proposed_theta[2])
  # theta12_star <- exp(proposed_theta[3])
  
  if(iter > iterAfterAdapting){
    Sigma_n <- cov(params_values[100:(iter-1),c(2,4,5)])
    Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
      beta_proposal * diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
  } else {
    Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
  }
  
  # Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
  
  proposed_theta <- mvrnorm(1, c(beta2, theta2, theta12),
                            Sigma = Sigma_proposal)
  beta2_star <- proposed_theta[1]
  theta2_star <- proposed_theta[2]
  theta12_star <- proposed_theta[3]
  
  if(theta2_star > 0 & theta12_star > 0){
    
    print("in 2")
    
    ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                x_all[,1,,], 
                                x_all[,2,,],
                                c(beta1, beta2,
                                  theta1, theta2, theta12),
                                c(beta1, beta2_star,
                                  theta1, theta2_star, theta12_star),
                                N_all,
                                mu_logtheta = 1, sd_logtheta = 100,
                                # a_theta, b_theta,
                                a_beta, b_beta)
    
    # ratio <- hastings_ratio(s1, s2, N1, N2,
    #                         x_all, N_all,
    #                         beta1, beta2,
    #                         beta1, beta2_star,
    #                         theta1, theta2, theta12,
    #                         theta1, theta2_star, theta12_star,
    #                         # mu_gamma, sigmasq_gamma, 
    #                         a_theta, b_theta,
    #                         a_beta, b_beta)
    # print(ratio)
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
                                      allTraps, R_traps,
                                      a1, b1, a2, b2)
        x_all[,1,,] <-  list_xall$x_all1
        x_all[,2,,] <-  list_xall$x_all2
        N_all <- list_xall$N_all
        
        # list_xall <- update_x_all_foreach(x_all, N_all, 50,
        #                                   theta1, theta2, theta12,
        #                                   beta1, beta2, diag(0.05, nrow = 2),
        #                                   Sigma_newpoint = diag(1, nrow = 2),
        #                                   allTraps, R_traps,
        #                                   a1, b1, a2, b2)
        # x_all <- list_xall$x_all
        # N_all <- list_xall$N_all
        
        print("accepted 2")
        # acceptances_theta[iter,1] <- 1
      }
      
    }  
    
  }
  
  params_values[iter,] <- c(beta1, beta2, theta1, theta2, theta12)
  params_output[iter,] <- c(beta1, beta2, log(theta1), log(theta2), log(theta12))
  
}

qplot(1:niter, params_output[,1]) + geom_hline(aes(yintercept = beta1_true))
qplot(1:niter, params_output[,2]) + geom_hline(aes(yintercept = beta2_true))
qplot(1:niter, params_output[,3]) + geom_hline(aes(yintercept = log(theta1_true))) #+ ylim(c(-20, -2))
qplot(1:niter, params_output[,4]) + geom_hline(aes(yintercept = log(theta2_true))) # ylim(c(-20, -2))
qplot(1:niter, params_output[,5]) + geom_hline(aes(yintercept = log(theta12_true)))  # ylim(c(-20, -2))

# RNC ---------------------------------------------------------------------

beta1 <- 400
beta2 <- 400
theta1 <- exp(-8)
theta2 <- exp(-8)
theta12 <- exp(-8)

M <- 40

x_all <- array(NA, dim = c(M, 2, Nmax, 2))
N_all <- matrix(NA, nrow = M, ncol = 2)

for (m in 1:M) {
  
  print(m)
  
  list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                        niter = 1500,
                                        Sigma_prop = diag(.01, nrow = 2), Nmax = Nmax, lambda = beta1,
                                        Sigma_newpoint = diag(1, nrow = 2),
                                        allTraps, R_traps,
                                        a1, b1, a2, b2)
  N1_sim <- list_sims$N1
  N2_sim <- list_sims$N2
  x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
  x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
  N_all[m,] <- c(N1_sim, N2_sim)
  
}

beta1_star <- 410
beta2_star <- 410
theta1_star <- exp(-8.1)
theta2_star <- exp(-8.1)
theta12_star <- exp(-8.1)

importance_sampling_estimate_foreach(x_all, N_all,
                                     beta1, beta2,
                                     beta1_star, beta2_star, 
                                     theta1, theta2, theta12,
                                     theta1_star, theta2_star, theta12_star)

7827980167 # 200
2032811774 # 200
1549312132 # 400
1792531337 # 400
1191847312 # 400
1417133569 # 40
2128148623 # 40
