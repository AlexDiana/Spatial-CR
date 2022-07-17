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

logproposal_p0 <- function(p0, CH, mean_c_i){
  sum(CH) * log(p0) + sum(S - CH) * log(1 - p0 * mean_c_i)
}

update_p <- function(p0_1, sigma_1, p0_2, sigma_2, 
                     s1, s2, CH1, CH2, trapsX, trapsY, 
                     S_usage, N1, N2, 
                     sigma_0_1, sd_sigma_1,
                     sigma_0_2, sd_sigma_2,
                     sigma_p0_1_prop, sigma_p0_2_prop,
                     sigma_sigma1_prop, sigma_sigma2_prop, a_p, b_p){
  
  # update p0
  
  # p0_1_star <- rnorm(1, p0_1, sd = sigma_p0_1_prop)
  # 
  # if(p0_1_star > 0 & p0_1_star < 1){
  #   
  #   loglikelihood_p0_star <- loglikelihood_p_cpp(p0_1_star, sigma_1, CH1, s1,
  #                                                trapsX, trapsY, S_usage, N1)
  #   
  #   loglikelihood_p0 <- loglikelihood_p_cpp(p0_1, sigma_1, CH1, s1,
  #                                           trapsX, trapsY, S_usage, N1)
  #   
  #   logprior_p0_star <- dbeta(p0_1_star, a_p, b_p, log = T)
  #   
  #   logprior_p0 <- dbeta(p0_1, a_p, b_p, log = T)
  #   
  #   logposterior_p0_star <- logprior_p0_star + loglikelihood_p0_star
  #   
  #   logposterior_p0 <- logprior_p0 + loglikelihood_p0
  #   
  #   # logproposal_p0_current <- logproposal_p0(p0_1, CH1, mean(c_i))
  #   # 
  #   # logproposal_p0_star <- logproposal_p0(p0_1_star, CH1, mean(c_i))
  #   
  #   MH_ratio <- exp(logposterior_p0_star - logposterior_p0)
  #   
  #   if(runif(1) < MH_ratio){
  #     p0_1 <- p0_1_star
  #   }
  #   
  # }
  
  # update sigma
  
  # sigma_1_star <- rnorm(1, sigma_1, sd = sigma_sigma1_prop)
  # 
  # if(sigma_1_star > 0){
  #   
  #   loglikelihood_sigma_star <- loglikelihood_p_cpp(p0_1, sigma_1_star, CH1, s1,
  #                                                   trapsX, trapsY, S_usage, N1)
  #   
  #   loglikelihood_sigma <- loglikelihood_p_cpp(p0_1, sigma_1, CH1, s1,
  #                                              trapsX, trapsY, S_usage, N1)
  #   
  #   logprior_sigma_star <- dnorm(sigma_1_star, sigma_0_1, sd_sigma_1, log = T)
  #   
  #   logprior_sigma <- dnorm(sigma_1, sigma_0_1, sd_sigma_1, log = T)
  #   
  #   logposterior_sigma_star <- loglikelihood_sigma_star + logprior_sigma_star
  #   
  #   logposterior_sigma <- loglikelihood_sigma + logprior_sigma
  #   
  #   if(runif(1) < exp(loglikelihood_sigma_star - loglikelihood_sigma)){
  #     sigma_1 <- sigma_1_star
  #   }
  #   
  # }
  
  # update p0 and sigma
  
  p0_1_star <- rnorm(1, p0_1, sd = sigma_p0_1_prop)
  sigma_1_star <- rnorm(1, sigma_1, sd = sigma_sigma1_prop)
  
  if(p0_1_star > 0 & p0_1_star < 1 & sigma_1_star > 0){
    loglikelihood_p0_star <- loglikelihood_p_cpp(p0_1_star, sigma_1_star, CH1, s1,
                                                 trapsX, trapsY, S_usage, N1)
    
    loglikelihood_p0 <- loglikelihood_p_cpp(p0_1, sigma_1, CH1, s1,
                                            trapsX, trapsY, S_usage, N1)
    
    logprior_p0_star <- dbeta(p0_1_star, a_p, b_p, log = T)
    
    logprior_p0 <- dbeta(p0_1, a_p, b_p, log = T)
    
    logprior_sigma_star <- dnorm(sigma_1_star, sigma_0_1, sd_sigma_1, log = T)
    
    logprior_sigma <- dnorm(sigma_1, sigma_0_1, sd_sigma_1, log = T)
    
    logposterior_p0_star <- logprior_p0_star + loglikelihood_p0_star + logprior_sigma_star
    
    logposterior_p0 <- logprior_p0 + loglikelihood_p0 + logprior_sigma
    
    MH_ratio <- exp(logposterior_p0_star - logposterior_p0)
    
    if(runif(1) < MH_ratio){
      p0_1 <- p0_1_star
      sigma_1 <- sigma_1_star
    }
  }
  
  # update p0
  
  p0_2_star <- rnorm(1, p0_2, sd = sigma_p0_2_prop)
  
  if(p0_2_star > 0 & p0_2_star < 1){
    
    loglikelihood_p0_star <- loglikelihood_p_cpp2(p0_2_star, sigma_2, CH2, s2,
                                                 trapsX, trapsY, S_usage, N2)
    
    loglikelihood_p0 <- loglikelihood_p_cpp2(p0_2, sigma_2, CH2, s2,
                                            trapsX, trapsY, S_usage, N2)
    
    logprior_p0_star <- dbeta(p0_2_star, a_p, b_p, log = T)
    
    logprior_p0 <- dbeta(p0_2, a_p, b_p, log = T)
    
    # logproposal_p0_current <- logproposal_p0(p0_2, CH2, mean(c_i))
    # 
    # logproposal_p0_star <- logproposal_p0(p0_2_star, CH2, mean(c_i))
    
    logposterior_p0_star <- logprior_p0_star + loglikelihood_p0_star
    
    logposterior_p0 <- logprior_p0 + loglikelihood_p0
    
    MH_ratio <- exp(logposterior_p0_star - logposterior_p0)
    
    if(runif(1) < MH_ratio){
      p0_2 <- p0_2_star
    }
    
  }
  
  # update sigma
  
  sigma_2_star <- rnorm(1, sigma_2, sd = sigma_sigma2_prop)
  
  if(sigma_2_star > 0){
    
    loglikelihood_sigma_star <- loglikelihood_p_cpp2(p0_2, sigma_2_star, CH2, s2,
                                                    trapsX, trapsY, S_usage, N2)
    
    loglikelihood_sigma <- loglikelihood_p_cpp2(p0_2, sigma_2, CH2, s2,
                                               trapsX, trapsY, S_usage, N2)
    
    logprior_sigma_star <- dnorm(sigma_2_star, sigma_0_2, sd_sigma_2, log = T)
    
    logprior_sigma <- dnorm(sigma_2, sigma_0_2, sd_sigma_2, log = T)
    
    logposterior_sigma_star <- loglikelihood_sigma_star + logprior_sigma_star
    
    logposterior_sigma <- loglikelihood_sigma + logprior_sigma
    
    if(runif(1) < exp(loglikelihood_sigma_star - loglikelihood_sigma)){
      sigma_2 <- sigma_2_star
    }
    
  }
  
  list("p0_1" = p0_1, "sigma_1" = sigma_1,
       "p0_2" = p0_2, "sigma_2" = sigma_2)
}

update_p_laplace <- function(p0, sigma, 
                             s, CH, trapsX, trapsY, 
                             S_usage, N, 
                             sigma_0, sd_sigma0,
                             sd_sigma_prop,
                             a_p, b_p){
  
  
  # update p0
  
  c_bar <- sum(CH)
  
  k_i <- matrix(NA, N, K)
  
  for (i in 1:N) {
    for (j in 1:K) {
      k_i[i, j] <- S_usage[j] - CH[i,j]
    }
  }
  
  c_i <- matrix(NA, N, K)
  
  for (i in 1:N) {
    for (j in 1:K) {
      c_i[i, j] <- exp(- (1 / (2 * sigma^2)) * ( (s[i,1] - trapsX[j])^2 + (s[i,2] - trapsY[j])^2))
    }
  }
  
  c_i <- as.vector(c_i)
  
  a_fp <- c_bar
  b_fp <-  sum(- c_i * k_i)
  c_fp <- sum(k_i * c_i^2)
  
  p_star <- (- b_fp - sqrt(b_fp^2 - 4 * a_fp * c_fp)) / (2 * c_fp)
  h_star <- h_star <- 1 / (a_fp / p_star^2)
  
  p_new <- rnorm(1, p_star, sqrt(h_star))
  
  if(p_new > 0 & p_new < 1){
    
    loglikelihood_p0_star <- loglikelihood_p_cpp2(p_new, sigma, CH, s,
                                                  trapsX, trapsY, S_usage, N)
    
    loglikelihood_p0 <- loglikelihood_p_cpp2(p0, sigma, CH, s,
                                             trapsX, trapsY, S_usage, N)
    
    logprior_p0_star <- dbeta(p_new, a_p, b_p, log = T)
    
    logprior_p0 <- dbeta(p0, a_p, b_p, log = T)
    
    # logproposal_p0_current <- logproposal_p0(p0_2, CH2, mean(c_i))
    # 
    # logproposal_p0_star <- logproposal_p0(p0_2_star, CH2, mean(c_i))
    
    logposterior_p0_star <- logprior_p0_star + loglikelihood_p0_star
    
    logposterior_p0 <- logprior_p0 + loglikelihood_p0
    
    MH_ratio <- exp(logposterior_p0_star - logposterior_p0)
    
    if(runif(1) < MH_ratio){
      p0 <- p_new
    }
    
  }
  
  # update sigma
  
  sigma_star <- rnorm(1, sigma, sd = sd_sigma_prop)
  
  if(sigma_star > 0){
    
    loglikelihood_sigma_star <- loglikelihood_p_cpp2(p0, sigma_star, CH, s,
                                                     trapsX, trapsY, S_usage, N)
    
    loglikelihood_sigma <- loglikelihood_p_cpp2(p0, sigma, CH, s,
                                                trapsX, trapsY, S_usage, N)
    
    logprior_sigma_star <- dnorm(sigma_star, sigma_0, sd_sigma0, log = T)
    
    logprior_sigma <- dnorm(sigma, sigma_0, sd_sigma0, log = T)
    
    logposterior_sigma_star <- loglikelihood_sigma_star + logprior_sigma_star
    
    logposterior_sigma <- loglikelihood_sigma + logprior_sigma
    
    if(runif(1) < exp(loglikelihood_sigma_star - loglikelihood_sigma)){
      sigma <- sigma_star
    }
    
  }
  
  list("p0" = p0, "sigma" = sigma)
}

update_s <- function(s1, s2, theta1, theta2, theta12, beta1, beta2, sigma_prop, N1, N2){
  
  log_hastings_ratio_den <- log_f_bivsoftcore_cpp(s1, s2, N1, N2,
                                                  beta1, beta2,
                                                  theta1, theta2, theta12)
  
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
                                                                 theta1, theta2, theta12)
      
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
                                                                 theta1, theta2, theta12)
      
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
                                                                theta1, theta2, theta12)
    
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
                                                                theta1, theta2, theta12)
    
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
                                                                theta1, theta2, theta12)
    
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
                                                               theta1, theta2, theta12)
    
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

update_s_poissonr <- function(s1, s2, theta1, theta2, theta12, beta1, beta2, sigma_prop, N1, N2,
                              lambda_points1, lambda_points2){

  q <- .333333333
  p <- .5
  
  # first process
  
  if(runif(1) < q){
    
  } else {
    
    if(runif(1) < p){
    
      pointsToPropose <- 1 + rpois(1, lambda_points1)
      x_new <- t(sapply(1:pointsToPropose, function(i){
        proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
      }))
      
      log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_addition(x_new, pointsToPropose,
                                                                     s1, s2, N1, N2,
                                                                     beta1, theta1, theta12)
      
      loglikelihood_xi_star <- 0
      for(l in 1:pointsToPropose){
        loglikelihood_xi_star <- loglikelihood_xi_star + 
          log(probAcceptingIndividual(s1, x_new[i,], p0_1, alpha_1, D1, N1 + l,
                                      K, S_usage, trapsX, trapsY))
      }
        
    } else {
      
      pointsToRemove = 1 + rpois(1, lambda_points1)
      
      if(N1 - pointsToRemove >= D1){
        
        idxPointsToRemove <- D1 + sample.int(N1 - D1, pointsToRemove)
        
        log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_removal(idxPointsToRemove, pointsToRemove,
                                                                     s1, s2, N1, N2,
                                                                     beta1, theta1, theta12)
        
        loglikelihood_xi_star <- 0
        for(l in seq_along(idxPointsToRemove)){
          loglikelihood_xi_star <- loglikelihood_xi_star + 
            log(probAcceptingIndividual(s1, s1[idxPointsToRemove[l],], p0_1, alpha_1, D1, N1 + l,
                                        K, S_usage, trapsX, trapsY))
        }
        
        logposterior_xistar <- log_hastings_ratio_num
        logposterior_xi <- loglikelihood_xi_star
        
        if(runif(1)  < exp(logposterior_xistar - logposterior_xi)){
          
        }
        
      }
      
    }
    
    
  }
  
  # second process
  
  if(runif(1) < q){
    
  } else {
    
    if(runif(1) < p){
      
      pointsToPropose <- 1 + rpois(1, lambda_points2)
      x_new <- t(sapply(1:pointsToPropose, function(i){
        proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
      }))
      
      log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_addition(x_new, pointsToPropose,
                                                                     s2, s1, N2, N1,
                                                                     beta2, theta2, theta12)
      
      loglikelihood_xi_star <- 0
      for(l in 1:pointsToPropose){
        loglikelihood_xi_star <- loglikelihood_xi_star + 
          log(probAcceptingIndividual(s1, x_new[i,], p0_1, alpha_1, D1, N1 + l,
                                      K, S_usage, trapsX, trapsY))
      }
      
    } else {
      
      pointsToRemove = 1 + rpois(1, lambda_points1)
      
      if(N1 - pointsToRemove >= D1){
        
        idxPointsToRemove <- D1 + sample.int(N1 - D1, pointsToRemove)
        
        log_hastings_ratio_num <- log_f_bivsoftcore_cpp_quick_removal(idxPointsToRemove, pointsToRemove,
                                                                      s1, s2, N1, N2,
                                                                      beta1, theta1, theta12)
        
        loglikelihood_xi_star <- 0
        for(l in seq_along(idxPointsToRemove)){
          loglikelihood_xi_star <- loglikelihood_xi_star + 
            log(probAcceptingIndividual(s1, s1[idxPointsToRemove[l],], p0_1, alpha_1, D1, N1 + l,
                                        K, S_usage, trapsX, trapsY))
        }
        
        logposterior_xistar <- log_hastings_ratio_num
        logposterior_xi <- loglikelihood_xi_star
        
        if(runif(1)  < exp(logposterior_xistar - logposterior_xi)){
          
        }
        
      }
      
    }
    
    
  }
  
  list("s1" = s1, "s2" = s2, 
       "N1" = N1, "N2" = N2)
}

log_prior <- function(theta1, theta2, theta12, beta1, beta2,
                      # a_theta, b_theta){
                      mu_logtheta, sd_logtheta, a_beta, b_beta){
  
  dnorm(log(theta1), mu_logtheta, sd_logtheta, log = T) + 
    dnorm(log(theta2), mu_logtheta, sd_logtheta, log = T) + 
    dnorm(log(theta12), mu_logtheta, sd_logtheta, log = T) + 
  # dgamma(theta1, a_theta, b_theta, log = T) + dgamma(theta2, a_theta, b_theta, log = T) + 
  #   dgamma(theta12, a_theta, b_theta, log = T) + 
    dgamma(beta1, a_beta, b_beta, log = T) + dgamma(beta2, a_beta, b_beta, log = T)
  
}

log_gamma_prior <- function(gamma, mu, sigmasq){
  
   - (1 / (2 * sigmasq)) * (sqrt(log(2) /  gamma) - mu)^2  +  (-3/2) * log(gamma)
  
}

hastings_ratio <- function(data1, data2, N1, N2, 
                           x_all, 
                           N_all,
                           beta1, beta2, 
                           beta1_star, beta2_star, 
                           theta1, theta2, theta12,
                           theta1_star, theta2_star, theta12_star,
                           # mu_gamma, sigmasq_gamma,
                           # a_theta, b_theta,
                           mu_logtheta, sd_logtheta,
                           a_beta, b_beta){
  
  x_all1 <- x_all[,1,,]
  x_all2 <- x_all[,2,,]
  
  ratio_constants <- 1 /  mean(importanceSamplingEstimateParallel(x_all1, x_all2, 
                                                             c(beta1, beta2, 
                                                               theta1, theta2, theta12),
                                                             c(beta1_star, beta2_star, 
                                                               theta1_star, theta2_star, theta12_star),
                                                             N_all))
  # ratio_constants <- 1 / importance_sampling_estimate_foreach(x_all, N_all,
  #                                                         beta1, beta2,
  #                                                         beta1_star, beta2_star, 
  #                                                         theta1, theta2, theta12,
  #                                                         theta1_star, theta2_star, theta12_star)
  
  den <- log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1, beta2, theta1, theta2, theta12) + 
    log_prior(theta1, theta2, theta12, beta1, beta2, mu_logtheta, sd_logtheta, a_beta, b_beta)
    # log_gamma_prior(theta1, mu_gamma, sigmasq_gamma) + 
    # log_gamma_prior(theta12, mu_gamma, sigmasq_gamma) +
    # dgamma(beta1, a_beta, b_beta, log = T) + dgamma(beta2, a_beta, b_beta, log = T)#
  
  
  num <- log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1_star, beta2_star, theta1_star, theta2_star, theta12_star) + 
    log_prior(theta1_star, theta2_star, theta12_star, beta1_star, beta2_star, mu_logtheta, sd_logtheta, a_beta, b_beta)
    # log_gamma_prior(theta1_star, mu_gamma, sigmasq_gamma) +
    # log_gamma_prior(theta12_star, mu_gamma, sigmasq_gamma) +
    # dgamma(beta1_star, a_beta, b_beta, log = T) + dgamma(beta2_star, a_beta, b_beta, log = T)
  
  exp(num - den) * ratio_constants
}

hastings_ratio_2 <- function(data1, data2, N1, N2, 
                             x_all1, x_all2, 
                             N_all,
                             beta1, beta2, 
                           beta1_star, beta2_star, 
                           theta1, theta2, theta12,
                           theta1_star, theta2_star, theta12_star,
                           # mu_gamma, sigmasq_gamma,
                           # a_theta, b_theta,
                           mu_logtheta, sd_logtheta,
                           a_beta, b_beta){
  
  x_all1 <- x_all[,1,,]
  x_all2 <- x_all[,2,,]
  
  ratio_constants <- 1 /  mean(importanceSamplingEstimateParallel(x_all1, x_all2, 
                                                             c(beta1, beta2, 
                                                               theta1, theta2, theta12),
                                                             c(beta1_star, beta2_star, 
                                                               theta1_star, theta2_star, theta12_star),
                                                             N_all))
  # ratio_constants <- 1 / importance_sampling_estimate_foreach(x_all, N_all,
  #                                                         beta1, beta2,
  #                                                         beta1_star, beta2_star, 
  #                                                         theta1, theta2, theta12,
  #                                                         theta1_star, theta2_star, theta12_star)
  
  den <- log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1, beta2, theta1, theta2, theta12) + 
    log_prior(theta1, theta2, theta12, beta1, beta2, mu_logtheta, sd_logtheta, a_beta, b_beta)
    # log_gamma_prior(theta1, mu_gamma, sigmasq_gamma) + 
    # log_gamma_prior(theta12, mu_gamma, sigmasq_gamma) +
    # dgamma(beta1, a_beta, b_beta, log = T) + dgamma(beta2, a_beta, b_beta, log = T)#
  
  
  num <- log_f_bivsoftcore_cpp(data1, data2, N1, N2, beta1_star, beta2_star, theta1_star, theta2_star, theta12_star) + 
    log_prior(theta1_star, theta2_star, theta12_star, beta1_star, beta2_star, mu_logtheta, sd_logtheta, a_beta, b_beta)
    # log_gamma_prior(theta1_star, mu_gamma, sigmasq_gamma) +
    # log_gamma_prior(theta12_star, mu_gamma, sigmasq_gamma) +
    # dgamma(beta1_star, a_beta, b_beta, log = T) + dgamma(beta2_star, a_beta, b_beta, log = T)
  
  exp(num - den) * ratio_constants
}

importance_sampling_estimate <- function(x_all, N_all,
                                         beta1, beta2, 
                                         beta1_star, beta2_star, 
                                         theta1, theta2, theta12,
                                         theta1_star, theta2_star, theta12_star){
  
  constant <- 0
  M <- dim(x_all)[1]
  
  for (m in 1:M) {
    
    N1 <- N_all[m,1]
    N2 <- N_all[m,2]
    x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
    x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
    
    num <- log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                 beta1_star, beta2_star, 
                                 theta1_star, theta2_star, theta12_star)
    
    den <- log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                 beta1, beta2, 
                                 theta1, theta2, theta12)
    
    constant <- constant + exp(num - den)
    print(exp(num - den))
  }
  
  constant / M
}

importance_sampling_estimate_foreach <- function(x_all, N_all,
                                                 beta1, beta2,
                                                 beta1_star, beta2_star,
                                                 theta1, theta2, theta12,
                                                 theta1_star, theta2_star, theta12_star){
  
  r <- foreach(m = 1:50, .combine=c, 
               .packages = "SpatialFunctionsCR") %dopar% {
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 x1 <- matrix(x_all[m,1,seq_len(N1),], nrow = N1)
                 x2 <- matrix(x_all[m,2,seq_len(N2),], nrow = N2)
                 
                 num <-  SpatialFunctionsCR::log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                              beta1_star, beta2, 
                                              theta1_star, theta2, theta12_star)
                 
                 den <-  SpatialFunctionsCR::log_f_bivsoftcore_cpp(x1, x2, N1, N2, 
                                              beta1, beta2, 
                                              theta1, theta2, theta12)
                 
                 exp(num - den)         
               }
  
  mean(r)
}

update_x_all_foreach <- function(x_all, N_all, niter,
                                 theta1, theta2, theta12, 
                                 beta1, beta2, Sigma_prop, Sigma_newpoint, 
                                 allTraps, R_traps,
                                 a1, b1, a2, b2){
  
  M <- dim(x_all)[1]
  
  r <- foreach(m = 1:M, .combine='c', .multicombine=F,
               .packages = "SpatialFunctionsCR") %dopar% {
                 
                 N1 <- N_all[m,1]
                 N2 <- N_all[m,2]
                 x1 <- x_all[m,1,,]
                 x2 <- x_all[m,2,,]
                 
                 list_sims <- simulate_bivsoftcore_from_startingpoint(x1, x2, N1, N2, 
                                                                      theta1, theta2, theta12,
                                                                      beta1, beta2,
                                                                      niter, Sigma_prop = diag(0.01, nrow = 2), 
                                                                      Sigma_newpoint, 
                                                                      allTraps, R_traps,
                                                                      a1, b1, a2, b2)
                 
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
