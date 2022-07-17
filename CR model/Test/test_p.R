simulate_CH_from_p <- function(p0_1, alpha_1, 
                               p0_2, alpha_2,
                               s1, N1, s2, N2){
  
  CH_1 <- matrix(NA, nrow = N1, ncol = K)
  
  for (i in 1:N1) {
    for (k in 1:K) {
      p <- p0_1 * exp(- alpha_1 * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
      CH_1[i,k] <- rbinom(1, S, p)
    }
  }
  
  CH_2 <- matrix(NA, nrow = N2, ncol = K)
  for (i in 1:N2) {
    for (k in 1:K) {
      p <- p0_2 * exp(- alpha_2 * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
      CH_2[i,k] <- rbinom(1, S, p)
    }
  }
  
  list("CH_1" = CH_1,
       "CH_2" = CH_2)
}

# FIXED VALUES -------------------------------

# other parameters

{
  theta1 <- exp(-100)
  theta2 <- exp(-100)
  theta12 <- exp(-100)
  beta1 <- 100
  beta2 <- 200
  
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 1000, beta1, 
                                        Sigma_newpoint = diag(1, nrow = 2),
                                        allTraps, R_traps,
                                        a1, b1, a2, b2)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  
  s1 <- s1[1:N1,]
  s2 <- s2[1:N2,]
  
}

# prior

{
  a_p <- 1
  b_p <- 100
  
  sigma_0_1 <- 2.5 / traps_meansd
  sd_sigma_1 <- .5 / traps_meansd
  
  sigma_0_2 <- 10 / traps_meansd
  sd_sigma_2 <- .5 / traps_meansd
  
  # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = alpha_0_1, sd = sd_alpha_1))
  # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = alpha_0_2, sd = sd_alpha_2))
  
  sd_sigma_prop <- .001 
}

# STARTING VALUES ---------------------------------------------------------

p0_1 <- a_p / (a_p + b_p)
sigma_1 <- sigma_0_1
p0_2 <- .02#.001
sigma_2 <- sigma_0_2

# MCMC --------------------------------------------------------------------

niter <- 10000

sigma_iter <- array(NA, dim = c(niter, 2))
p0_iter <- array(NA, dim = c(niter, 2))

for (iter in 1:niter) {
  
  if(iter %% 50 == 0 ){
    print(iter)  
  }
  
  list_CH <- simulate_CH_from_p_cpp(p0_1, sigma_1, 
                                p0_2, sigma_2,
                                s1, N1, s2, N2,
                                S_usage, trapsX, trapsY)
  CH1 <- list_CH$CH_1
  CH2 <- list_CH$CH_2
  # list_CH <- simulate_CH_from_p(p0_1, alpha_1, 
  #                               p0_2, alpha_2,
  #                               s1, N1, s2, N2)
  # CH1 <- list_CH$CH_1
  # CH2 <- list_CH$CH_2
  
  # indexesCaughtIndividuals1 <- which(apply(CH1,1,sum) != 0)
  # indexesUncaughtIndividuals1 <- which(apply(CH1,1,sum) == 0)
  # D1 <- length(indexesCaughtIndividuals1)
  # CH1 <- CH1[c(indexesCaughtIndividuals1, indexesUncaughtIndividuals1),]
  # 
  # s1 <- s1[c(indexesCaughtIndividuals1, indexesUncaughtIndividuals1),]
  # 
  # indexesCaughtIndividuals2 <- which(apply(CH2,1,sum) != 0)
  # indexesUncaughtIndividuals2 <- which(apply(CH2,1,sum) == 0)
  # D2 <- length(indexesCaughtIndividuals2)
  # CH2 <- CH2[c(indexesCaughtIndividuals2, indexesUncaughtIndividuals2),]
  # 
  # s2 <- s2[c(indexesCaughtIndividuals2, indexesUncaughtIndividuals2),]
  
  # list_p <- update_p(p0_1, alpha_1, p0_2, alpha_2, s1, s2, CH1, CH2, trapsX, trapsY,
  #                    S, N1, N2, 
  #                    alpha_0_1, sd_alpha_1,
  #                    alpha_0_2, sd_alpha_2,
  #                    sigma_p0_prop, sigma_alpha_prop,
  #                    a_p, b_p)
  # p0_1 <- list_p$p0_1
  # p0_2 <- list_p$p0_2
  # alpha_1 <- list_p$alpha_1
  # alpha_2 <- list_p$alpha_2
  
  list_p <- update_psigma_cpp(p0_1, sigma_1, CH1, s1, N1, 
                              K, S_usage, trapsX, trapsY,
                              sigma_0_1, sd_sigma0 = exp(-120), sd_sigma_prop = 100,
                              a_p, b_p)
  p0_1 <- list_p[1]
  sigma_1 <- list_p[2]
  
  # list_p2 <- update_psigma_cpp(p0_2, sigma_2, CH2, s2, N2, 
  #                              K, S_usage, trapsX, trapsY,
  #                              sigma_0_2, sd_sigma_2, sd_sigma2_prop,
  #                              a_p, b_p)
  # p0_2 <- list_p2[1]
  # sigma_2 <- list_p2[2]
  
  p0_iter[iter,] <- c(p0_1, p0_2)
  sigma_iter[iter,] <- c(sigma_1, sigma_2)
  
}

qplot(1:niter, p0_iter[,1]) + geom_hline(aes(yintercept = (a_p / (a_p + b_p))))
qplot(1:niter, sigma_iter[,1])
