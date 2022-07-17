
meanBeta1 <- 250
varBeta1 <- 10000
b_beta1 <- meanBeta1 / varBeta1
a_beta1 <- meanBeta1 * b_beta1  


meanBeta2 <- 250
varBeta2 <- 10000
b_beta2 <- meanBeta2 / varBeta2
a_beta2 <- meanBeta2 * b_beta2  

mu_logtheta <- -30
sd_logtheta <- 20

nsims <- 2000
N_output <- matrix(NA, nsims, 2)

for (sim in 1:nsims) {
  
  beta1 <- rgamma(1, shape = a_beta1, rate = b_beta1)
  beta2 <- rgamma(1, shape = a_beta2, rate = b_beta2)
  
  theta1 <- exp(dnorm(1, mu_logtheta, sd_logtheta))
  theta2 <- exp(dnorm(1, mu_logtheta, sd_logtheta))
  theta12 <- exp(dnorm(1, mu_logtheta, sd_logtheta))
  
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
  
  N_output[iter,] <- c(N1, N2)
}



