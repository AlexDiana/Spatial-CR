simulate_CH_from_s <- function(p0_1, sigma_1, 
                               p0_2, sigma_2,
                               s1, N1, s2, N2){
  
  K <- nrow(traps)
  
  list_CH <- simulate_CH_from_p_cpp(p0_1, sigma_1, p0_2, sigma_2,
                                    s1, N1, s2, N2,
                                    S_usage, trapsX, trapsY)
  CH_1 <- list_CH$CH_1
  CH_2 <- list_CH$CH_2
  
  indexesCaughtIndividuals1 <- which(apply(CH_1,1,sum) != 0)
  indexesUncaughtIndividuals1 <- which(apply(CH_1,1,sum) == 0)
  D1 <- length(indexesCaughtIndividuals1)
  CH_1 <- CH_1[c(indexesCaughtIndividuals1, indexesUncaughtIndividuals1),]
  
  indexesCaughtIndividuals2 <- which(apply(CH_2,1,sum) != 0)
  indexesUncaughtIndividuals2 <- which(apply(CH_2,1,sum) == 0)
  D2 <- length(indexesCaughtIndividuals2)
  CH_2 <- CH_2[c(indexesCaughtIndividuals2, indexesUncaughtIndividuals2),]
  
  s1_0 <- matrix(NA, nrow = Nmax, ncol = 2)
  s2_0 <- matrix(NA, nrow = Nmax, ncol = 2)
  
  s1_0[seq_len(D1),] <- s1[indexesCaughtIndividuals1,]
  if(N1 > D1){
    s1_0[(D1 + 1):N1,] <- s1[indexesUncaughtIndividuals1,]
  }
  s2_0[seq_len(D2),] <- s2[indexesCaughtIndividuals2,]
  if(N2 > D2){
    s2_0[(D2 + 1):N2,] <- s2[indexesUncaughtIndividuals2,]  
  }
  
  list("CH_1" = CH_1,
       "CH_2" = CH_2,
       "s1" = s1_0,
       "s2" = s2_0,
       "D1" = D1,
       "D2" = D2)
}

# FIXED VALUES -------------------------------

# other parameters

{
  p0_1 <- .005#.001
  p0_2 <- .005#.001
  sigma_1 <- 1.5 / trapsX_sd
  sigma_2 <- 1.5 / trapsX_sd
}

# prior

{
  theta1 <- exp(-100)
  theta2 <- exp(-100)
  theta12 <- exp(-100)
  beta1 <- 250
  beta2 <- 250
}

# proposals
{
  sigma_prop <- .01
  Sigma_newpoint <- cov(cbind(trapsX, trapsY))
  
  normConstUnif <- mean(sapply(1:M, function(j){
    x <- c(runif(1, a1, b1), 
           runif(1, a2, b2))
    checkPointIsInRegionTraps(x, allTraps, R_traps)
  }))
  
  S_area <- (b1 - a1) * (b2 - a2) * normConstUnif
  
  lambda_newPoints1 <- 5
  lambda_newPoints2 <- 5
  
  # mixture prop
  {
    best_k <- 2
    
    mixtureMeans1 <- matrix(c(0,0.5,
                              0,0.5), nrow = best_k, ncol = 2)
    mixtureSd1 <- matrix(c(.5,.5,
                           .5,.5), nrow = best_k, ncol = 2)
    normConstNormals1 <- rep(NA, best_k)
    
    for (i in 1:best_k) {
      
      normConstNormals1[i] <- mean(sapply(1:M, function(j){
        x <- c(rnorm(1, mixtureMeans1[i,1], mixtureSd1[i,1]), 
               rnorm(1, mixtureMeans1[i,2], mixtureSd1[i,2]))
        checkPointIsInRegionTraps(x, allTraps, R_traps)
      }))
    }
    
    normConst1 <- c(normConstUnif, normConstNormals1)  
    
    mixtureMeans2 <- mixtureMeans1
    mixtureSd2 <- mixtureSd1
    normConst2 <- normConst1
  }
  
  probMixture <- c(0.5, 0.5)
}

# STARTING VALUES ---------------------------------------------------------

list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                      Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                      Sigma_newpoint = diag(1, nrow = 2), allTraps, .2,
                                      a1, b1, a2, b2)
s1 <- list_s1s2$data1
s2 <- list_s1s2$data2
N1 <- list_s1s2$N1
N2 <- list_s1s2$N2

print(N1)
print(N2)

# MCMC --------------------------------------------------------------------

niter <- 10000

N_iter <- matrix(NA, nrow = niter, ncol = 2)

for (iter in 1:niter) {
  
  print(iter)
  print(paste0("N1 = ",N1," - N2 = ", N2))
  
  # SIMULATE CH --------------------------------------
  
  list_CH <- simulate_CH_from_s(p0_1, sigma_1, 
                                p0_2, sigma_2,
                                s1, N1, s2, N2)
  CH1 <- list_CH$CH_1
  CH2 <- list_CH$CH_2
  s1 <- list_CH$s1
  s2 <- list_CH$s2
  D1 <- list_CH$D1
  D2 <- list_CH$D2
  
  # UPDATE HOME CENTERS -----------------------------------------------
  
  # list_s <- update_s_cpp_new(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
  #                            CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
  #                            Sigma_newpoint, lambda_movedPoints1 = 20, lambda_movedPoints2 = 20,
  #                            lambda_newPoints1 = 0, lambda_newPoints2 = 0, probMixture[1],
  #                            probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
  #                            normConst1, normConst2, allTraps, R_traps, S_area, S_usage, trapsX,
  #                            trapsY, a1, b1, a2, b2)
  # s1 <- list_s$s1
  # s2 <- list_s$s2
  # N1 <- list_s$N1
  # N2 <- list_s$N2
  
  # UPDATE N ---------------
  
  # list_s <- update_N_cpp_new(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
  #                            CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, 
  #                            Sigma_prop = diag(.1, nrow = 2),
  #                            Sigma_newpoint, lambda_movedPoints1 = 20, lambda_movedPoints2 = 20,
  #                            lambda_newPoints1, lambda_newPoints2, probMixture[1],
  #                            probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
  #                            normConst1, normConst2, allTraps, R_traps, S_area, S_usage, trapsX,
  #                            trapsY, a1, b1, a2, b2)
  
  list_s <-  update_s_cpp_withacceptance(s1, s2,  theta1, theta2, theta12, beta1, beta2, CH1,
                                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
                                         allTraps, R_traps,
                                         a1, b1, a2, b2)
  
  s1 <- list_s$s1
  s2 <- list_s$s2
  N1 <- list_s$N1
  N2 <- list_s$N2
  
  N_iter[iter,] <- c(N1,N2)  
  
}

apply(N_iter, 2, mean)
apply(N_iter, 2, var)

ggplot(data = NULL, aes(x = 1:niter, y = N_iter[1:niter,1])) + geom_point() + 
  geom_hline(aes(yintercept = beta1))
ggplot(data = NULL, aes(x = 1:niter, y = N_iter[1:niter,2])) + geom_point() + 
  geom_hline(aes(yintercept = beta1))
