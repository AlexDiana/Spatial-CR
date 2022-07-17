library(MASS); library(reshape2); library(ggplot2)
library(foreach); library(doParallel)
library(Rcpp); library(RcppArmadillo)
library(here)

sourceCpp("CR model/Simulations/General/code.cpp")
source('CR model/functions.r', echo=TRUE)

# SETTINGS ---

a1 <- 0
b1 <- 1
a2 <- 0
b2 <- 1

S_area1 <- 1
S_area2 <- 1

numTraps_grid <- c(15, 10)
N_grid <- c(100, 200)
interRatioStrength_grid <- c(exp(-6), exp(-7.5))
p_grid <- c(.3, .03)

nruns <- 10
usingTrueData <- F

logparams_ratio <- array(NA, dim = c(nruns, 
                                     length(numTraps_grid), 
                                     length(N_grid),
                                     length(interRatioStrength_grid),
                                     length(p_grid),
                                     2, # theta 1/12 or theta 2/theta12
                                     2))

updateTheta <- T

run <- 1
idx_N_ind <- 1
idx_interRatioStrength <- 1
idx_p <- 1
idx_n_traps <- 1

for (run in 1:nruns) {
  
  for (idx_N_ind in seq_len(N_grid)) {
    
    N_ind <- N_grid[idx_N_ind]
    
    for (idx_interRatioStrength in seq_len(interRatioStrength_grid)) {
      
      interRatioStrength <- interRatioStrength_grid[idx_interRatioStrength]
      
      # SIMULATE INDIVIDUALS -------
      
      N1 <- N_ind
      N2 <- N_ind
      
      theta1 <- exp(-100)
      theta2 <- exp(-100)
      theta12 <- interRatioStrength
      
      theta1_true <- theta1
      theta2_true <- theta2
      theta12_true <- theta12
      
      list_s1s2 <- simulate_cond_bivsoftcore_cpp(theta1, theta2, theta12, N1, N2, 2000, 
                                                 Sigma_prop = diag(0.01, nrow = 2), 2000, N1, 
                                                 Sigma_newpoint = diag(1, nrow = 2), 
                                                 a1, b1, a2, b2)
      s1 <- list_s1s2$data1
      s2 <- list_s1s2$data2
      
      ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") +
        geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
      
      for (idx_n_traps in seq_along(numTraps_grid)) {
        
        numTraps <- numTraps_grid[idx_n_traps]
                
        # CREATE NEW TRAPS -------------------
        
        # trapsX_grid <- seq(a1, b1, length.out = numTraps + 2)[-(numTraps+2)][-1]
        # trapsY_grid <- seq(a2, b2, length.out = numTraps + 2)[-(numTraps+2)][-1]
        trapsX_grid <- seq(a1, b1, length.out = numTraps + 1)
        trapsX_grid <- (trapsX_grid[-1] + trapsX_grid[-length(trapsX_grid)]) / 2
        trapsY_grid <- seq(a2, b2, length.out = numTraps + 1)
        trapsY_grid <- (trapsY_grid[-1] + trapsY_grid[-length(trapsY_grid)]) / 2
        traps <- expand.grid(trapsX_grid, trapsY_grid)
        
        ggplot(data = NULL, aes(x = traps[,1], y = traps[,2])) + geom_point() + 
          xlim(c(0,1)) + ylim(c(0,1))
        
        trapsX <- traps[,1]
        trapsY <- traps[,2]
        
        trapsX_mean <- mean(trapsX)
        trapsY_mean <- mean(trapsY)
        trapsX_sd <- sd(trapsX)
        trapsY_sd <- sd(trapsY)
        
        K <- length(trapsX)
        
        S_usage <- rep(10, K)
        
        traps_meansd <- (sd(trapsX) + sd(trapsY)) / 2
        
        trapsX <- trapsX
        trapsY <- trapsY
        
        allTraps <- cbind(trapsX, trapsY)
        R_traps <- .75
        
        for (idx_p in p_grid) {
          
          p <- p_grid[idx_p]
          
          p0_1 <- p
          p0_2 <- p
          
          p0_1_true <- p0_1
          p0_2_true <- p0_2
          
          # SIMULATE CAPTURE HISTORIES -----
          
          sigma_1 <- .1
          sigma_2 <- .1
          
          sigma_1_true <- sigma_1
          sigma_2_true <- sigma_2
          
          K <- nrow(traps)
          
          CH_1 <- matrix(NA, nrow = N1, ncol = K)
          for (i in 1:N1) {
            for (k in 1:K) {
              p <- p0_1 * exp(-(1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
              CH_1[i,k] <- rbinom(1, S_usage[k], p)
            }
          }
          indexesCaughtIndividuals1 <- which(apply(CH_1,1,sum) != 0)
          indexesUncaughtIndividuals1 <- which(apply(CH_1,1,sum) == 0)
          D1 <- length(indexesCaughtIndividuals1)
          CH_1 <- CH_1[indexesCaughtIndividuals1,]
          
          CH_2 <- matrix(NA, nrow = N2, ncol = K)
          for (i in 1:N2) {
            for (k in 1:K) {
              p <- p0_2 * exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
              CH_2[i,k] <- rbinom(1, S_usage[k], p)
            }
          }
          indexesCaughtIndividuals2 <- which(apply(CH_2,1,sum) != 0)
          indexesUncaughtIndividuals2 <- which(apply(CH_2,1,sum) == 0)
          D2 <- length(indexesCaughtIndividuals2)
          CH_2 <- CH_2[indexesCaughtIndividuals2,]
          
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
          
          print(D1)
          print(D2)
          
          # captured <- factor(c(rep(1, D1),rep(0, N1 - D1)))
          # ggplot() +  geom_point(data = NULL, aes(x = s1_0[,1], y = s1_0[,2], 
          #                                         color = captured)) + 
          #   geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "black", size = 2)
          # 
          # captured <- rep(0,N1)
          # captured[indexesCaughtIndividuals1] <- 1
          # # captured <- factor(c(rep(1, D1),rep(0, N1 - D1)))
          # ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2], 
          #                                         color = factor(captured))) + 
          #   geom_point(data = NULL, aes(x = trapsX, y = trapsY), color = "red", size = 2)
          
          # CLEAN DATA --------------------------------------------------------------
          
          D1 <- nrow(CH_1)
          D2 <- nrow(CH_2)
          
          CH1 <- CH_1
          CH2 <- CH_2
          
          
          
          # PRIOR -------------------------------------------------------------------
          
          Nmax <- 1200
          
          # proposal interaction parameters
          {
            epsilon_beta <- .5^2
            epsilon_theta1 <- .000005^2
            epsilon_theta2 <- .000005^2
            epsilon_theta12 <- .000005^2 
          }
          
          # intensity parameter
          {
            meanBeta <- N1
            varBeta <- 10000
            b_beta <- meanBeta / varBeta
            a_beta <- meanBeta * b_beta  
            
            mu_logtheta <- -30
            sd_logtheta <- 20
            
            # ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,400))
          }
          
          # proposal value
          
          sigma_prop <- .01
          Sigma_newpoint <- cov(cbind(trapsX, trapsY))
          lambda_newPoints1 <- 1
          lambda_newPoints2 <- 1
          
          # capture probability params
          {
            a_p <- 1
            b_p <- 1
            
            # leopards (.5 to 2.5)
            sigma_0_1 <- sigma_1
            sd_sigma_1 <- 2 / traps_meansd
            
            # tigers ( 1.5 to 3.5)
            sigma_0_2 <- sigma_2
            sd_sigma_2 <- 2 / traps_meansd
            
            sigma_p0_1_prop <- .0002
            sigma_sigma1_prop <- .005
            
            sigma_p0_2_prop <- .0002
            sigma_sigma2_prop <- .01
          }
          
          # MCMC --------------------------------------------------------------------
          
          pointsToAccept <- 1000#nburn <- 0
          niter <- 1000
          nchain <- 1
          nthin <- 1
          nburn <- 500
          iterAfterAdapting <- 200
          updateTheta <- T
          iterToDiscard <- 1
          beta_proposal <- .1
          
          # variables for output
          {
            sigma_output <- array(NA, dim = c(nchain, niter, 2))
            s_output <- array(NA, dim = c(nchain, niter, 2, Nmax, 2))
            N_output <- array(NA, dim = c(nchain, niter, 2))
            p0_output <- array(NA, dim = c(nchain, niter, 2))
            params_output <- array(NA, dim = c(nchain, niter, 5))
            
            # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
            # params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
            
            # acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
          }
          
          for(chain in 1:nchain) {
            
            # home centers
            {
              N1 <- N1_0
              N2 <- N2_0
              
              CH1 <- rbind(CH1[1:D1,],
                           matrix(0, nrow = N1 - D1, ncol = K))
              
              CH2 <- rbind(CH2[1:D2,],
                           matrix(0, nrow = N2 - D2, ncol = K))
              
              s1 <- matrix(NA, nrow = Nmax, ncol = 2)
              s2 <- matrix(NA, nrow = Nmax, ncol = 2)
              
              s1[1:N1,] <- s1_0
              s2[1:N2,] <- s2_0
              
            }
            
            # capture probability parameter
            {
              p0_1 <- p0_1_true
              sigma_1 <- sigma_1_true
              p0_2 <- p0_2_true
              sigma_2 <- sigma_2_true
              
            }
            
            # point process parameters
            {
              theta1 <- theta1_true
              theta2 <- theta2_true
              theta12 <- theta12_true
              beta1 <- N1
              beta2 <- N2
              M <- 50
              
              x_all <- array(NA, dim = c(M, 2, Nmax, 2))
              N_all <- matrix(NA, nrow = M, ncol = 2)
              
              for (m in 1:M) {
                
                print(m)
                
                list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                                      niter = 2000,
                                                      Sigma_prop = diag(.01, nrow = 2), 
                                                      Nmax = Nmax, lambda = beta1,
                                                      Sigma_newpoint = diag(1, nrow = 2),
                                                      a1, b1, a2, b2)
                N1_sim <- list_sims$N1
                N2_sim <- list_sims$N2
                x_all[m,1,1:N1_sim,] <- list_sims$data1[1:N1_sim,]
                x_all[m,2,1:N2_sim,] <- list_sims$data2[1:N2_sim,]
                N_all[m,] <- c(N1_sim, N2_sim)
                
              }
              
              params_values <- matrix(NA, nrow = pointsToAccept * 10 + niter, ncol = 5)
              
              
            }
            
            # output variables
            {
              sigma_iter <- array(NA, dim = c(niter, 2))
              s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
              N_iter <- matrix(NA, nrow = niter, ncol = 2)
              p0_iter <- array(NA, dim = c(niter, 2))
              params_iter <- matrix(NA, nrow = niter, ncol = 5)
              
              # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
              params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
              
              acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
              theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
            }
            
            for(iter in 1:(nburn + niter)){
              
              if(iter < nburn){
                print(paste0("Chain = ",chain," / Burn-in Iteration = ", iter)) 
              } else {
                print(paste0("Chain = ",chain," / Iteration = ", iter - nburn))
              }
              
              print(paste0("N1 = ",N1," - N2 = ",N2))
              
              # UPDATE P ----------------------------------------------------------------
              
              list_p <- update_psigma_cpp(p0_1, sigma_1, CH1, s1, N1, 
                                          K, S_usage, trapsX, trapsY,
                                          sigma_0_1, sd_sigma_1, sigma_sigma1_prop,
                                          a_p, b_p)
              p0_1 <- list_p[1]
              sigma_1 <- list_p[2]
              
              list_p2 <- update_psigma_cpp(p0_2, sigma_2, CH2, s2, N2, 
                                           K, S_usage, trapsX, trapsY,
                                           sigma_0_2, sd_sigma_2, sigma_sigma2_prop,
                                           a_p, b_p)
              p0_2 <- list_p2[1]
              sigma_2 <- list_p2[2]
              
              # list_p <- update_p(p0_1, sigma_1, p0_2, sigma_2,
              #                    s1, s2, CH1, CH2, trapsX, trapsY,
              #                    S_usage, N1, N2,
              #                    sigma_0_1, sd_sigma_1,
              #                    sigma_0_2, sd_sigma_2,
              #                    sigma_p0_1_prop, sigma_p0_2_prop,
              #                    sigma_sigma1_prop, sigma_sigma2_prop,
              #                    a_p, b_p)
              # p0_1 <- list_p$p0_1
              # p0_2 <- list_p$p0_2
              # sigma_1 <- list_p$sigma_1
              # sigma_2 <- list_p$sigma_2
              
              # UPDATE HOME CENTERS -----------------------------------------------
              
              list_s <-  update_s_cpp(s1, s2, theta1, theta2, theta12, beta1, beta2, 
                                      CH1, CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                      Sigma_prop = diag(.005, nrow = 2),
                                      Sigma_newpoint,
                                      lambda_movedPoints1 = 20, 
                                      lambda_movedPoints2 = 20,
                                      lambda_newPoints1 = 0, 
                                      lambda_newPoints2 = 0, 
                                      a1, b1, a2, b2,
                                      allTraps, R_traps, S_area1, S_area2,
                                      S_usage, trapsX, trapsY)
              # s1[1,]
              s1 <- list_s$s1
              s2 <- list_s$s2
              N1 <- list_s$N1
              N2 <- list_s$N2
              
              CH1 <- rbind(CH1[1:D1,],
                           matrix(0, nrow = N1 - D1, ncol = K))
              
              CH2 <- rbind(CH2[1:D2,],
                           matrix(0, nrow = N2 - D2, ncol = K))
              
              
              # UPDATE N --------
              
              list_s <- update_N_cpp(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
                                     CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                                     Sigma_newpoint, lambda_movedPoints1 = 20, lambda_movedPoints2 = 20,
                                     lambda_newPoints1, lambda_newPoints2, 
                                     a1, b1, a2, b2,
                                     allTraps, R_traps,
                                     S_area1, S_area2, 
                                     S_usage, trapsX,
                                     trapsY)
              s1 <- list_s$s1
              s2 <- list_s$s2
              N1 <- list_s$N1
              N2 <- list_s$N2
              mh_ratio1 <- list_s$mh_ratio1
              mh_ratio2 <- list_s$mh_ratio2
              
              CH1 <- rbind(CH1[1:D1,],
                           matrix(0, nrow = N1 - D1, ncol = K))
              
              CH2 <- rbind(CH2[1:D2,],
                           matrix(0, nrow = N2 - D2, ncol = K))
              
              # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
              
              if(updateTheta){
                
                if(iter > iterAfterAdapting){
                  Sigma_n <- cov(params_values[iterToDiscard:(iter-1),c(1,3,5)])
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
                
                if(theta1_star > 0 & theta12_star > 0){
                  
                  ratio <- hastings_ratio_cpp(s1, s2, N1, N2,
                                              x_all[,1,,], 
                                              x_all[,2,,],
                                              c(beta1, beta2,
                                                theta1, theta2, theta12),
                                              c(beta1_star, beta2,
                                                theta1_star, theta2, theta12_star),
                                              N_all,
                                              mu_logtheta, sd_logtheta,
                                              a_beta, b_beta)
                  
                  if(!is.na(ratio) & theta1_star > 0){
                    
                    if(runif(1) < ratio){
                      beta1 <- beta1_star
                      theta1 <- theta1_star
                      theta12 <- theta12_star
                      
                      # list_xall <- update_x_all_foreach(x_all, N_all, 100,
                      #                                   theta1, theta2, theta12,
                      #                                   beta1, beta2, diag(0.05, nrow = 2),
                      #                                   Sigma_newpoint = diag(1, nrow = 2),
                      #                                   allTraps, R_traps,
                      #                                   a1, b1, a2, b2)
                      # x_all <- list_xall$x_all
                      # N_all <- list_xall$N_all
                      
                      x_all1 <- x_all[,1,,]
                      x_all2 <- x_all[,2,,]
                      
                      list_xall <- update_x_all_cpp(x_all1, x_all2, N_all, 100,
                                                    theta1, theta2, theta12,
                                                    beta1, beta2, diag(0.05, nrow = 2),
                                                    Sigma_newpoint = diag(1, nrow = 2),
                                                    a1, b1, a2, b2)
                      x_all[,1,,] <-  list_xall$x_all1
                      x_all[,2,,] <-  list_xall$x_all2
                      N_all <- list_xall$N_all
                      
                    }
                    
                  }
                  
                }
                
                params_values[iter,c(1,3,5)] <- c(beta1, theta1, theta12)
                
              } else {
                
                beta1 <- rgamma(1, a_beta1 + N1, b_beta1 + 1)
                
              }
              
              # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
              
              if(updateTheta){
                
                if(iter > iterAfterAdapting){
                  Sigma_n <- cov(params_values[iterToDiscard:(iter-1),c(2,4,5)])
                  Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
                    beta_proposal * diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
                } else {
                  Sigma_proposal <- diag(c(epsilon_beta, epsilon_theta2, epsilon_theta12), nrow = 3) 
                }
                
                proposed_theta <- mvrnorm(1, c(beta2, theta2, theta12),
                                          Sigma = Sigma_proposal)
                beta2_star <- proposed_theta[1]
                theta2_star <- proposed_theta[2]
                theta12_star <- proposed_theta[3]
                
                if(theta2_star > 0 & theta12_star > 0){
                  
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
                  
                  # print(paste0("mh_ratio_beta2 = ",ratio))
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
                                                    a1, b1, a2, b2)
                      x_all[,1,,] <-  list_xall$x_all1
                      x_all[,2,,] <-  list_xall$x_all2
                      N_all <- list_xall$N_all
                      
                    }
                    
                  }
                  
                }
                
                params_values[iter,c(2,4,5)] <- c(beta2, theta2, theta12)  
                
              } else {
                
                beta2 <- rgamma(1, a_beta2 + N2, b_beta2 + 1)
                
              }
              
              
              
              # beta1 <- rgamma(1, a_beta + N1, b_beta + 1)
              # beta2 <- rgamma(1, a_beta + N2, b_beta + 1)
              
              # WRITE RESULTS ---------------
              
              if(iter > nburn){
                
                # trueIter <- (iter - nburn)/nthin
                trueIter <- iter - nburn
                
                s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
                s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
                p0_iter[trueIter,] <- c(p0_1, p0_2)
                sigma_iter[trueIter,] <- c(sigma_1, sigma_2)
                N_iter[trueIter,] <- c(N1,N2)  
                params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
                
                iter <- iter + 1
              }
              
              
            }
            
            # Write results in MCMC output
            {
              s_output[chain,,,,] <- s_iter 
              sigma_output[chain,,] <- sigma_iter
              N_output[chain,,] <- N_iter
              p0_output[chain,,] <- p0_iter
              params_output[chain,,] <- params_iter
            }
            
          }
          
          # quick diagnostics
          {
            i <- 2
            j <- 1
            dim_s <- 1
            qplot(1:niter, s_iter[,j,i,dim_s])
          }
          
          logparams_ratio[run, idx_n_traps, idx_N_ind, 
                          idx_interRatioStrength, idx_p, 1, ] <- 
            as.numeric(quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975))) 
          logparams_ratio[run, idx_n_traps, idx_N_ind, 
                          idx_interRatioStrength, idx_p, 2, ] <- 
            as.numeric(quantile(log(params_iter[,4] / params_iter[,5]), probs = c(0.025, 0.975))) 
          
        }
        
      }
      
    }
    
  }
  
}

