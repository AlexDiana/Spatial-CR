repetitions <- 10

trueN <- matrix(NA, repetitions, 2)
meanN_inter <- matrix(NA, repetitions, 2)
sdN_inter <- matrix(NA, repetitions, 2)
meanN_nointer <- matrix(NA, repetitions, 2)
sdN_nointer <- matrix(NA, repetitions, 2)

for (repet in 1:repetitions) {
  
  # simulate new data
  {
    usingTrueData <- F
    
    # starting params
    {
      theta1 <- exp(-5)
      theta2 <- exp(-5)
      theta12 <- exp(-10)
      beta1 <- 500
      beta2 <- 500
      
      beta1_true <- beta1
      beta2_true <- beta2
      theta1_true <- theta1
      theta2_true <- theta2
      theta12_true <- theta12
      
      p0_1 <- .005#.005
      p0_2 <- .005#.005
      sigma_1 <- 1 / trapsX_sd
      sigma_2 <- 1 / trapsX_sd
      
      p0_1_true <- p0_1
      p0_2_true <- p0_2
      
      sigma_1_true <- sigma_1
      sigma_2_true <- sigma_2
      
    }
    
    # data simulation
    {
      list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                            Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                            Sigma_newpoint = diag(1, nrow = 2), allTraps, R_traps,
                                            a1, b1, a2, b2)
      s1 <- list_s1s2$data1
      s2 <- list_s1s2$data2
      N1 <- list_s1s2$N1
      N2 <- list_s1s2$N2
      
      s1 <- s1[1:N1,]
      s2 <- s2[1:N2,]
      
      print(N1)
      print(N2)
      # ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") + 
      # geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
      
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
      
      CH_leopards <- CH_1
      CH_tigers <- CH_2
      
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
    }
  }
  
  trueN[repet,] <- c(N1_0, N2_0)
  
  # clean data
  {
    # CLEAN DATA --------------------------------------------------------------
    
    D1 <- nrow(CH_leopards)
    D2 <- nrow(CH_tigers)
    
    K <- ncol(CH_leopards) # number of traps
    
    CH1 <- CH_leopards
    CH2 <- CH_tigers
    
  }
  
  # fit interaction model
  {
    # PRIOR -------------------------------------------------------------------
    
    Nmax <- 1200
    
    # interaction parameter
    {
      a_theta <- 1
      b_theta <- 1
      # mu_gamma <- 1
      # sigmasq_gamma <- 20
      ggplot(data = NULL) + xlim(c(0,0.6)) +
        stat_function(fun = dgamma, args = list(shape = a_theta, rate = b_theta))
    }
    
    # proposal interaction parameters
    {
      epsilon_beta <- 70^2
      epsilon_logtheta1 <- .75
      epsilon_logtheta2 <- .75
      epsilon_logtheta12 <- .25
      # epsilon_beta <- 1000
      # epsilon_logtheta1 <- 1000
      # epsilon_logtheta2 <- 500
      # epsilon_logtheta12 <- .02  
    }
    
    # adaptation parameters
    {
      beta_proposal <- .1
      iterAfterAdapting <- 1000  
    }
    
    # intensity parameter
    {
      meanBeta <- 200
      varBeta <- 1000
      b_beta <- meanBeta / varBeta
      a_beta <- meanBeta * b_beta  
      
      ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,400))
    }
    
    # proposal value
    sigma_prop <- .01
    Sigma_newpoint <- cov(cbind(trapsX, trapsY))
    lambda_newPoints1 <- 7
    lambda_newPoints2 <- 7
    
    # capture probability params
    
    {
      a_p <- 1
      b_p <- 1
      
      # leopards (.5 to 2.5)
      sigma_0_1 <- 1.5 / trapsX_sd
      sd_sigma_1 <- 2 / trapsX_sd
      
      # tigers ( 1.5 to 3.5)
      sigma_0_2 <- 1.5 / trapsY_sd
      sd_sigma_2 <- 2 / trapsY_sd
      
      # ggplot(data = NULL) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1)) + xlim(c(0.01, .5))
      # ggplot(data = NULL) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1)) + xlim(c(0.01, .5))
      
      sigma_p0_1_prop <- .0002
      sigma_p0_2_prop <- .0001
      sigma_sigma1_prop <- .005
      sigma_sigma2_prop <- .01
    }
    
    # points proposal parameters
    {
      probMixture <- c(.95,.95)
    }
    
    # MCMC --------------------------------------------------------------------
    
    pointsToAccept <- 1000#nburn <- 0
    niter <- 10000
    nchain <- 1
    nthin <- 1
    nburn <- 0
    
    # grid fornew points density
    {
      gridLength_x <- seq(min(trapsX) - .05, 
                          max(trapsY) + .05, length.out = 100)
      gridLength_y <- seq(min(trapsY) - .05, 
                          max(trapsY) + .05, length.out = 100)
    }
    
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
      
      # starting values
      # if(nburn != 0){
      
      # home centers
      {
        if(usingTrueData){
          N1 <- 2 * D1
          N2 <- 2 * D2
        } else {
          N1 <- N1_0
          N2 <- N2_0
        }
        
        CH1 <- rbind(CH1[1:D1,],
                     matrix(0, nrow = N1 - D1, ncol = K))
        
        CH2 <- rbind(CH2[1:D2,],
                     matrix(0, nrow = N2 - D2, ncol = K))
        
        if(usingTrueData){
          
          s1 <- matrix(NA, nrow = Nmax, ncol = 2)
          for (i in 1:N1) {
            if(i <= D1){ # assign the home center to the average location of the traps where the indiv. was captured
              s1[i,1] <- mean(trapsX[CH1[i,] != 1])
              s1[i,2] <- mean(trapsY[CH1[i,] != 1])#apply(trapsX[CH1[i,] == 1, 2:3], 2, mean)
              s1[i,] <- mvrnorm(1, mu = s1[i,], diag(.005, nrow = 2))
            } else { # otherwise a random point
              s1[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
            }
          }
          
          s2 <- matrix(NA, nrow = Nmax, ncol = 2)
          for (i in 1:N2) {
            if(i <= D2){ # assign the home center to the average location of the traps where the indiv. was captured
              s2[i,1] <- mean(trapsX[CH2[i,] != 1])
              s2[i,2] <- mean(trapsY[CH2[i,] != 1])
              s2[i,] <- mvrnorm(1, mu = s2[i,], diag(.005, nrow = 2))
            } else { # otherwise a random point
              s2[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
            }
          }
          
        } else {
          
          s1 <- matrix(NA, nrow = Nmax, ncol = 2)
          s2 <- matrix(NA, nrow = Nmax, ncol = 2)
          
          s1[1:N1,] <- s1_0
          s2[1:N2,] <- s2_0
          
        }
        
      }
      
      # capture probability parameter
      {
        if(usingTrueData){
          
          sigma_1 <- sigma_0_1
          
          c_i <- sapply(1:N1, function(i){
            sapply(1:K, function(k){
              exp(- (1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
            })
          })
          
          p0_1 <- rbeta(1, sum(CH1) + 1, sum(mean(S_usage) - CH1) + 1) / mean(c_i)
          
          sigma_2 <- sigma_0_2
          
          c_i <- sapply(1:N2, function(i){
            sapply(1:K, function(k){
              exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
            })
          })
          
          p0_2 <- rbeta(1, sum(CH2) + 1, sum(mean(S_usage) - CH2) + 1) / mean(c_i)
          
        } else {
          p0_1 <- p0_1_true
          sigma_1 <- sigma_1_true
          p0_2 <- p0_2_true
          sigma_2 <- sigma_2_true
        }
        
      }
      
      # point process parameters
      {
        if(usingTrueData){
          theta1 <- .002
          theta2 <- .002
          theta12 <- .002  
          beta1 <- N1
          beta2 <- N2
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
        
        M <- 50
        
        x_all <- array(NA, dim = c(M, 2, Nmax, 2))
        N_all <- matrix(NA, nrow = M, ncol = 2)
        
        for (m in 1:M) {
          
          print(m)
          
          list_sims <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2,
                                                niter = 2000,
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
      
      # new proposal parameters
      {
        
        numAcceptedPoints <- c(0,0)
        acceptedPoints1 <- matrix(NA, pointsToAccept, 2)
        acceptedPoints2 <- matrix(NA, pointsToAccept, 2)
        
        # compute normalizing constant for uniform
        {
          M <- 10000
          
          normConstUnif <- mean(sapply(1:M, function(i){
            x <- c(runif(1, a1, b1), runif(1, a2, b2))
            checkPointIsInRegionTraps(x, allTraps, R_traps)
          }))
          
          S_area <- (b1 - a1) * (b2 - a2) * normConstUnif
          
          normConst1 <- c(normConstUnif)
          normConst2 <- c(normConstUnif)
          
          mixtureMeans1 <- matrix(0, nrow = 0, ncol = 2)
          mixtureSd1 <- matrix(0, nrow = 0, ncol = 2)
          mixtureMeans2 <- matrix(0, nrow = 0, ncol = 2)
          mixtureSd2 <- matrix(0, nrow = 0, ncol = 2)
          
        }
        
      }
      
      # }
      
      # ready starting values
      # if(nburn == 0) {
      # load("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/CR Model/newstarintg_values.rda")
      # }
      
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
      
      iter <- 1
      outOfBurnIn <- F
      mixture1Updated <- F
      mixture2Updated <- F
      
      while(iter <= niter){
        
        if(!outOfBurnIn){
          print(paste0("Chain = ",chain," / AcceptedPoints1 = ",numAcceptedPoints[1]," / AcceptedPoints2 = ",numAcceptedPoints[2])) 
        } else if(((iter - nburn)/nthin) %% 5 == 0){
          print(paste0("Chain = ",chain," / Iteration = ",(iter - nburn)/nthin))
        }
        
        print(paste0("N1 = ",N1," - N2 = ",N2))
        # print(paste0("p0_1 = ",p0_1," - p0_2 = ",p0_2))
        # print(paste0("sigma_1 = ",sigma_1," - sigma_2 = ",sigma_2))
        
        # UPDATE P ----------------------------------------------------------------
        
        list_p <- update_p(p0_1, sigma_1, p0_2, sigma_2,
                           s1, s2, CH1, CH2, trapsX, trapsY,
                           S_usage, N1, N2,
                           sigma_0_1, sd_sigma_1,
                           sigma_0_2, sd_sigma_2,
                           sigma_p0_1_prop, sigma_p0_2_prop,
                           sigma_sigma1_prop, sigma_sigma2_prop,
                           a_p, b_p)
        p0_1 <- list_p$p0_1
        p0_2 <- list_p$p0_2
        sigma_1 <- list_p$sigma_1
        sigma_2 <- list_p$sigma_2
        
        # UPDATE HOME CENTERS -----------------------------------------------
        
        if(outOfBurnIn){
          
          # list_s <-  update_s_cpp(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
          #                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
          #                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
          #                         allTraps, R_traps,
          #                         a1, b1, a2, b2)
          list_s <- update_s_cpp_new(s1, s2, theta1, theta2, theta12, beta1, beta2, CH1,
                                     CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                                     Sigma_newpoint, lambda_newPoints1, lambda_newPoints2, probMixture[1],
                                     probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
                                     normConst1, normConst2, allTraps, R_traps, S_area, S_usage, trapsX,
                                     trapsY, a1, b1, a2, b2)
        } else {
          list_s <-  update_s_cpp_withacceptance(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
                                                 CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                                 S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
                                                 allTraps, R_traps,
                                                 a1, b1, a2, b2)
        }
        
        s1 <- list_s$s1
        s2 <- list_s$s2
        N1 <- list_s$N1
        N2 <- list_s$N2
        
        CH1 <- rbind(CH1[1:D1,],
                     matrix(0, nrow = N1 - D1, ncol = K))
        
        CH2 <- rbind(CH2[1:D2,],
                     matrix(0, nrow = N2 - D2, ncol = K))
        
        if(!outOfBurnIn){
          if(numAcceptedPoints[1] < pointsToAccept & list_s$isPoint1Accepted){
            acceptedPoints1[numAcceptedPoints[1] + 1,] <- list_s$acceptedPoint1
            numAcceptedPoints[1] <- numAcceptedPoints[1] + 1
          }
          
          if(numAcceptedPoints[2] < pointsToAccept & list_s$isPoint2Accepted){
            acceptedPoints2[numAcceptedPoints[2] + 1,] <- list_s$acceptedPoint2
            numAcceptedPoints[2] <- numAcceptedPoints[2] + 1
          }  
        }
        
        # UPDATE BETA 1, THETA 1 AND THETA 12 ---------------------------------------------------------
        
        # # if(iter > iterAfterAdapting){
        # # Sigma_n <- cov(params_values[1:(iter-1),c(1,3,5)])
        # # Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
        # # beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
        # # } else {
        # Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta1, epsilon_logtheta12), nrow = 3)
        # # }
        # 
        # proposed_theta <- mvrnorm(1, c(beta1, log(theta1), log(theta12)),
        #                           Sigma = Sigma_proposal)
        # beta1_star <- proposed_theta[1]
        # theta1_star <- exp(proposed_theta[2])
        # theta12_star <- exp(proposed_theta[3])
        # 
        # ratio <- hastings_ratio(s1, s2, N1, N2,
        #                         x_all, N_all,
        #                         beta1, beta2,
        #                         beta1_star, beta2,
        #                         theta1, theta2, theta12,
        #                         theta1_star, theta2, theta12_star,
        #                         # mu_gamma, sigmasq_gamma,
        #                         a_theta, b_theta,
        #                         a_beta, b_beta)
        # # print(ratio)
        # if(!is.na(ratio) & theta1_star > 0){
        #   
        #   if(runif(1) < ratio){
        #     beta1 <- beta1_star
        #     theta1 <- theta1_star
        #     theta12 <- theta12_star
        #     
        #     list_xall <- update_x_all_foreach(x_all, N_all, 100,
        #                                       theta1, theta2, theta12,
        #                                       beta1, beta2, diag(0.05, nrow = 2),
        #                                       Sigma_newpoint = diag(1, nrow = 2),
        #                                       allTraps, R_traps,
        #                                       a1, b1, a2, b2)
        #     x_all <- list_xall$x_all
        #     N_all <- list_xall$N_all
        #     # print("accepted 1")
        #   }
        #   
        # }
        
        # params_values[iter,c(1,3,5)] <- c(beta1, log(theta1), log(theta12))
        
        # UPDATE BETA 2, THETA 2 AND THETA 12 ---------------------------------------------------------
        
        # # if(iter > iterAfterAdapting){
        # # Sigma_n <- cov(params_values[1:(iter-1),c(2,4,5)])
        # # Sigma_proposal <- (1 - beta_proposal) * (2.38) * Sigma_n / 3 +
        # # beta_proposal * (((0.1)^2) / 3) * diag(1, nrow = 3)
        # # } else {
        # Sigma_proposal <- diag(c(epsilon_beta, epsilon_logtheta2, epsilon_logtheta12), nrow = 3)
        # # }
        # 
        # proposed_theta <- mvrnorm(1, c(beta2, log(theta2), log(theta12)),
        #                           Sigma = Sigma_proposal)
        # beta2_star <- proposed_theta[1]
        # theta2_star <- exp(proposed_theta[2])
        # theta12_star <- exp(proposed_theta[3])
        # 
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
        # if(!is.na(ratio) & theta2_star > 0){
        #   
        #   if(runif(1) < ratio){
        #     beta2 <- beta2_star
        #     theta2 <- theta2_star
        #     theta12 <- theta12_star
        #     
        #     list_xall <- update_x_all_foreach(x_all, N_all, 100,
        #                                       theta1, theta2, theta12,
        #                                       beta1, beta2, diag(0.05, nrow = 2),
        #                                       Sigma_newpoint = diag(1, nrow = 2),
        #                                       allTraps, R_traps,
        #                                       a1, b1, a2, b2)
        #     x_all <- list_xall$x_all
        #     N_all <- list_xall$N_all
        #     # print("accepted 2")
        #     # acceptances_theta[iter,1] <- 1
        #   }
        #   
        # }
        # 
        # params_values[iter,c(2,4,5)] <- c(beta2, log(theta2), log(theta12))
        
        # WRITE RESULTS ---------------
        
        outOfBurnIn <- numAcceptedPoints[1] >= pointsToAccept & numAcceptedPoints[2] >= pointsToAccept
        
        if(outOfBurnIn){
          
          # trueIter <- (iter - nburn)/nthin
          trueIter <- iter
          
          # s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
          # s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
          p0_iter[trueIter,] <- c(p0_1, p0_2)
          sigma_iter[trueIter,] <- c(sigma_1, sigma_2)
          N_iter[trueIter,] <- c(N1,N2)  
          params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
          # papangelou_density_iter[trueIter,,,] <- computeNewPointsDensity(gridLength_x, gridLength_y,
          # s1, s2, N1, N2,
          # theta1,  theta2,  theta12,
          # beta1,  beta2)
          
          iter <- iter + 1
        }
        
        # update proposals
        {
          if(outOfBurnIn & !mixture1Updated){
            
            numCenters_all <- 5:6
            kmeans_totss <- rep(length(numCenters_all))
            
            for (l in seq_along(numCenters_all)) {
              
              numCenters <- numCenters_all[l]
              
              kmeansfit <- kmeans(acceptedPoints1, centers = numCenters)
              
              kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints1)) 
            }
            
            best_k <- numCenters_all[which.min(kmeans_totss)]
            kmeansfit <- kmeans(acceptedPoints1, centers = best_k)
            
            mixtureMeans1 <- matrix(NA, nrow = best_k, ncol = 2)
            mixtureSd1 <- matrix(NA, nrow = best_k, ncol = 2)
            normConstNormals1 <- rep(NA, best_k)
            
            for (i in 1:best_k) {
              mixtureMeans1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, mean)
              mixtureSd1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, sd)
              
              normConstNormals1[i] <- mean(sapply(1:M, function(j){
                x <- c(rnorm(1, mixtureMeans1[i,1], mixtureSd1[i,1]), 
                       rnorm(1, mixtureMeans1[i,2], mixtureSd1[i,2]))
                checkPointIsInRegionTraps(x, allTraps, R_traps)
              }))
            }
            
            normConst1 <- c(normConstUnif, normConstNormals1)
            
            mixture1Updated <- T
          }
          
          if(outOfBurnIn & !mixture2Updated){
            
            numCenters_all <- 5:6
            kmeans_totss <- rep(length(numCenters_all))
            
            for (l in seq_along(numCenters_all)) {
              
              numCenters <- numCenters_all[l]
              
              kmeansfit <- kmeans(acceptedPoints2, centers = numCenters)
              
              kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints2)) 
            }
            
            best_k <- numCenters_all[which.min(kmeans_totss)]
            kmeansfit <- kmeans(acceptedPoints2, centers = best_k)
            
            mixtureMeans2 <- matrix(NA, nrow = best_k, ncol = 2)
            mixtureSd2 <- matrix(NA, nrow = best_k, ncol = 2)
            normConstNormals2 <- rep(NA, best_k)
            
            for (i in 1:best_k) {
              mixtureMeans2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, mean)
              mixtureSd2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, sd)
              
              normConstNormals2[i] <- mean(sapply(1:M, function(j){
                x <- c(rnorm(1, mixtureMeans2[i,1], mixtureSd2[i,1]), 
                       rnorm(1, mixtureMeans2[i,2], mixtureSd2[i,2]))
                checkPointIsInRegionTraps(x, allTraps, R_traps)
              }))
            }
            
            normConst2 <- c(normConstUnif, normConstNormals2)
            
            mixture2Updated <- T
          }    
        }
        
        if(iter %% 250 == 0){
          plotCurrent1 <- qplot(1:iter, params_iter[1:iter,1], geom =  "line")
          plotCurrent2 <- qplot(1:iter, log(params_iter[1:iter,3]), geom =  "line")
          plotCurrent3 <- qplot(1:iter, log(params_iter[1:iter,5]), geom =  "line")
          ggsave(filename = paste0("plot_param1_iter",iter,".jpg"), plotCurrent1)
          ggsave(filename = paste0("plot_param2_iter",iter,".jpg"), plotCurrent2)
          ggsave(filename = paste0("plot_param3_iter",iter,".jpg"), plotCurrent3)
        }
        
      }
      
      # Write results in MCMC output
      {
        sigma_output[chain,,] <- sigma_iter
        s_output[chain,,,,] <- s_iter#array(NA, dim = c(nchain, niter, 2, Nmax, 2))
        N_output[chain,,] <- N_iter#array(NA, dim = c(nchain, niter, 2))
        p0_output[chain,,] <- p0_iter#array(NA, dim = c(nchain, niter, 2))
        params_output[chain,,] <- params_iter#array(NA, dim = c(nchain, niter, 5))
      }
      
    }
    
    meanN_inter[repet,] <- apply(N_iter,2,mean)
    sdN_inter[repet,] <- apply(N_iter,2,sd)
  }
  
  # fit no interaction model
  {
    # PRIOR -------------------------------------------------------------------
    
    Nmax <- 1200
    
    # proposal interaction parameters
    {
      epsilon_beta <- 70^2
    }
    
    # intensity parameter
    {
      meanBeta <- N1
      varBeta <- 10000
      b_beta <- meanBeta / varBeta
      a_beta <- meanBeta * b_beta  
      
      # ggplot(data = NULL) + stat_function(fun = dgamma, args = list(shape = a_beta, rate = b_beta)) + xlim(c(0,400))
    }
    
    # proposal value
    sigma_prop <- .01
    Sigma_newpoint <- cov(cbind(trapsX, trapsY))
    lambda_newPoints1 <- 11
    lambda_newPoints2 <- 11
    
    # capture probability params
    {
      a_p <- 1
      b_p <- 1
      
      # leopards (.5 to 2.5)
      sigma_0_1 <- 1.5 / traps_meansd
      sd_sigma_1 <- 2 / traps_meansd
      
      # tigers ( 1.5 to 3.5)
      sigma_0_2 <- 2.5 / traps_meansd
      sd_sigma_2 <- 2 / traps_meansd
      
      # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_1, sd = sd_sigma_1))
      # ggplot(data = NULL, aes(x = c(-2,5))) + stat_function(fun = dnorm, args = list(mean = sigma_0_2, sd = sd_sigma_2))
      
      sigma_p0_1_prop <- .0002
      sigma_sigma1_prop <- .005
      
      sigma_p0_2_prop <- .0001
      sigma_sigma2_prop <- .01
    }
    
    # points proposal parameters
    {
      probMixture <- c(.95,.95)
    }
    
    # MCMC --------------------------------------------------------------------
    
    pointsToAccept <- 1000#nburn <- 0
    niter <- 20000
    nchain <- 1
    nthin <- 1
    nburn <- 0
    
    # grid fornew points density
    {
      gridLength_x <- seq(min(trapsX) - .05, 
                          max(trapsY) + .05, length.out = 100)
      gridLength_y <- seq(min(trapsY) - .05, 
                          max(trapsY) + .05, length.out = 100)
    }
    
    # variables for output
    {
      sigma_output <- array(NA, dim = c(nchain, niter, 2))
      # s_output <- array(NA, dim = c(nchain, niter, 2, Nmax, 2))
      N_output <- array(NA, dim = c(nchain, niter, 2))
      p0_output <- array(NA, dim = c(nchain, niter, 2))
      params_output <- array(NA, dim = c(nchain, niter, 5))
      
      # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
      # params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
      
      # acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
    }
    
    for(chain in 1:nchain) {
      
      # starting values
      # if(nburn != 0){
      
      # home centers
      {
        if(usingTrueData){
          N1 <- 2 * D1
          N2 <- 2 * D2
        } else {
          N1 <- N1_0
          N2 <- N2_0
        }
        
        CH1 <- rbind(CH1[1:D1,],
                     matrix(0, nrow = N1 - D1, ncol = K))
        
        CH2 <- rbind(CH2[1:D2,],
                     matrix(0, nrow = N2 - D2, ncol = K))
        
        if(usingTrueData){
          
          s1 <- matrix(NA, nrow = Nmax, ncol = 2)
          for (i in 1:N1) {
            if(i <= D1){ # assign the home center to the average location of the traps where the indiv. was captured
              s1[i,1] <- mean(trapsX[CH1[i,] != 1])
              s1[i,2] <- mean(trapsY[CH1[i,] != 1])#apply(trapsX[CH1[i,] == 1, 2:3], 2, mean)
              s1[i,] <- mvrnorm(1, mu = s1[i,], diag(.005, nrow = 2))
            } else { # otherwise a random point
              s1[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
            }
          }
          
          s2 <- matrix(NA, nrow = Nmax, ncol = 2)
          for (i in 1:N2) {
            if(i <= D2){ # assign the home center to the average location of the traps where the indiv. was captured
              s2[i,1] <- mean(trapsX[CH2[i,] != 1])
              s2[i,2] <- mean(trapsY[CH2[i,] != 1])
              s2[i,] <- mvrnorm(1, mu = s2[i,], diag(.005, nrow = 2))
            } else { # otherwise a random point
              s2[i,] <- proposeNewPoint(a1, b1, a2, b2, allTraps, R_traps)
            }
          }
          
        } else {
          
          s1 <- matrix(NA, nrow = Nmax, ncol = 2)
          s2 <- matrix(NA, nrow = Nmax, ncol = 2)
          
          s1[1:N1,] <- s1_0
          s2[1:N2,] <- s2_0
          
        }
        
      }
      
      # capture probability parameter
      {
        if(usingTrueData){
          
          sigma_1 <- sigma_0_1
          
          c_i <- sapply(1:N1, function(i){
            sapply(1:K, function(k){
              exp(- (1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
            })
          })
          
          p0_1 <- rbeta(1, sum(CH1) + 1, sum(mean(S_usage) - CH1) + 1) / mean(c_i)
          
          sigma_2 <- sigma_0_2
          
          c_i <- sapply(1:N2, function(i){
            sapply(1:K, function(k){
              exp(- (1 / (2 * sigma_2^2)) * ((s2[i,1] - trapsX[k])^2 + (s2[i,2] - trapsY[k])^2))
            })
          })
          
          p0_2 <- rbeta(1, sum(CH2) + 1, sum(mean(S_usage) - CH2) + 1) / mean(c_i)
          
        } else {
          p0_1 <- p0_1_true
          sigma_1 <- sigma_1_true
          p0_2 <- p0_2_true
          sigma_2 <- sigma_2_true
        }
        
      }
      
      # point process parameters
      {
        if(usingTrueData){
          theta1 <- .002
          theta2 <- .002
          theta12 <- .002  
          beta1 <- N1
          beta2 <- N2
        } else {
          theta1 <- theta1_true
          theta2 <- theta2_true
          theta12 <- theta12_true
          beta1 <- N1
          beta2 <- N2
        }
        
      }
      
      # new proposal parameters
      {
        
        numAcceptedPoints <- c(0,0)
        acceptedPoints1 <- matrix(NA, pointsToAccept, 2)
        acceptedPoints2 <- matrix(NA, pointsToAccept, 2)
        
        # compute normalizing constant for uniform
        {
          M <- 10000
          
          normConstUnif <- mean(sapply(1:M, function(i){
            x <- c(runif(1, a1, b1), runif(1, a2, b2))
            checkPointIsInRegionTraps(x, allTraps, R_traps)
          }))
          
          S_area <- (b1 - a1) * (b2 - a2) * normConstUnif
          
          normConst1 <- c(normConstUnif)
          normConst2 <- c(normConstUnif)
          
          mixtureMeans1 <- matrix(0, nrow = 0, ncol = 2)
          mixtureSd1 <- matrix(0, nrow = 0, ncol = 2)
          mixtureMeans2 <- matrix(0, nrow = 0, ncol = 2)
          mixtureSd2 <- matrix(0, nrow = 0, ncol = 2)
          
        }
        
      }
      
      # }
      
      # output variables
      {
        sigma_iter <- array(NA, dim = c(niter, 2))
        # s_iter <- array(NA, dim = c(niter, 2, Nmax, 2))
        N_iter <- matrix(NA, nrow = niter, ncol = 2)
        p0_iter <- array(NA, dim = c(niter, 2))
        params_iter <- matrix(NA, nrow = niter, ncol = 5)
        
        # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
        params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
        
        acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
        theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
      }
      
      iter <- 1
      outOfBurnIn <- F
      mixture1Updated <- F
      mixture2Updated <- F
      
      while(iter <= niter){
        
        if(!outOfBurnIn){
          print(paste0("Chain = ",chain," / AcceptedPoints1 = ",numAcceptedPoints[1]," / AcceptedPoints2 = ",numAcceptedPoints[2])) 
        } else if(((iter - nburn)/nthin) %% 5 == 0){
          print(paste0("Chain = ",chain," / Iteration = ",(iter - nburn)/nthin))
        }
        
        print(paste0("N1 = ",N1," - N2 = ",N2))
        
        # UPDATE P ----------------------------------------------------------------
        
        list_p <- update_p(p0_1, sigma_1, p0_2, sigma_2,
                           s1, s2, CH1, CH2, trapsX, trapsY,
                           S_usage, N1, N2,
                           sigma_0_1, sd_sigma_1,
                           sigma_0_2, sd_sigma_2,
                           sigma_p0_1_prop, sigma_p0_2_prop,
                           sigma_sigma1_prop, sigma_sigma2_prop,
                           a_p, b_p)
        p0_1 <- list_p$p0_1
        p0_2 <- list_p$p0_2
        sigma_1 <- list_p$sigma_1
        sigma_2 <- list_p$sigma_2
        
        # UPDATE HOME CENTERS -----------------------------------------------
        
        # if(iter < nburn){
        # list_s <-  update_s_cpp(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
        #                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
        #                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
        #                         allTraps, R_traps,
        #                         a1, b1, a2, b2)  
        
        if(outOfBurnIn){
          
          # list_s <-  update_s_cpp(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
          #                         CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
          #                         S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
          #                         allTraps, R_traps,
          #                         a1, b1, a2, b2)
          list_s <- update_s_cpp_new(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
                                     CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2, Sigma_prop = diag(.1, nrow = 2),
                                     Sigma_newpoint, lambda_newPoints1, lambda_newPoints2, probMixture[1],
                                     probMixture[2], mixtureMeans1, mixtureMeans2, mixtureSd1, mixtureSd2,
                                     normConst1, normConst2, allTraps, R_traps, S_area, S_usage, trapsX,
                                     trapsY, a1, b1, a2, b2)
        } else {
          list_s <-  update_s_cpp_withacceptance(s1, s2, exp(-100), exp(-100), exp(-100), beta1, beta2, CH1,
                                                 CH2, N1, N2, D1, D2, p0_1, sigma_1, p0_2, sigma_2,
                                                 S_usage, trapsX, trapsY, Sigma_prop = diag(.1, nrow = 2), Sigma_newpoint,
                                                 allTraps, R_traps,
                                                 a1, b1, a2, b2)
        }
        
        s1 <- list_s$s1
        s2 <- list_s$s2
        N1 <- list_s$N1
        N2 <- list_s$N2
        
        CH1 <- rbind(CH1[1:D1,],
                     matrix(0, nrow = N1 - D1, ncol = K))
        
        CH2 <- rbind(CH2[1:D2,],
                     matrix(0, nrow = N2 - D2, ncol = K))
        
        if(!outOfBurnIn){
          if(numAcceptedPoints[1] < pointsToAccept & list_s$isPoint1Accepted){
            acceptedPoints1[numAcceptedPoints[1] + 1,] <- list_s$acceptedPoint1
            numAcceptedPoints[1] <- numAcceptedPoints[1] + 1
          }
          
          if(numAcceptedPoints[2] < pointsToAccept & list_s$isPoint2Accepted){
            acceptedPoints2[numAcceptedPoints[2] + 1,] <- list_s$acceptedPoint2
            numAcceptedPoints[2] <- numAcceptedPoints[2] + 1
          }  
        }
        
        # UPDATE BETA ---------
        
        beta1 <- rgamma(1, a_beta + N1, b_beta + 1)
        beta2 <- rgamma(1, a_beta + N2, b_beta + 1)
        
        # WRITE RESULTS ---------------
        
        outOfBurnIn <- numAcceptedPoints[1] >= pointsToAccept & numAcceptedPoints[2] >= pointsToAccept
        
        if(outOfBurnIn){
          
          # trueIter <- (iter - nburn)/nthin
          trueIter <- iter
          
          # s_iter[trueIter,1,1:N1,] <- s1[1:N1,]
          # s_iter[trueIter,2,1:N2,] <- s2[1:N2,]
          p0_iter[trueIter,] <- c(p0_1, p0_2)
          sigma_iter[trueIter,] <- c(sigma_1, sigma_2)
          N_iter[trueIter,] <- c(N1,N2)  
          params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
          # papangelou_density_iter[trueIter,,,] <- computeNewPointsDensity(gridLength_x, gridLength_y,
          # s1, s2, N1, N2,
          # theta1,  theta2,  theta12,
          # beta1,  beta2)
          
          iter <- iter + 1
        }
        
        # update proposals
        {
          if(outOfBurnIn & !mixture1Updated){
            
            numCenters_all <- 5:30
            kmeans_totss <- rep(length(numCenters_all))
            
            for (l in seq_along(numCenters_all)) {
              
              numCenters <- numCenters_all[l]
              
              kmeansfit <- kmeans(acceptedPoints1, centers = numCenters)
              
              kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints1)) 
            }
            
            best_k <- numCenters_all[which.min(kmeans_totss)]
            kmeansfit <- kmeans(acceptedPoints1, centers = best_k)
            
            mixtureMeans1 <- matrix(NA, nrow = best_k, ncol = 2)
            mixtureSd1 <- matrix(NA, nrow = best_k, ncol = 2)
            normConstNormals1 <- rep(NA, best_k)
            
            for (i in 1:best_k) {
              mixtureMeans1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, mean)
              mixtureSd1[i,] <- apply(acceptedPoints1[which(kmeansfit$cluster == i),,drop = F], 2, sd)
              
              normConstNormals1[i] <- mean(sapply(1:M, function(j){
                x <- c(rnorm(1, mixtureMeans1[i,1], mixtureSd1[i,1]), 
                       rnorm(1, mixtureMeans1[i,2], mixtureSd1[i,2]))
                checkPointIsInRegionTraps(x, allTraps, R_traps)
              }))
            }
            
            normConst1 <- c(normConstUnif, normConstNormals1)
            
            mixture1Updated <- T
          }
          
          if(outOfBurnIn & !mixture2Updated){
            
            numCenters_all <- 5:30
            kmeans_totss <- rep(length(numCenters_all))
            
            for (l in seq_along(numCenters_all)) {
              
              numCenters <- numCenters_all[l]
              
              kmeansfit <- kmeans(acceptedPoints2, centers = numCenters)
              
              kmeans_totss[l] <- kmeansfit$tot.withinss + .5*numCenters*2*log(nrow(acceptedPoints2)) 
            }
            
            best_k <- numCenters_all[which.min(kmeans_totss)]
            kmeansfit <- kmeans(acceptedPoints2, centers = best_k)
            
            mixtureMeans2 <- matrix(NA, nrow = best_k, ncol = 2)
            mixtureSd2 <- matrix(NA, nrow = best_k, ncol = 2)
            normConstNormals2 <- rep(NA, best_k)
            
            for (i in 1:best_k) {
              mixtureMeans2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, mean)
              mixtureSd2[i,] <- apply(acceptedPoints2[which(kmeansfit$cluster == i),,drop = F], 2, sd)
              
              normConstNormals2[i] <- mean(sapply(1:M, function(j){
                x <- c(rnorm(1, mixtureMeans2[i,1], mixtureSd2[i,1]), 
                       rnorm(1, mixtureMeans2[i,2], mixtureSd2[i,2]))
                checkPointIsInRegionTraps(x, allTraps, R_traps)
              }))
            }
            
            normConst2 <- c(normConstUnif, normConstNormals2)
            
            mixture2Updated <- T
          }    
        }
        
        # intermediate plots
        {
          if(iter %% 250 == 0){
            # plotCurrent1 <- qplot(1:iter, params_iter[1:iter,1], geom =  "line")
            # plotCurrent2 <- qplot(1:iter, log(params_iter[1:iter,3]), geom =  "line")
            # plotCurrent3 <- qplot(1:iter, log(params_iter[1:iter,5]), geom =  "line")
            # ggsave(filename = paste0("plot_param1_iter",iter,".jpg"), plotCurrent1)
            # ggsave(filename = paste0("plot_param2_iter",iter,".jpg"), plotCurrent2)
            # ggsave(filename = paste0("plot_param3_iter",iter,".jpg"), plotCurrent3)
          }  
        }
        
      }
      
      # Write results in MCMC output
      {
        sigma_output[chain,,] <- sigma_iter
        N_output[chain,,] <- N_iter
        p0_output[chain,,] <- p0_iter
        params_output[chain,,] <- params_iter
      }
      
    }
    
    
    
    meanN_nointer[repet,] <- apply(N_iter,2,mean)
    sdN_nointer[repet,] <- apply(N_iter,2,sd)
  }
  
}

save(trueN, meanN_inter, sdN_inter, meanN_nointer, sdN_nointer, file = "allResults.rda")


