library(MASS); library(reshape2); library(ggplot2)
library(foreach); library(doParallel)
# library(SpatialFunctionsCR);
library(Rcpp); library(RcppArmadillo)
library(here)

#See how many cores we have
# ncl<- detectCores()
# cl <- makeCluster(ncl)
# registerDoParallel(cl)
# stopCluster(cl)

sourceCpp("Parameter Estimation/code.cpp")
source('Parameter Estimation/functions.r', echo=TRUE)

# SIMULATED DATA ----------------------------------------------------------

usingTrueData <- F

# starting params
{
  theta1 <- exp(-100)
  theta2 <- exp(-100)
  theta12 <- exp(-7.5) # exp(-6)
  beta1 <- N1 <- 150
  beta2 <- N2 <- 150
  
  beta1_true <- beta1
  beta2_true <- beta2
  theta1_true <- theta1
  theta2_true <- theta2
  theta12_true <- theta12
  
  a1 <- 0
  b1 <- 1
  a2 <- 0
  b2 <- 1
}

# data simulation
{
  list_s1s2 <- simulate_cond_bivsoftcore_cpp(theta1, theta2, theta12, N1, N2, 5000, 
                                             Sigma_prop = diag(0.01, nrow = 2), N1, beta1, 
                                             Sigma_newpoint = diag(1, nrow = 2), 
                                             a1, b1, a2, b2)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
 
  ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") +
  geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black")
  
}

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
  params_output <- array(NA, dim = c(nchain, niter, 5))
}

for(chain in 1:nchain) {
 
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
    
    params_iter <- matrix(NA, nrow = niter, ncol = 5)
    
    # papangelou_density_iter <- array(0, dim = c(niter, 2, length(gridLength_x), length(gridLength_y)))  
    params_values <- matrix(NA, nrow = nburn + nthin*niter, ncol = 5)
    
    acceptances_theta <- matrix(0, nrow = nburn + niter * nthin, ncol = 3)
    theta_all <- matrix(NA, nrow = nburn + niter * nthin, ncol = 5)
  }
  
  for(iter in 1:(nburn + niter)){
    
    if(iter < nburn){
      print(paste0("Chain = ",chain," / Iteration = ", iter)) 
    } else {
      print(paste0("Chain = ",chain," / Iteration = ", iter - nburn))
    }
    
    print(paste0("beta1 = ",beta1," - beta2 = ",beta2))
    
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
        
        # ratio <- hastings_ratio_2(s1, s2, N1, N2,
        #                           x_all[,1,,], 
        #                           x_all[,2,,],
        #                           N_all,
        #                         beta1, beta2,
        #                         beta1, beta2_star,
        #                         theta1, theta2, theta12,
        #                         theta1, theta2_star, theta12_star,
        #                         mu_logtheta, sd_logtheta,
        #                         # a_theta, b_theta,
        #                         a_beta, b_beta)
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
            
            # list_xall <- update_x_all_foreach(x_all, N_all, 100,
            #                                   theta1, theta2, theta12,
            #                                   beta1, beta2, diag(0.05, nrow = 2),
            #                                   Sigma_newpoint = diag(1, nrow = 2),
            #                                   allTraps, R_traps,
            #                                   a1, b1, a2, b2)
            # x_all <- list_xall$x_all
            # N_all <- list_xall$N_all
            
            # print("accepted 2")
            # acceptances_theta[iter,1] <- 1
          }
          
        }
        
      }
      
      params_values[iter,c(2,4,5)] <- c(beta2, theta2, theta12)  
      
    } else {
      
      beta2 <- rgamma(1, a_beta2 + N2, b_beta2 + 1)
      
    }
    
    
    # WRITE RESULTS ---------------
    
    if(iter > nburn){
      
      # trueIter <- (iter - nburn)/nthin
      trueIter <- iter - nburn
       
      params_iter[trueIter,] <- c(beta1, beta2, theta1, theta2, theta12)
      
    }
    
  }
  
  # Write results in MCMC output
  {
    params_output[chain,,] <- params_iter
  }
  
}

qplot(1:niter, log(params_iter[,3] / params_iter[,5])) + geom_hline(aes(yintercept = log(theta1_true / theta12_true)))
quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975))
qplot(1:niter, log(params_iter[,4] / params_iter[,5])) + geom_hline(aes(yintercept = log(theta2_true / theta12_true)))
quantile(log(params_iter[,4] / params_iter[,5]), probs = c(0.025, 0.975))

# CHAIN DIAGNOSTICS -------------------------------------------------------

# plot chains
{
  plotVar <- function(variable_output){
    
    nchain <- nrow(variable_output)
    niter <- ncol(variable_output)
    
    # beta_psi_output <- matrix(beta_psi_output[,,1], nrow = nchain, ncol = niter)
    
    variable_output_long <- melt(variable_output)
    
    ggplot() + geom_line(data = variable_output_long, aes(x = Var2, y = value, group = Var1, 
                                                          color = factor(Var1))) + 
      xlab("Iterations") + ylab("Value") +
      theme(plot.title = element_text(hjust = 0.5, size = 17),
            axis.title = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold", hjust = 1),
            axis.line = element_line(colour="black", size=0.15),
            # panel.grid.minor = element_line(colour="grey", size=0.15),
            panel.grid.major = element_line(colour="grey", size=0.15),
            panel.background = element_rect(fill = "white", color = "black"),
            legend.position = "none") 
    
  }
  
  
}

# diagnostics p
{
  plotVar(p0_output[,,1,drop = F])
  plotVar(p0_output[,,2])
  plotVar(sigma_output[,,1])
  plotVar(sigma_output[,,2])
}

# diagnostics N
{
  plotVar(N_output[,,1]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
  plotVar(N_output[,,2]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
}

# diagnostics theta
{
  plotVar(params_output[,,1]) 
  plotVar(params_output[,,2]) 
  plotVar(log(params_output[,,3]))
  plotVar(log(params_output[,,4]))
  plotVar(log(params_output[,,5]))
}
# DIAGNOSTICS -------------------------------------------------------------

# diagnostics p
{
  plot(1:niter, p0_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = p0_1_true))
  plot(1:niter, p0_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = p0_2_true))
  plot(1:niter, sigma_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = sigma_1_true))
  plot(1:niter, sigma_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = sigma_2_true))
}

# diagnostics N
{
  plot(1:niter, N_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = N1_0))
  plot(1:niter, N_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = N2_0))
}

# diagnostics params
{
  plot(1:niter, params_iter[,1]) #+ geom_hline(data = NULL, aes(yintercept = beta1_true))
  plot(1:niter, params_iter[,2]) #+ geom_hline(data = NULL, aes(yintercept = beta2_true))
  plot(1:niter, params_iter[,3]) #+ geom_hline(data = NULL, aes(yintercept = theta1_true))
  plot(1:niter, params_iter[,4]) #+ geom_hline(data = NULL, aes(yintercept = theta2_true))
  plot(1:niter, params_iter[,5]) #+ geom_hline(data = NULL, aes(yintercept = theta12_true))
}

params_iter_long <- melt(params_iter)

CI_params <- apply(params_iter, 2, function(x){
  quantile(x, probs = c(.05,.5,.95), na.rm = T)
})
c(D1,D2)
intParamsPlot <- ggplot(data = NULL, aes(x = c("theta1","theta2","theta12"), y = CI_params[2,3:5], 
                                         ymin = CI_params[1,3:5], 
                                         ymax = CI_params[3,3:5])) + 
  geom_point(stat = "identity", size = 4, shape = 1) +
  geom_errorbar() + 
  # geom_point(aes(x = c("theta1","theta2","theta12"), y = c(theta1_true, theta2_true, theta12_true)), shape = 4, size = 4) + 
  scale_x_discrete(name = "Parameter", labels = c(expression(theta[1]),expression(theta[12]),expression(theta[2]))) +
  scale_y_continuous(name = "", limits = c(0, 0.15)) +
  geom_text(data = NULL, aes(x = 2, y = .145), 
            label = expression(paste(D[1]," = ",32, " - ", D[2]," = ",37)),
            size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold", angle = 0, hjust = 1),
        axis.line = element_line(colour="black", size=0.15),
        panel.grid.major = element_line(colour="grey", size=0.15),
        panel.background = element_rect(fill = "white", color = "black")) 
intParamsPlot
setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Simulation")

ggsave(filename = paste0("beta = ",beta1_true, "- p0 = ",p0_true,".jpeg"), intParamsPlot,
       width = 6.36, height = 5.25, units = "in")


qplot(1:niter, params_iter[,1])
qplot(1:niter, params_iter[,2])
qplot(1:niter, params_iter[,5])

s1 <- s1[1:N1,]
s2 <- s2[1:N2,]

qplot(1:niter, p0_iter)

qplot(1:niter, N_iter[,1])
qplot(1:niter, N_iter[,2])

qplot(1:niter, s_iter[,1,1,1])
qplot(1:niter, s_iter[,1,21,2])
qplot(1:niter, s_iter[,2,11,2])

ggplot(data = params_iter, aes(x = c))

# PLOTS -------------------------------------------------------------------

# plots for gamma
{
  CI_params <- apply(params_iter, 2, function(x){
    quantile(x, probs = c(.05,.5,.95))
  })
  
  plotgamma <- ggplot(data = NULL, aes(x = reorder(factor(c("theta1","theta2","theta12")),c(1,3,2)), y = CI_params[2,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[3,])) + 
    geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
    scale_x_discrete(labels = c(expression(theta[1]),expression(theta[1][2]),expression(theta[2])),
                     name = "Interaction parameters") +
    scale_y_continuous(name = "") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13, face = "bold"),
          axis.text = element_text(size = 13, face = "bold"))
  
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "gamma.jpeg",  plotgamma,
         width = 6.36, height = 5.25, units = "in")
}

# plots for N
{
  plotN <- ggplot(data = NULL, aes(x = N_iter[,1])) + 
    geom_bar(aes(y = (..count..)/sum(..count..)), 
             fill = "cornsilk", color = "black", stat = "count", width = .2, size = .5) +
    scale_x_continuous(breaks = sort(unique(N_iter[,1])), 
                       labels =  round(sort(unique(N_iter[,1])) /  925,3),
                       minor_breaks = sort(unique(N_iter[,1])),   
                       name = expression(paste("Number of individuals per ",km^2))) +
    scale_y_continuous(name = "Density") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13, face = "bold"))
  plotN
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "Popsize1.jpeg",  plotN,
         width = 6.36, height = 5.25, units = "in")
  
  plotN2 <- ggplot(data = NULL, aes(x = N_iter[,2])) + 
    geom_bar(aes(y = (..count..)/sum(..count..)), 
             fill = "cornsilk", color = "black", stat = "count", width = .2, size = .5) +
    scale_x_continuous(breaks = sort(unique(N_iter[,2])), 
                       labels =  round(sort(unique(N_iter[,2])) /  925,3),
                       minor_breaks = sort(unique(N_iter[,2])),   
                       name = expression(paste("Number of individuals per ",km^2))) +
    scale_y_continuous(name = "Density") + 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(size = rel(0.5)), 
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), 
          legend.key = element_rect(fill = "white", 
                                    colour = NA),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13, face = "bold"))
  plotN2
  setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")
  
  ggsave(filename = "Popsize2.jpeg",  plotN2,
         width = 6.36, height = 5.25, units = "in")
  
  
}


# CREATE BOUNDARIES -------------------------------------------------------

pointsBoundaries <- matrix(NA, nrow = 10, ncol = 2)
pointsBoundaries[1,] <- c(-2.4, 0.5)
pointsBoundaries[2,] <- c(-1.3, 2.65)
pointsBoundaries[3,] <- c(.85, -2.65)
pointsBoundaries[4,] <- c(2.2, -1.35)
pointsBoundaries[5,] <- c(1.8, .6)
pointsBoundaries <- pointsBoundaries[complete.cases(pointsBoundaries),]

coefsLine <- function(x1, y1, x2, y2){
  
  m <- (y2 - y1) / (x2 - x1)
  
  q <- y1 - m * x1
  
  list("m" = m, "q" = q)
}

line1 <- function(x){
  
  m <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[2,1], pointsBoundaries[2,2])$m
  
  q <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[2,1], pointsBoundaries[2,2])$q
  
  return(m * x + q)
}

line2 <- function(x){
  
  m <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$m
  
  q <- coefsLine(pointsBoundaries[1,1], pointsBoundaries[1,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$q
  
  return(m * x + q)
}

line3 <- function(x){
  
  m <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$m
  
  q <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[3,1], pointsBoundaries[3,2])$q
  
  return(m * x + q)
}

line4 <- function(x){
  
  m <- coefsLine(pointsBoundaries[2,1], pointsBoundaries[2,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$m
  
  q <- coefsLine(pointsBoundaries[2,1], pointsBoundaries[2,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$q
  
  return(m * x + q)
}

line5 <- function(x){
  
  m <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$m
  
  q <- coefsLine(pointsBoundaries[4,1], pointsBoundaries[4,2], 
                 pointsBoundaries[5,1], pointsBoundaries[5,2])$q
  
  return(m * x + q)
}

{
  ggplot() +
    geom_point(data = NULL, aes(x = traps$x, y = traps$y)) + 
    scale_x_continuous(name = "Longitude") + 
    scale_y_continuous(name = "Latitude") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                  fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                  fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
    stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                  fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                  fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
    stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                  fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  
  
}




# HOME RANGES -------------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

sd_s1 <- apply(s_iter[,1,1:D1,], c(2,3), sd)
sd_s2 <- apply(s_iter[,2,1:D2,], c(2,3), sd)

sd_s1_mean <- apply(sd_s1, 1, mean)
sd_s2_mean <- apply(sd_s2, 1, mean)

ggplot() + geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2], size = as.numeric(4 * sd_s1_mean)), color = "red", alpha = .3) + 
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2], size = 4 * sd_s2_mean), color = "black", alpha = .3) + 
  scale_size_continuous() +
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  



# DENSITY PLOT ------------------------------------------------------------

N1 <- N_iter[niter,1]
N2 <- N_iter[niter,2]
s1_mean <- apply(s_iter[,1,1:N1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:N2,], c(2,3), mean)

ggplot() + geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]), color = "red") + 
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]), color = "black") 

gridLength_x <- seq(min(as.vector(s_iter[,,,1]), na.rm = T)  - .05, 
                    max(as.vector(s_iter[,,,1]), na.rm = T) + .05, length.out = 25)
gridLength_y <- seq(min(as.vector(s_iter[,,,2]), na.rm = T) - .05, 
                    max(as.vector(s_iter[,,,2]), na.rm = T) + .05, length.out = 25)
pointsInGrid <- array(0, dim = c(niter, 2, length(gridLength_x) - 1, length(gridLength_y) - 1))

for (iter in 1:niter) {
  print(iter)
  N1 <- N_iter[niter,1]
  for (i in 1:N1) {
    cellX <- findInterval(s_iter[iter,1,i,1], gridLength_x)
    cellY <- findInterval(s_iter[iter,1,i,2], gridLength_y)
    pointsInGrid[iter, 1, cellX, cellY] <- pointsInGrid[iter, 1, cellX, cellY] + 1
  }
  N2 <- N_iter[niter,2]
  for (i in 1:N2) {
    cellX <- findInterval(s_iter[iter,2,i,1], gridLength_x)
    cellY <- findInterval(s_iter[iter,2,i,2], gridLength_y)
    pointsInGrid[iter, 2, cellX, cellY] <- pointsInGrid[iter, 2, cellX, cellY] + 1
  }
}

binwidth <- gridLength_x[2] - gridLength_x[1]
midPoints_x <- (gridLength_x[-1] + gridLength_x[-length(gridLength_x)]) / 2 
midPoints_y <- (gridLength_y[-1] + gridLength_y[-length(gridLength_y)]) / 2 
allGrid <- expand.grid(midPoints_x, midPoints_y)
meanDensityGrid <- apply(pointsInGrid, c(2,3,4), mean)
meanDensityGrid_vec <- as.vector(meanDensityGrid[1,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensityGrid_vec), 
             size = 100 * binwidth, shape = 15) + 
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  


# PAPANGELOU DENSITY ------------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

papangelou_density_mean <- apply(papangelou_density_iter, c(2,3,4), mean)

allGrid <- expand.grid(gridLength_x, gridLength_y)
binwidth <- gridLength_x[2] - gridLength_x[1]

meanDensity1Grid_vec <- as.vector(papangelou_density_mean[1,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensity1Grid_vec), 
             size = 100 * binwidth, shape = 15) + 
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 50 * binwidth, color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 50 * binwidth, color = "red") +
  scale_x_continuous(name = "Longitude") + 
  scale_color_continuous(low = "black", high = "white") +
  scale_y_continuous(name = "Latitude") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "blue") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "blue") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "blue")  

meanDensity2Grid_vec <- as.vector(papangelou_density_mean[2,,])
ggplot() +
  geom_point(data = NULL, aes(x = allGrid[,1], y = allGrid[,2], color = meanDensity2Grid_vec), 
             size = 100 * binwidth, shape = 15) + 
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 50 * binwidth, color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 50 * binwidth, color = "red") +
  scale_x_continuous(name = "Longitude") + 
  scale_y_continuous(name = "Latitude") + 
  scale_color_continuous(low = "black", high = "white") +
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "red") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "red") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "red")  

# HOME CENTERS  -----------------------------------------------------------

s1_mean <- apply(s_iter[,1,1:D1,], c(2,3), mean)
s2_mean <- apply(s_iter[,2,1:D2,], c(2,3), mean)

s1_sd <- apply(s_iter[,1,1:D1,], c(2,3), sd)
s2_sd <- apply(s_iter[,2,1:D2,], c(2,3), sd)

binwidth <- gridLength_x[2] - gridLength_x[1]

meanTrapsX <- 0
sdTrapsX <-  0.1035411
meanTrapsY <- 0
sdTrapsY <- 0.07364998

plotHomeRange <- ggplot() +
  geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2]),
             size = 35 * binwidth,  color = "blue") +
  geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2]),
             size = 35 * binwidth,  color = "red") +
  geom_point(data = NULL, aes(x = trapsX, y = trapsY), size = 1, shape = 3) +
  # geom_point(data = NULL, aes(x = s1_mean[,1], y = s1_mean[,2], size = 4 * s1_sd[,1]),
  # shape = 1, color = "blue") +
  # geom_point(data = NULL, aes(x = s2_mean[,1], y = s2_mean[,2], size = 4 * s2_sd[,1]),
  # shape = 1, color = "red") +
  scale_x_continuous(name = "Longitude", breaks = -2:2, labels = round((-2:2)*sdTrapsX + meanTrapsX,2)) + 
  scale_y_continuous(name = "Latitude", breaks = -2:2, labels = round((-2:2)*sdTrapsY + meanTrapsY,2)) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1]), size = 2, color = "black") +  
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1]), size = 2, color = "black") + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]), size = 2, color = "black") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"),
        panel.grid.major = element_line(colour="grey", size=0.015),
        panel.background = element_rect(fill = "white", color = "black")) 

plotHomeRange
setwd("C:/Users/alexd/Dropbox/PhD/Latex/Papers/Spatial model/Img/Application")

ggsave(filename = "homerangesplot.jpeg", plotHomeRange,
       width = 5, height = 5, units = "in")

