library(ggplot2); library(Rcpp); library; library(MASS); library(ggplot2); library(reshape2)
library(foreach); library(doParallel)

#See how many cores we have
ncl<- detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
# stopCluster(cl)

setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
sourceCpp("code.cpp")

a <- 0
b <- 1

theta1 <- exp(-20)
theta2 <- exp(-20)
theta12 <- exp(-20)
beta1 <- 50
beta2 <- 100

niter <- 100

N_output <- matrix(NA, nrow = niter, ncol = 2)

for (iter in 1:niter) {
  
  print(iter)
  
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 500, 100, 2,
                                        a, b)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  # print(N1)
  # print(N2)
  
  N_output[iter,] <- c(N1, N2)
}

apply(N_output, 2, mean)
apply(N_output, 2, var)
