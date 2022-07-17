# FUNCTIONS ---------------------------------------------------------------

setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial")
# 
library(Rcpp); library(RcppArmadillo); library(ggplot2)
# 
sourceCpp("code.cpp")

# CONDITIONAL SIMULATION ---------------

N1 <- 100
N2 <- 100

theta1 <- .01#0000
theta2 <- .01#0000
theta12 <- .1#0000

niter <- 5000

list_conditional <- simulate_conditional_bivsoftcore_cpp(N1, N2, theta1, theta2, theta12, niter,
                                                         Sigma_prop = diag((.01)^2, nrow = 2),
                                                         a = 0, b = 1)

data1 <- list_conditional$data1
data2 <- list_conditional$data2

data1 <- data1[1:N1,]
data2 <- data2[1:N2,]

ggplot() + geom_point(data = NULL, aes(x = data1[,1], y = data1[,2]), color = "grey") + 
  geom_point(data = NULL, aes(x = data2[,1], y = data2[,2]), color = "black") + theme_bw() + xlab("") + ylab("")



