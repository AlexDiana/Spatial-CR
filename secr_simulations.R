# read data
library(secr)
setwd("/cluster/home/osr/ad625/Spatial CR/SECRdata")

CH <- read.capthist(captfile = "tig.txt", trapfile = "trapnew.txt", fmt = "trapID", 
                    detector = "proximity")
trapsUsage <- attr(attr(CH, "traps"),"usage")
nOcc <- ncol(trapsUsage)

thab<-rgdal::readOGR(dsn = "/cluster/home/osr/ad625/Spatial CR/SECRdata/TigerMask", layer = "ctrmask")
smask <- make.mask(traps(CH), 
                   spacing = 500,
                   type = 'trapbuffer',
                   # type = 'polygon',
                   buffer = 8000,
                   poly = thab)
plot(smask)
plot(traps(CH), add = T)

# fit on real data

fit <- secr.fit(CH2,  
                mask=smask, 
                # buffer = 8000,
                trace = T)
region.N(fit)

# fit on simulated data

attrCH <- attributes(CH)

# starting params
{
  theta1 <- exp(-10)
  theta2 <- exp(-10)
  theta12 <- exp(-10)
  beta1 <- 350
  beta2 <- 350
  
  beta1_true <- beta1
  beta2_true <- beta2
  theta1_true <- theta1
  theta2_true <- theta2
  theta12_true <- theta12
  
  p0_1 <- .05#.005
  p0_2 <- .02#.005
  sigma_1 <- 2.5 / traps_meansd
  sigma_2 <- 2.5 / trapsX_sd
  
  p0_1_true <- p0_1
  p0_2_true <- p0_2
  
  sigma_1_true <- sigma_1
  sigma_2_true <- sigma_2
  
}

# data simulation
{
  list_s1s2 <- simulate_bivsoftcore_cpp(theta1, theta2, theta12, beta1, beta2, 2000, 
                                        Sigma_prop = diag(0.01, nrow = 2), 2000, beta1, 
                                        Sigma_newpoint = diag(1, nrow = 2), 
                                        leop_polycoord, leop_polyhole, leop_polystart, 
                                        leop_polyboundaries, tiger_polycoord, 
                                        tiger_polyhole, tiger_polystart, 
                                        tiger_polyboundaries)
  s1 <- list_s1s2$data1
  s2 <- list_s1s2$data2
  N1 <- list_s1s2$N1
  N2 <- list_s1s2$N2
  
  s1 <- s1[1:N1,]
  s2 <- s2[1:N2,]
  
  print(N1)
  print(N2)
  ggplot() +  geom_point(data = NULL, aes(x = s1[,1], y = s1[,2]), color = "red") +
    geom_point(data = NULL, aes(x = s2[,1], y = s2[,2]), color = "black") +
    geom_point(data = NULL, aes(x = trapsX,
                                y = trapsY), color = "green")
  
  
  K <- nrow(traps)
  
  CH_new <- array(NA, dim = c(N1, nOcc, K))
  for (i in 1:N1) {
    for (k in 1:K) {
      for (l in 1:nOcc) {
        if(trapsUsage[k,l] == 1){
          p <- p0_1 * exp(-(1 / (2 * sigma_1^2)) * ((s1[i,1] - trapsX[k])^2 + (s1[i,2] - trapsY[k])^2))
          CH_new[i,l,k] <- rbinom(1, 1, p)    
        } else {
          CH_new[i,l,k] <- 0
        }
      }
      
    }
  }
  indexesCaughtIndividuals1 <- which(apply(CH_new,1,sum) != 0)
  indexesUncaughtIndividuals1 <- which(apply(CH_new,1,sum) == 0)
  D1 <- length(indexesCaughtIndividuals1)
  CH_new <- CH_new[indexesCaughtIndividuals1,,]
  
  CH <- CH_new  
  
  attrCH$dim[1] <- D1
  attrCH$dimnames[[1]] <- as.character(1:D1)
  attributes(CH) <- attrCH
}

fit <- secr.fit(CH,  
                mask=smask, 
                # buffer = 8000,
                model = list(g0 ~ 1),
                trace = T)
region.N(fit)
