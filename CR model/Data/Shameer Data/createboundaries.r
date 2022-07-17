load("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/Data/Shameer Data/data.rda")

traps$X <- (traps$X - mean(traps$X)) / sd(traps$X)
traps$Y <- (traps$Y - mean(traps$Y)) / sd(traps$Y)

pointsBoundaries <- matrix(NA, nrow = 10, ncol = 2)
pointsBoundaries[1,] <- c(-2.6, -0.3)
pointsBoundaries[2,] <- c(-.39, 2.05)
pointsBoundaries[3,] <- c(.85, -2.65)
pointsBoundaries[4,] <- c(1.9, .45)
pointsBoundaries[5,] <- c(1.5, 1.9)
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
  
ggplot() + geom_point(data = traps, aes(x = X, y = Y, color = ACTIVE_LEOPARD)) + 
  theme(legend.title = element_blank()) + scale_color_discrete() + theme_bw() +
  geom_point(data = NULL, aes(x = pointsBoundaries[,1], y = pointsBoundaries[,2]), 
             color = "black", size = 2) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[2,1])), 
                fun = line1, xlim = c(pointsBoundaries[1,1],pointsBoundaries[2,1])) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[1,1],pointsBoundaries[3,1])), 
                fun = line2, xlim = c(pointsBoundaries[1,1],pointsBoundaries[3,1])) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[3,1])), 
                fun = line3, xlim = c(pointsBoundaries[4,1],pointsBoundaries[3,1])) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[2,1],pointsBoundaries[5,1])), 
                fun = line4, xlim = c(pointsBoundaries[2,1],pointsBoundaries[5,1])) + 
  stat_function(data = NULL, aes(x = c(pointsBoundaries[4,1],pointsBoundaries[5,1])), 
                fun = line5, xlim = c(pointsBoundaries[4,1],pointsBoundaries[5,1]))

# DIVIDE INTO A GRID


