s_output_all <- apply(s_output,c(3,4,5),c)
s_output_all <- s_output_all[,,1:max(N_output),]

N_output_all <- apply(N_output,3,c)

library(MASS)

grid_points <- 250
density_output <- matrix(NA, nrow = niter*nchain, grid_points * grid_points)
idx_species <- 1
for (iter in 1:(niter*nchain)) {
  
  if(iter %% 1000 == 0){
    print(iter)  
  }
  
  N_iter <- N_output_all[iter,]
  kd2fit <- kde2d(x = s_output_all[iter,idx_species,
                                   1:N_iter[idx_species],1], 
                  y = s_output_all[iter,idx_species,
                                   1:N_iter[idx_species],2], 
                  n = grid_points,
                  lims = tiger_polyboundaries)
  
  density_output[iter,] <- as.vector(kd2fit$z)
  
}

density_masses <- apply(density_output, 2, mean)

allGrid <- expand.grid(kd2fit$x, kd2fit$y)

pointInside <- apply(allGrid, 1, function(x){
  checkPointIsInRegionPolygonsAndTraps(x, leop_polycoord, leop_polyhole,
                                       leop_polystart, allTraps, R_traps)
  # checkPointIsInRegionPolygonsAndTraps(x, tiger_polycoord, tiger_polyhole,
                                       # tiger_polystart, allTraps, R_traps)
})

newGridStep <- allGrid[2,1] - allGrid[1,1]

density_masses <- density_masses[pointInside]
allGrid <- allGrid[pointInside,]

npoints <- nrow(allGrid)
# X_tilde_rescaled <- X_tilde_star
# X_tilde_rescaled[,1] <- meanDataEast + sdDataEast * X_tilde_rescaled[,1]
# X_tilde_rescaled[,2] <- meanDataNorth + sdDataNorth * X_tilde_rescaled[,2]
# 
# View(X_tilde_rescaled)
# 
# npoints <- nrow(X_tilde_rescaled)
datapoly <- data.frame(
  id = rep(1:npoints, each = 4),
  x = rep(allGrid[,1], each =  4),
  y = rep(allGrid[,2], each =  4),
  value = rep(density_masses, each =  4)
)

sdDataEast <- 1
sdDataNorth <- 1

variations_x <- rep(c(newGridStep * sdDataEast, -newGridStep * sdDataEast,
                      -newGridStep * sdDataEast, newGridStep * sdDataEast), times = npoints)
variations_y <- rep(c(-newGridStep * sdDataNorth, -newGridStep * sdDataNorth,
                      newGridStep * sdDataNorth, newGridStep * sdDataNorth), times = npoints)

datapoly$x <- datapoly$x + variations_x
datapoly$y <- datapoly$y + variations_y

breaks_x <- seq(-3, 3, by = .75)
breaks_y <- seq(-2.5, 2.5, by = .75)
labels_x <- round(trapsX_mean + traps_meansd * breaks_x, 0)
labels_y <- round(trapsY_mean + traps_meansd * breaks_y, 0)

p <- ggplot() 

if(idx_species == 1){
  p <- p + 
    geom_polygon(data = NULL, aes(x = leop_polycoord[(leop_polystart[2] + 1):(leop_polystart[3] + 0),1],
                                  y = leop_polycoord[(leop_polystart[2] + 1):(leop_polystart[3] + 0),2]),
                 fill = "grey") +
    geom_polygon(data = NULL, aes(x = leop_polycoord[(leop_polystart[3] + 1):(leop_polystart[4] + 0),1],
                                  y = leop_polycoord[(leop_polystart[3] + 1):(leop_polystart[4] + 0),2]),
                 fill = "grey") +
    geom_polygon(data = NULL, aes(x = leop_polycoord[(leop_polystart[4] + 1):(leop_polystart[5]),1],
                                  y = leop_polycoord[(leop_polystart[4] + 1):(leop_polystart[5]),2]),
                 fill = "grey")  +
  geom_polygon(data = NULL, aes(x = leop_polycoord[1:(leop_polystart[2] + 0),1],
                                y = leop_polycoord[1:(leop_polystart[2] + 0),2]),
               fill = "blue")
} else {
  p <- p + 
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[4] + 1):(tiger_polystart[5] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[4] + 1):(tiger_polystart[5] + 0),2]),
                 fill = "grey")   +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[5] + 1):(tiger_polystart[6] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[5] + 1):(tiger_polystart[6] + 0),2]),
                 fill = "grey") +   
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[6] + 1):(tiger_polystart[7] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[6] + 1):(tiger_polystart[7] + 0),2]),
                 fill = "grey") +   
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[7] + 1):(tiger_polystart[8] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[7] + 1):(tiger_polystart[8] + 0),2]),
                 fill = "grey") +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[8] + 1):(tiger_polystart[9] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[8] + 1):(tiger_polystart[9] + 0),2]),
                 fill = "grey") +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[9] + 1):(tiger_polystart[10] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[9] + 1):(tiger_polystart[10] + 0),2]),
                 fill = "grey") +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[10] + 1):(tiger_polystart[11] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[10] + 1):(tiger_polystart[11] + 0),2]),
                 fill = "grey")  +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[11] + 1):(tiger_polystart[12] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[11] + 1):(tiger_polystart[12] + 0),2]),
                 fill = "grey")  +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[1:(tiger_polystart[2] + 1),1],
                                  y = tiger_polycoord[1:(tiger_polystart[2] + 1),2]),
                 fill = "white") +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[2] + 2):(tiger_polystart[3] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[2] + 2):(tiger_polystart[3] + 0),2]),
                 fill = "white") +
    geom_polygon(data = NULL, aes(x = tiger_polycoord[(tiger_polystart[3] + 2):(tiger_polystart[4] + 0),1],
                                  y = tiger_polycoord[(tiger_polystart[3] + 2):(tiger_polystart[4] + 0),2]),
                 fill = "blue") 
    
  p
}

p <- p +
  geom_polygon(data = datapoly, aes(x = x, y = y,
                                    fill = value, group = id),
               alpha = .7) +
  scale_fill_gradientn(limits = range(datapoly$value),
                       colors = c("yellow","red")) + 
  geom_point(data = NULL, aes(x = allTraps[,1],
                              y = allTraps[,2]), color = "black", size = .125, alpha = .8) + 
  theme_grey(base_size = 11, base_family = "", 
             base_line_size = 11 / 22, base_rect_size = 11 / 22) %+replace% 
  theme(legend.position = "none", 
        panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.border = element_rect(fill = NA, 
                                                                                  colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85", 
                                        colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                      colour = NA), complete = TRUE) + 
  scale_x_continuous(breaks = breaks_x, labels = labels_x, name = "Kilometers") +
  scale_y_continuous(breaks = breaks_y, labels = labels_y, name = "Kilometers") 


p

#
