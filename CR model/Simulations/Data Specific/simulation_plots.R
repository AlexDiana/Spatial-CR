# FIRST INTERACTION REGIME -----

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p100_beta_500_theta1_15_theta12_4.RData")
CI_params_p100_theta1_15 <- quantile(params_iter[,5] - params_iter[,3], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p5_beta_500_theta1_15_theta12_4.RData")
CI_params_p5_theta1_15 <- quantile(params_iter[,5] - params_iter[,3], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p1_beta_500_theta1_15_theta12_4.RData")
CI_params_p1_theta1_15 <- quantile(params_iter[,5] - params_iter[,3], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p01_beta_500_theta1_15_theta12_4.RData")
CI_params_p01_theta1_15 <- quantile(params_iter[,5] - params_iter[,3], probs = c(0.025, 0.975, 0.5))

CI_params <- cbind(CI_params_p100_theta1_15,
                   CI_params_p5_theta1_15,
                   CI_params_p1_theta1_15,
                   CI_params_p01_theta1_15)

plotgamma_1 <- ggplot(data = NULL, aes(x = reorder(factor(c(".99",".5",".1",".01")),c(1,2,3,4)), 
                                     y = CI_params[3,], 
                                     ymin = CI_params[1,], 
                                     ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  scale_x_discrete(labels = c("p = .99","p = .5","p = .1","p = .01"),
                   name = "Capture Probabilities") +
  scale_y_continuous(name = expression(gamma[1][2] - gamma[1]), limits = c(-0.3, 0.15)) +  
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
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"))

plotgamma_1

setwd("C:/Users/Alex/Dropbox/PhD/Latex/Papers/Spatial model/Img/Simulation")

ggsave(filename = "theta1_15_theta2_4.jpeg",  plotgamma_1,
       width = 6.36, height = 5.25, units = "in")

# SECOND INTERACTION REGIME -----

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p100_beta_600_theta1_4_theta12_15.RData")
CI_params_p100_theta1_15 <- quantile(params_iter[,3] - params_iter[,5], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p5_beta_600_theta1_4_theta12_15.RData")
CI_params_p5_theta1_15 <- quantile(params_iter[,3] - params_iter[,5], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p1_beta_600_theta1_4_theta12_15.RData")
CI_params_p1_theta1_15 <- quantile(params_iter[,3] - params_iter[,5], probs = c(0.025, 0.975, 0.5))
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/sim_p01_beta_600_theta1_4_theta12_15.RData")
CI_params_p01_theta1_15 <- quantile(params_iter[,3] - params_iter[,5], probs = c(0.025, 0.975, 0.5))

CI_params <- cbind(CI_params_p100_theta1_15,
                   CI_params_p5_theta1_15,
                   CI_params_p1_theta1_15,
                   CI_params_p01_theta1_15)

plotgamma_2 <- ggplot(data = NULL, aes(x = reorder(factor(c(".99",".5",".1",".01")),c(1,2,3,4)), 
                                     y = CI_params[3,], 
                                     ymin = CI_params[1,], 
                                     ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  scale_x_discrete(labels = c("p = .99","p = .5","p = .1","p = .01"),
                   name = "Capture Probabilities") +
  scale_y_continuous(name = expression(gamma[1] - gamma[1][2]), limits = c(-0.3, 0.15)) + 
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
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"))

plotgamma_2

setwd("C:/Users/Alex/Dropbox/PhD/Latex/Papers/Spatial model/Img/Simulation")

ggsave(filename = "theta1_4_theta2_15.jpeg",  plotgamma_2,
       width = 6.36, height = 5.25, units = "in")
