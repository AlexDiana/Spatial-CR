library(scales)


# PLOT FUNCTION -----------------------------------------------------------

require(scales)
# p + scale_y_continuous(trans = log2_trans(),
#                        breaks = trans_breaks("log2", function(x) 2^x),
#                        labels = trans_format("log2", math_format(2^.x)))

# logm_trans = function() trans_new("logm", function(x) x / 100, function(x) 100 * x)

scaleFactor <- 1
scaleFactor2 <- 2

f1 <- function(x) {
  ifelse(x > - 1, 
         x / scaleFactor2,
         -log(- scaleFactor * x) + log(scaleFactor) - 1 / scaleFactor2
  )
}

f2 <- function(x) {
  ifelse(x > - 1 / scaleFactor2, 
         x * scaleFactor2, 
         - 1 / scaleFactor * exp(-(x - log(scaleFactor) + 1 / scaleFactor2)) )
}

x_grid <- seq(-0.001, -10, length.out = 100)
y_grid <- Vectorize(f1)(x_grid)
ggplot(data = NULL) + xlim(c(-.001, -3)) + stat_function(fun = f2)

logm_trans = function() trans_new("logm", f1, 
                                  f2)


# FIRST INTERACTION REGIME -----

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_003_N_300_int_20.RData")
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_003_N_300_int_20_again.RData")
CI_params_p03_N300_int20 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_003_N_300_int_10.RData")
CI_params_p03_N300_int10 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_03_N_300_int_20.RData")
CI_params_p3_N300_int20 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_03_N_300_int_10.RData")
CI_params_p3_N300_int10 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_003_N_200_int_20.RData")
CI_params_p03_N200_int20 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_003_N_200_int_10.RData")
CI_params_p03_N200_int10 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_03_N_200_int_20.RData")
load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_03_N_200_int_20_again.RData")
CI_params_p3_N200_int20 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

load("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Simulation Results/p_03_N_200_int_10.RData")
CI_params_p3_N200_int10 <- quantile(log(params_iter[,3] / params_iter[,5]), probs = c(0.025, 0.975, 0.5))

# N - 300 / INT = 20 --------------

CI_params <- cbind(CI_params_p3_N300_int20,
                   CI_params_p03_N300_int20)
# CI_params <- cbind(CI_params_p03_N300_int20,
#                    CI_params_p03_N300_int10,
#                    CI_params_p3_N300_int20,
#                    CI_params_p3_N300_int10,
#                    CI_params_p03_N200_int20,
#                    CI_params_p03_N200_int10,
#                    CI_params_p3_N200_int20,
#                    CI_params_p3_N200_int10)


plotgamma_1 <- ggplot(data = NULL, aes(x = reorder(factor(c(".3",".03")),c(1,2)), 
                                       y = CI_params[3,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_hline(aes(yintercept = log(1 / 20)), color = "blue") + 
  scale_x_discrete(labels = c("p = .3","p = .03"),
                   name = "Capture Probabilities") +
  scale_y_continuous(name = expression(log(gamma[1] / gamma[1][2])),
                     trans = "logm",
                     breaks = c(0, -1 , -2, -5, -10, -20, -40),
                     limits = c(-60, 1)) +  
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

ggsave(filename = "N300_int20.jpeg",  plotgamma_1,
       width = 6.36, height = 5.25, units = "in")

# N - 300 / INT = 10 --------------

CI_params <- cbind(CI_params_p3_N300_int10,
                   CI_params_p03_N300_int10)
# CI_params <- cbind(CI_params_p03_N300_int20,
#                    CI_params_p03_N300_int10,
#                    CI_params_p3_N300_int20,
#                    CI_params_p3_N300_int10,
#                    CI_params_p03_N200_int20,
#                    CI_params_p03_N200_int10,
#                    CI_params_p3_N200_int20,
#                    CI_params_p3_N200_int10)

plotgamma_1 <- ggplot(data = NULL, aes(x = reorder(factor(c(".3",".03")),c(1,2)), 
                                       y = CI_params[3,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_hline(aes(yintercept = log(1 / 10)), color = "blue") + 
  scale_x_discrete(labels = c("p = .3","p = .03"),
                   name = "Capture Probabilities") +
  #scale_y_continuous(name = expression(gamma[1][2] / gamma[1]), limits = c(0, 1.5)) +  
  scale_y_continuous(name = expression(log(gamma[1] / gamma[1][2])),
                     trans = "logm",
                     breaks = c(0, -1 , -2, -5, -10, -20, -40),
                     limits = c(-60, 1)) +  
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

ggsave(filename = "N300_int10.jpeg",  plotgamma_1,
       width = 6.36, height = 5.25, units = "in")

# N - 200 / INT = 20 --------------

CI_params <- cbind(CI_params_p3_N200_int20,
                   CI_params_p03_N200_int20)
# CI_params <- cbind(CI_params_p03_N300_int20,
#                    CI_params_p03_N300_int10,
#                    CI_params_p3_N300_int20,
#                    CI_params_p3_N300_int10,
#                    CI_params_p03_N200_int20,
#                    CI_params_p03_N200_int10,
#                    CI_params_p3_N200_int20,
#                    CI_params_p3_N200_int10)

plotgamma_1 <- ggplot(data = NULL, aes(x = reorder(factor(c(".3",".03")),c(1,2)), 
                                       y = CI_params[3,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_hline(aes(yintercept = log(1 / 20)), color = "blue") + 
  scale_x_discrete(labels = c("p = .3","p = .03"),
                   name = "Capture Probabilities") +
  #scale_y_continuous(name = expression(gamma[1][2] / gamma[1]), limits = c(0, 2)) +  
  scale_y_continuous(name = expression(log(gamma[1] / gamma[1][2])),
                     trans = "logm",
                     breaks = c(0, -1 , -2, -5, -10, -20, -40),
                     limits = c(-60, 1)) +   
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

ggsave(filename = "N200_int20.jpeg",  plotgamma_1,
       width = 6.36, height = 5.25, units = "in")

# N - 200 / INT = 10 --------------

CI_params <- cbind(CI_params_p3_N200_int10,
                   CI_params_p03_N200_int10)
# CI_params <- cbind(CI_params_p03_N300_int20,
#                    CI_params_p03_N300_int10,
#                    CI_params_p3_N300_int20,
#                    CI_params_p3_N300_int10,
#                    CI_params_p03_N200_int20,
#                    CI_params_p03_N200_int10,
#                    CI_params_p3_N200_int20,
#                    CI_params_p3_N200_int10)

plotgamma_1 <- ggplot(data = NULL, aes(x = reorder(factor(c(".3",".03")),c(1,2)), 
                                       y = CI_params[3,], 
                                       ymin = CI_params[1,], 
                                       ymax = CI_params[2,])) + 
  geom_point(stat = "identity", size = 3) + geom_errorbar(size = 1, width = .5) + 
  geom_hline(aes(yintercept = 0), color = "red") + 
  geom_hline(aes(yintercept = log( 1/ 10)), color = "blue") + 
  scale_x_discrete(labels = c("p = .3","p = .03"),
                   name = "Capture Probabilities") +
  scale_y_continuous(name = expression(log(gamma[1] / gamma[1][2])),
                     trans = "logm",
                     breaks = c(0, -1 , -2, -5, -10, -20, -40),
                     limits = c(-60, 1)) +  
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

ggsave(filename = "N200_int10.jpeg",  plotgamma_1,
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
