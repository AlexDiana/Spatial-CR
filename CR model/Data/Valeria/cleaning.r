setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/CR Model/Data/Valeria")

captureHistory_Ocelot_MagMedio <- read.table("CaptureHistoryMagMedioOcelot.txt", header = T)
captureHistory_Jaguar_MagMedio <- read.csv("Jaguar Data Mag Medio.csv")

captureHistory_Ocelot_Aurora <- read.csv("CaptureHistoryAuroraOcelot.csv")
captureHistory_Jaguar_Aurora <- read.csv("Jaguar Data Aurora.csv")

trapinfo_MagMedio <- read.csv("TrapInfoMagMedio.csv")
trapinfo_Aurora <- read.csv("TrapInfoAuroraOcelot.csv")

library(ggplot2)

ggplot(data = trapinfo_MagMedio, aes(x = X.coordinate,y = Y.coordinate)) + geom_point()


# CAPTURE HISTORIES AURORA -------------------------------------------------------

traps <- trapinfo_Aurora

all_jaguars <- unique(captureHistory_Jaguar_Aurora$Individuals)
CH_jaguars <- matrix(0, nrow = length(all_jaguars), ncol = nrow(traps))
rownames(CH_jaguars) <- all_jaguars
colnames(CH_jaguars) <- traps$LOC_ID
for (i in 1:nrow(captureHistory_Jaguar_Aurora)) {
  index_trap <- which(colnames(CH_jaguars) == captureHistory_Jaguar_Aurora$Camera.trap.stations[i])
  index_individual <- which(rownames(CH_jaguars) == captureHistory_Jaguar_Aurora$Individuals[i])
  CH_jaguars[index_individual, index_trap] <- CH_jaguars[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}

jaguarsactive <- apply(CH_jaguars, 2, sum)

all_ocelots <- unique(captureHistory_Ocelot_Aurora$ANIMAL_ID)
CH_ocelots <- matrix(0, nrow = length(all_ocelots), ncol = nrow(traps))
rownames(CH_ocelots) <- all_ocelots
colnames(CH_ocelots) <- traps$LOC_ID
for (i in 1:nrow(captureHistory_Ocelot_Aurora)) {
  index_trap <- which(colnames(CH_ocelots) == captureHistory_Ocelot_Aurora$LOC_ID[i])
  index_individual <- which(rownames(CH_ocelots) == captureHistory_Ocelot_Aurora$ANIMAL_ID[i])
  CH_ocelots[index_individual, index_trap] <- CH_ocelots[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}

ocelotsactive <- apply(CH_ocelots, 2, sum)

ggplot(data = trapinfo_Aurora, aes(x = X, y = Y, color = jaguarsactive > 0)) + geom_point()
ggplot(data = trapinfo_Aurora, aes(x = X, y = Y, color = ocelotsactive > 0)) + geom_point()

S <- max(captureHistory_Ocelot_Aurora$SO, captureHistory_Jaguar_Aurora$Sampling.Occasions)

save(CH_jaguars, CH_ocelots, traps, S, file = "valeria_data_aurora.rda")

# CAPTURE HISTORIES MAG MEDIO -------------------------------------------------------

traps <- trapinfo_MagMedio

all_jaguars <- unique(captureHistory_Jaguar_MagMedio$Individuals)
CH_jaguars <- matrix(0, nrow = length(all_jaguars), ncol = nrow(traps))
rownames(CH_jaguars) <- all_jaguars
colnames(CH_jaguars) <- traps$Camera.trap.stations
for (i in 1:nrow(captureHistory_Jaguar_MagMedio)) {
  index_trap <- which(colnames(CH_jaguars) == captureHistory_Jaguar_MagMedio$Camera.trap.stations[i])
  index_individual <- which(rownames(CH_jaguars) == captureHistory_Jaguar_MagMedio$Individuals[i])
  CH_jaguars[index_individual, index_trap] <- CH_jaguars[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}

jaguarsactive <- apply(CH_jaguars, 2, sum)

all_ocelots <- unique(captureHistory_Ocelot_MagMedio$ANIMAL_ID)
CH_ocelots <- matrix(0, nrow = length(all_ocelots), ncol = nrow(traps))
rownames(CH_ocelots) <- all_ocelots
colnames(CH_ocelots) <- traps$Camera.trap.stations
for (i in 1:nrow(captureHistory_Ocelot_MagMedio)) {
  index_trap <- which(colnames(CH_ocelots) == captureHistory_Ocelot_MagMedio$LOC_ID[i])
  index_individual <- which(rownames(CH_ocelots) == captureHistory_Ocelot_MagMedio$ANIMAL_ID[i])
  CH_ocelots[index_individual, index_trap] <- CH_ocelots[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}

ocelotsactive <- apply(CH_ocelots, 2, sum)

ggplot(data = trapinfo_MagMedio, aes(x = X.coordinate, y = Y.coordinate, color = jaguarsactive > 0)) + geom_point()
ggplot(data = trapinfo_MagMedio, aes(x = X.coordinate, y = Y.coordinate, color = ocelotsactive > 0)) + geom_point()

S <- max(captureHistory_Ocelot_MagMedio$SO, captureHistory_Jaguar_MagMedio$Sampling.Occasions)

save(CH_jaguars, CH_ocelots, traps, S, file = "valeria_data_magmedio.rda")
