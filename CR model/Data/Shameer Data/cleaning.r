setwd("C:/Users/alexd/Dropbox/R Folder/PhD/Spatial/CR Model/Data/Shameer Data")

traps <- read.csv(file = "trap.csv", stringsAsFactors = FALSE)

workingmatrix <- read.table(file = "working_matrix.txt", header = T, stringsAsFactors = FALSE)
nonWorkingTraps <- workingmatrix$P1[which(workingmatrix$X1 == 0)]

traps <- traps[!(traps$CODE %in% nonWorkingTraps),]

library(ggplot2)

ggplot(data = traps, aes(x = X, y = Y)) + geom_point()

capture_leopard <- read.table("capture_leopard.txt", header = F, stringsAsFactors = F)
colnames(capture_leopard) <- c("Session", "ID", "Occasion", "Trap")

sapply(capture_leopard$Trap, function(x){
  x %in% traps$CODE
})

capture_tiger <- read.table("capture_tiger.txt", header = F, stringsAsFactors = F)
colnames(capture_tiger) <- c("Session", "ID", "Occasion", "Trap")

sapply(capture_tiger$Trap, function(x){
  x %in% traps$CODE
})

# create  CR data

all_leopards <- unique(capture_leopard$ID)
CH_leopards <- matrix(0, nrow = length(all_leopards), ncol = nrow(traps))
rownames(CH_leopards) <- all_leopards
colnames(CH_leopards) <- traps$CODE
for (i in 1:nrow(capture_leopard)) {
  index_trap <- which(colnames(CH_leopards) == capture_leopard$Trap[i])
  index_individual <- which(rownames(CH_leopards) == capture_leopard$ID[i])
  CH_leopards[index_individual, index_trap] <- 1
  print(paste0(index_individual," - ",index_trap))
}

all_tigers <- unique(capture_tiger$ID)
CH_tigers <- matrix(0, nrow = length(all_tigers), ncol = nrow(traps))
rownames(CH_tigers) <- all_tigers
colnames(CH_tigers) <- traps$CODE
for (i in 1:nrow(capture_tiger)) {
  index_trap <- which(colnames(CH_tigers) == capture_tiger$Trap[i])
  index_individual <- which(rownames(CH_tigers) == capture_tiger$ID[i])
  CH_tigers[index_individual, index_trap] <- 1
  print(paste0(index_individual," - ",index_trap))
}

activeTrapsLeopard <- colnames(CH_leopards)[apply(CH_leopards,2,sum) > 0]
activeTrapsTiger <- colnames(CH_tigers)[apply(CH_tigers,2,sum) > 0]

traps$ACTIVE_LEOPARD <- F
traps$ACTIVE_LEOPARD[traps$CODE %in% activeTrapsLeopard] <- T

traps$ACTIVE_TIGER <- F
traps$ACTIVE_TIGER[traps$CODE %in% activeTrapsTiger] <- T

ggplot(data = traps, aes(x = X, y = Y, color = ACTIVE_LEOPARD)) + geom_point() + 
  theme(legend.title = element_blank()) + scale_color_discrete()
ggplot(data = traps, aes(x = X, y = Y)) + geom_point(aes(color = as.factor(ACTIVE_TIGER))) + 
  theme(legend.title = element_blank()) 

save(traps, CH_leopards, CH_tigers, file = "data.rda")


# PLOT MAPS ---------------------------------------------------------------

require(maptools); library(dplyr); library(lubridate); library(ggplot2); library(gridExtra)
mapsINDIA <- readRDS("gadm36_IND_3_sf.rds")
india = fortify(mapsINDIA)

ggplot() + geom_polygon(data=india, aes(x=long, y=lat,group=group),color="black",alpha=0.6) + 
  coord_cartesian(xlim = c(76.5, 78), ylim = c(8.5, 10)) + 
  geom_point(data = traps, aes(x = X, y = Y, color = "red"))
