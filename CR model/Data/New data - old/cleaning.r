setwd("C:/Users/Alex/Dropbox/R Folder/PhD/Spatial/CR Model/Data/New data")

traps <- read.csv(file = "traps.csv")
tiger <- read.csv(file = "tiger.csv")
leopards <- read.csv(file = "leopards.csv")

tiger <- tiger[!duplicated(tiger),]
leopards <- leopards[!duplicated(leopards),]

library(ggplot2)


ggplot(data = NULL, aes(x = traps$x, y = traps$y)) + geom_point()


all_leopards <- unique(leopards$leopard_id)
CH_leopards <- matrix(0, nrow = length(all_leopards), ncol = nrow(traps))
rownames(CH_leopards) <- all_leopards
colnames(CH_leopards) <- traps$trapid
for (i in 1:nrow(leopards)) {
  index_trap <- which(colnames(CH_leopards) == leopards$trapid[i])
  index_individual <- which(rownames(CH_leopards) == leopards$leopard_id[i])
  CH_leopards[index_individual, index_trap] <- CH_leopards[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}


all_tigers <- unique(tiger$tigerid)
CH_tigers <- matrix(0, nrow = length(all_tigers), ncol = nrow(traps))
rownames(CH_tigers) <- all_tigers
colnames(CH_tigers) <- traps$trapid
for (i in 1:nrow(tiger)) {
  index_trap <- which(colnames(CH_tigers) == tiger$trapid[i])
  index_individual <- which(rownames(CH_tigers) == tiger$tigerid[i])
  CH_tigers[index_individual, index_trap] <- CH_tigers[index_individual, index_trap] + 1
  print(paste0(index_individual," - ",index_trap))
}

S_usage <- apply(traps[,-(1:3)], 1, sum)

table(apply(CH_leopards, 1, sum))
table(apply(CH_tigers, 1, sum))

# save(CH_leopards, CH_tigers, traps, S_usage, file = "yad_data_new.rda")

