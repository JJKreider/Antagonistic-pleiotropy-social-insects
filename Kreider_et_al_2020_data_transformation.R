# Publication: Antagonistic pleiotropy and the evolution of ageing in social insects
# Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
# data transformation R-script

# this script summarizes the results from several replicate simulations of the same simulation scenario
library(tidyr)
replicates <- c(1:20)

# read data
setwd("") # set working directory to folder with simulation output from the same simulation scenario with extrinsic mortality = 0.0
fileNames00 <- Sys.glob("data_*.csv")
fileList00 <- lapply(fileNames00, function(x) read.csv(x, header = TRUE))

setwd("") # set working directory to folder with simulation output from the same simulation scenario with extrinsic mortality = 0.2
fileNames02 <- Sys.glob("data_*.csv")
fileList02 <- lapply(fileNames02, function(x) read.csv(x, header = TRUE))

lifespans <- data.frame()
lifespanRatios <- data.frame()

for (i in replicates) {
  lifespans[i, 1] <- fileList00[[i]]$avg.queen.lifespan[nrow(fileList00[[i]])]
  lifespans[i + length(replicates), 1] <- fileList00[[i]]$avg.worker.lifespan[nrow(fileList00[[i]])]
  lifespans[i, 2] <- "0.0"
  lifespans[i + length(replicates), 2] <- "0.0"
  lifespans[i, 3] <- "queen"
  lifespans[i + length(replicates), 3] <- "worker"

  lifespans[i + 2 * length(replicates), 1] <- fileList02[[i]]$avg.queen.lifespan[nrow(fileList02[[i]])]
  lifespans[i + 3 * length(replicates), 1] <- fileList02[[i]]$avg.worker.lifespan[nrow(fileList02[[i]])]
  lifespans[i + 2 * length(replicates), 2] <- "0.2"
  lifespans[i + 3 * length(replicates), 2] <- "0.2"
  lifespans[i + 2 * length(replicates), 3] <- "queen"
  lifespans[i + 3 * length(replicates), 3] <- "worker"
  
  lifespanRatios[i, 1] <- fileList00[[i]]$avg.queen.lifespan[nrow(fileList00[[i]])] / fileList00[[i]]$avg.worker.lifespan[nrow(fileList00[[i]])]
  lifespanRatios[i + length(replicates), 1] <- fileList02[[i]]$avg.queen.lifespan[nrow(fileList02[[i]])] / fileList02[[i]]$avg.worker.lifespan[nrow(fileList02[[i]])]
  lifespanRatios[i, 2] <- "0.0"
  lifespanRatios[i + length(replicates), 2] <- "0.2"
}

colnames(lifespans) <- c("lifespan", "ext", "caste")
colnames(lifespanRatios) <- c("lifespanRatio", "ext")

intrinsicSurvival00 <- data.frame()
fecundity00 <- data.frame()

for (i in replicates) {
  for (j in c(1:40)) {
    intrinsicSurvival00[i, j] <- fileList00[[i]][nrow(fileList00[[i]]), 7 + j]
    fecundity00[i, j] <- fileList00[[i]][nrow(fileList00[[i]]), 47 + j]
  }
}

intrinsicSurvival02 <- data.frame()
fecundity02 <- data.frame()

for (i in replicates) {
  for (j in c(1:40)) {
    intrinsicSurvival02[i, j] <- fileList02[[i]][nrow(fileList02[[i]]), 7 + j]
    fecundity02[i, j] <- fileList02[[i]][nrow(fileList02[[i]]), 47 + j]
  }
}

fecundity00 <- fecundity00[ ,1:20]
fecundity02 <- fecundity02[ ,1:20]
fecundity00 <- t(fecundity00)
fecundity02 <- t(fecundity02)

fecundity00 <- cbind(fecundity00, rep("0.0", 20))
fecundity02 <- cbind(fecundity02, rep("0.2", 20))

fecundityData <- data.frame()
fecundity00 <- as.data.frame(fecundity00)
fecundity02 <- as.data.frame(fecundity02)
colnames(fecundity00) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext")
colnames(fecundity02) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext")
fecundity00 <- gather(fecundity00, replicate, fecundity, -ext) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(fecundity))
fecundity02 <- gather(fecundity02, replicate, fecundity, -ext) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(fecundity))
fecundity00$fecundity <- as.numeric(fecundity00$fecundity)
fecundity02$fecundity <- as.numeric(fecundity02$fecundity)
fecundityData <- rbind(fecundity00, fecundity02)

lxQu00 <- data.frame()
lxQu02 <- data.frame()
lxW00 <- data.frame()
lxW02 <- data.frame()

intrinsicSurvival00 <- as.data.frame(t(intrinsicSurvival00))
intrinsicSurvival02 <- as.data.frame(t(intrinsicSurvival02))

lxQu00 <- intrinsicSurvival00[1, 1:20]
lxQu02 <- intrinsicSurvival02[1, 1:20]
lxW00 <- intrinsicSurvival00[21, 1:20]
lxW02 <- intrinsicSurvival02[21, 1:20]

for (i in replicates) {
  for (j in c(2:20)) {
    lxQu00[j, i] <- lxQu00[j - 1, i] * intrinsicSurvival00[j, i]
    lxQu02[j, i] <- lxQu02[j - 1, i] * intrinsicSurvival02[j, i]
    lxW00[j, i] <- lxW00[j - 1, i] * intrinsicSurvival00[j + 20, i]
    lxW02[j, i] <- lxW02[j - 1, i] * intrinsicSurvival02[j + 20, i]
  }
}

survivalData <- data.frame()
lxQu00 <- as.data.frame(lxQu00)
lxQu02 <- as.data.frame(lxQu02)
lxW00 <- as.data.frame(lxW00)
lxW02 <- as.data.frame(lxW02)

lxQu00 <- cbind(lxQu00, rep("0.0", 20), rep("queen", 20))
lxQu02 <- cbind(lxQu02, rep("0.2", 20), rep("queen", 20))
lxW00 <- cbind(lxW00, rep("0.0", 20), rep("worker", 20))
lxW02 <- cbind(lxW02, rep("0.2", 20), rep("worker", 20))

colnames(lxQu00) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext", "caste")
colnames(lxQu02) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext", "caste")
colnames(lxW00) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext", "caste")
colnames(lxW02) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "ext", "caste")

lxQu00 <- gather(lxQu00, replicate, survival, -ext, -caste) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(survival))
lxQu02 <- gather(lxQu02, replicate, survival, -ext, -caste) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(survival))
lxW00 <- gather(lxW00, replicate, survival, -ext, -caste) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(survival))
lxW02 <- gather(lxW02, replicate, survival, -ext, -caste) %>% dplyr::group_by(replicate) %>% dplyr::mutate(ageclass = 1:length(survival))

survivalData <- rbind(lxQu00, lxQu02, lxW00, lxW02)

setwd("") # set output directory
write.csv(lifespans, "lifespans.csv")
write.csv(lifespanRatios, "lifespanRatios.csv")
write.csv(survivalData, "survivalData.csv")
write.csv(fecundityData, "fecundityData.csv")

