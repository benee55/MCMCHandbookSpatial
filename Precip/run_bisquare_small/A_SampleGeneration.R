########################################################
# A_SampleGeneration_Cluster.R 
# Partitions spatial domain into clusters
########################################################
rm(list=ls())
library(fields);library(classInt)
setwd("~/Dropbox/GMU_BenLee/HandbookMCMC/Precip/run_bisquare_small/")
load("../samples/precipData.RData")
set.seed(12345)

#IDs for smaller samples
n<-20000
nCV<-4000
idModSmall<-sample(idMod, size =n , replace = FALSE )
idCVSmall<-sample(idMod, size =nCV , replace = FALSE )


###### Data #####
# Generate Split into Training and Testing Samples
gridLocationMod<-comboLocation[idModSmall,]
gridLocationCV<-comboLocation[idCVSmall,]
XMatFull<-cbind(1,comboLocation[,1])
XMat<-XMatFull[idModSmall,]
XMatCV<-XMatFull[idCVSmall,]

#Observations
obsMod<-obsFull[idModSmall]
obsCV<-obsFull[idCVSmall]

rm(comboLocation, XMatFull,obsFull, idCV, idMod)
# Save Data
save(gridLocationCV, gridLocationMod, 
     XMat, XMatCV, 
     idModSmall, idCVSmall, 
     obsCV, obsMod,
     file="../samples/precipDataSmall.RData")





