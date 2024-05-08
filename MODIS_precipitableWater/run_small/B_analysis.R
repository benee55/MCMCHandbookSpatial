rm(list=ls())
library(fields);library(classInt); library(nimble)
####################################
load("../samples/precipData.RData")
load("../samples/BisquareBasisFunctions.RData")
source("../source/sharedFunctions.R") # File with useful functions
n = nrow(gridLocationMod)
n.pred = nrow(gridLocationCV)
CVgridLocation<-gridLocationCV
rm(gridLocationCV,gridLocationMod,comboLocation, XMatFull, obsFull)
distA_CV<-rdist(CVgridLocation,knotsA) 
distB_CV<-rdist(CVgridLocation,knotsB)
distC_CV<-rdist(CVgridLocation,knotsC)
basisA_CV<-bisquareFunc(dist = distA_CV , eps = epsA) # Basis Functions - Resolution A
basisB_CV<-bisquareFunc(dist = distB_CV , eps = epsB) # Basis Functions - Resolution B
basisC_CV<-bisquareFunc(dist = distC_CV , eps = epsC) # Basis Functions - Resolution C
rm(distA_CV,distB_CV,distC_CV)
M_CV<-cbind(basisA_CV,basisB_CV,basisC_CV)
rm(basisA_CV,basisB_CV,basisC_CV)


source("../source/batchmeans.r")# batchmeans source
load("../output/mcmcPrecip_small.RData")
betaMat<-chain_output[,1:2]
deltaMat<-chain_output[,3:86]
linearPred<-XMatCV%*%apply(betaMat,2,mean) + M_CV%*%apply(deltaMat,2,mean)

predCV<-expectedValue<-exp(linearPred)

# Prediction - Root Mean squared error
rmspe<-sqrt(mean((predCV-obsCV)^2))

# Summary: Seems like the global option performs better 
save(rmspe,predCV,file="../output/resultsMCMC_small.RData")
