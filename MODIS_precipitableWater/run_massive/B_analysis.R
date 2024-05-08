rm(list=ls())
library(fields);library(nimble)
setwd("~/MCMCHandbook/PrecipitationExample/output/")
####################################
load("../samples/precipData.RData")
load("../samples/BisquareBasisFunctions.RData")
source("../source//sharedFunctions.R") # File with useful functions
n = nrow(gridLocationMod)
n.pred = nrow(gridLocationCV)
CVgridLocation<-gridLocationCV
rm(gridLocationCV,gridLocationMod)
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
load("mcmcPrecip_bisquare.RData")

useInterval<-round(seq(20000,nrow(chain_output),length.out=20))
meanVal<-mcseVal<-matrix(NA, nrow=length(useInterval), ncol=2)
rmspeVect<-vector("numeric")

for(k in 1:length(useInterval)){
  print(k)
  foo<-bmmat(chain_output[1:useInterval[k],1:2])
  meanVal[k,]<-foo[,1] ; mcseVal[k,]<-foo[,2]

  
  # RMSPE
  betaMat_nimble<-chain_output[1:useInterval[k],1:2]
  deltaMat_nimble<-chain_output[1:useInterval[k],3:86]

  linearPred<-XMatCV%*%apply(betaMat_nimble,2,mean) + M_CV%*%apply(deltaMat_nimble,2,mean)
  # linearPred<-XMatCV%*%t(betaMat_nimble) + M_CV%*%apply(deltaMat_nimble,2,mean)
  predCV<-exp(linearPred)
  # Prediction - Root Mean squared error
  rmspeVect[k]<-sqrt(mean((predCV-obsCV)^2))

}

par(mfrow=c(1,1))
keepInd<-seq(0,length(useInterval),by=5)
keepUseInterval<-seq(0,nrow(chain_output)/1000,length.out=5)
plot(x=1:length(rmspeVect),y=rmspeVect,typ="l", 
     cex.lab=1.4,lwd=2,
     ylim=range(rmspeVect-0.002,rmspeVect+0.002), 
     ylab="RMSPE", xlab="iterations for nimble (1000's)", xaxt="n")
ticks = keepInd
axis(side = 1, at = ticks , 
     labels = keepUseInterval)




plot(x=1:length(mcseVal[,1]),y=mcseVal[,1],typ="l", ylim=range(mcseVal,mcseVal_stan))
lines(x=1:length(mcseVal_stan[,1]),y=mcseVal_stan[,1], col="red")

plot(x=1:length(mcseVal[,2]),y=mcseVal[,2],typ="l", ylim=range(mcseVal,mcseVal_stan))
lines(x=1:length(mcseVal_stan[,2]),y=mcseVal_stan[,2], col="red")

plot.ts(chain_output[,1])


# Summary: Seems like the global option performs better 
# save(rmspe,predCV,file="../output/resultsMCMC_small.RData")
