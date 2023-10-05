################################################################################################
#
# Code to fit count spatial model using Bisquare Basis Functions
#
################################################################################################
################################################################################################
rm(list=ls())

library(fields) ; library(mvtnorm) ; library(classInt); 
library(rstan)
options(mc.cores = parallel::detectCores())
setwd("~/Dropbox/GMU_BenLee/HandbookMCMC/Precip/run_bisquare_small/") # Set Directory
source("sharedFunctions.R") # Helper file with functions
load("../samples/precipDataSmall.RData") # Dataset
load("../samples/BisquareBasisFunctions.RData") # Dataset
################################################################################################
# Preliminaries for Stan
################################################################################################
# Inputs for Stan
k=ncol(XMat) # Design Matrix dimensions 
# basis Functions
distA<-rdist(gridLocationMod,knotsA) 
distB<-rdist(gridLocationMod,knotsB)
distC<-rdist(gridLocationMod,knotsC)
basisA<-bisquareFunc(dist = distA , eps = epsA) # Basis Functions - Resolution A
basisB<-bisquareFunc(dist = distB , eps = epsB) # Basis Functions - Resolution B
basisC<-bisquareFunc(dist = distC , eps = epsC) # Basis Functions - Resolution C
rm(distA,distB,distC)
M<-cbind(basisA,basisB,basisC)
rm(basisA,basisB,basisC)


distA<-rdist(gridLocationCV,knotsA) 
distB<-rdist(gridLocationCV,knotsB)
distC<-rdist(gridLocationCV,knotsC)
basisA<-bisquareFunc(dist = distA , eps = epsA) # Basis Functions - Resolution A
basisB<-bisquareFunc(dist = distB , eps = epsB) # Basis Functions - Resolution B
basisC<-bisquareFunc(dist = distC , eps = epsC) # Basis Functions - Resolution C
rm(distA,distB,distC)
MCV<-cbind(basisA,basisB,basisC)
rm(basisA,basisB,basisC)

p<-ncol(M) # Number of basis functions
n<-length(obsMod)
################################################################################################
# Stan Setup
################################################################################################
# Fit Model Using Stan
iter=2000
gamma_samples<-stan(file="gammaBasis.stan",
             data = list(N=n,
                         K=k,
                         P=p,
                         y=obsMod,
                         X=XMat,
                         M=M),
             pars = c("betas","deltas","tau","inverse_phi"),
             iter=iter,chains = 1)

#Extract Pertinent Information from output
totTime<-sum(get_elapsed_time(gamma_samples)) # Model Fitting Time
parMat<-cbind(rstan::extract(gamma_samples, c('betas'))[[1]], # Extract samples for Model Parameters
              rstan::extract(gamma_samples, c('tau'))[[1]],
              rstan::extract(gamma_samples, c('inverse_phi'))[[1]])
deltaMat<-rstan::extract(gamma_samples, c('deltas'))[[1]] # Extract samples for basis coefficients

print(gamma_samples)
stan_trace(gamma_samples,pars = c('betas','tau','inverse_phi'))
################################################################################################
# Summarize Samples
################################################################################################
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(parMat,
                                       totTime=totTime),3)
summaryMat[[2]]<-round(summaryFunction(deltaMat,
                                       totTime=totTime),3)
summaryMat[[1]] # Summary Table
apply(summaryMat[[2]],1,mean) # Mean results for the basis coefficients

################################################################################################
# Trace Plots
################################################################################################
par(mfrow=c(4,2),mar=c(2,2,2,2))  
sampInd<-floor(seq(1,nrow(parMat),length.out = 1000))
for(i in 1:3){
  plot.ts(parMat[sampInd,i], 
          main=paste("beta",i,sep="")); 
  abline(h=beta[i],col="red",lwd=2)
  plot(density(parMat[sampInd,i]),
       main=paste("beta",i,sep="")); 
  abline(v=beta[i],col="red",lwd=2)
}

plot.ts(parMat[sampInd,5], 
        main="alpha"); 
abline(h=alpha_gamma,col="red",lwd=2)
plot(density(parMat[sampInd,5]),
     main="alpha"); 
abline(v=alpha_gamma,col="red",lwd=2)
################################################################################################
# Compute Mean Squared Prediction Error
################################################################################################
mcmcDat<-list()
mcmcDat[[1]]<-parMat; colnames(mcmcDat[[1]])<-c("beta1","beta2","beta3","tau","alpha")
mcmcDat[[2]]<-deltaMat
cvSummary<-cvFunction.gamma(mcmcDat=mcmcDat,mBase=MCV,
                           XMatCV=XMatCV,obsCV=obsCVGamma)
print(cvSummary[[1]]) #MSPE

par(mfrow=c(1,2))
plotRF(dat = cvSummary[[3]], location = CVgridLocation)
plotRF(dat = pWCVGamma, location = CVgridLocation)



# basis Functions
distA<-rdist(figLocation,knotsA) 
distB<-rdist(figLocation,knotsB)
distC<-rdist(figLocation,knotsC)
basisA<-bisquareFunc(dist = distA , eps = epsA) # Basis Functions - Resolution A
basisB<-bisquareFunc(dist = distB , eps = epsB) # Basis Functions - Resolution B
basisC<-bisquareFunc(dist = distC , eps = epsC) # Basis Functions - Resolution C
rm(distA,distB,distC)
M_Fig<-cbind(basisA,basisB,basisC)
rm(basisA,basisB,basisC)

dim(M_Fig)
dim(deltaMat)
predW_Fig<-apply(M_Fig%*%t(deltaMat),1,mean)
logMean_Fig<-apply(XMatFig%*%t(parMat[,1:3])+M_Fig%*%t(deltaMat),1,mean)


library(raster) ; library(rgdal) ; library(rmutil); library(fields)
library(RColorBrewer); library(colorspace)



par(mfrow=c(2,2))
x <- raster()
x <- raster(ncol=sqrt(nFig), nrow=sqrt(nFig), xmn=0, xmx=1, ymn=0, ymx=1)
values(x) <- as.numeric(gpWFig)
plot(x, main="True Random Effects",cex.main=2,
     zlim=c(-2.82,2.432),
     col=tim.colors(1000),axes=FALSE, box=FALSE)

x <- raster()
x <- raster(ncol=sqrt(nFig), nrow=sqrt(nFig), xmn=0, xmx=1, ymn=0, ymx=1)
values(x) <- as.numeric(pred_Fig)
plot(x, main="Prediction",cex.main=2,
     zlim=c(-2.82,2.432),
     col=tim.colors(1000),axes=FALSE, box=FALSE)

range(c(log(pWFigGamma), logMean_Fig))
x <- raster()
x <- raster(ncol=sqrt(nFig), nrow=sqrt(nFig), xmn=0, xmx=1, ymn=0, ymx=1)
values(x) <- as.numeric(log(pWFigGamma))
plot(x, main="True Log Mean",cex.main=2,
     zlim=c(-4.17,3.45),
     col=tim.colors(1000),axes=FALSE, box=FALSE)

x <- raster()
x <- raster(ncol=sqrt(nFig), nrow=sqrt(nFig), xmn=0, xmx=1, ymn=0, ymx=1)
values(x) <- as.numeric(logMean_Fig)
plot(x, main="Prediction: Log Mean",cex.main=2,
     zlim=c(-4.17,3.45),
     col=tim.colors(1000),axes=FALSE, box=FALSE)



# Table



################################################################################################
# Save Data
################################################################################################
save(summaryMat,parMat,deltaMat,totTime,cvSummary,file="../output/mcmc_Stan_small.RData") # Save Data
