rm(list=ls())
library(fields);library(classInt);library(nimble)# Load packages
library(nimble);
####################################
source("source/sharedFunctions.R") # File with useful functions
load("../samples/precipData.RData")
load("../samples//BisquareBasisFunctions.RData")
source("source/batchmeans.R")# batchmeans source

n = nrow(gridLocationMod)
n.pred = nrow(gridLocationCV)

CVgridLocation<-gridLocationCV
gridLocation<-gridLocationMod
rm(gridLocationCV,gridLocationMod,comboLocation, XMatFull, obsFull)
##################
# Generate Basis Functions
##################
distA<-rdist(gridLocation,knotsA) 
distB<-rdist(gridLocation,knotsB)
distC<-rdist(gridLocation,knotsC)
basisA<-bisquareFunc(dist = distA , eps = epsA) # Basis Functions - Resolution A
basisB<-bisquareFunc(dist = distB , eps = epsB) # Basis Functions - Resolution B
basisC<-bisquareFunc(dist = distC , eps = epsC) # Basis Functions - Resolution C
rm(distA,distB,distC)
M<-cbind(basisA,basisB,basisC)
rm(basisA,basisB,basisC)
distA_CV<-rdist(CVgridLocation,knotsA) 
distB_CV<-rdist(CVgridLocation,knotsB)
distC_CV<-rdist(CVgridLocation,knotsC)
basisA_CV<-bisquareFunc(dist = distA_CV , eps = epsA) # Basis Functions - Resolution A
basisB_CV<-bisquareFunc(dist = distB_CV , eps = epsB) # Basis Functions - Resolution B
basisC_CV<-bisquareFunc(dist = distC_CV , eps = epsC) # Basis Functions - Resolution C
rm(distA_CV,distB_CV,distC_CV)
M_CV<-cbind(basisA_CV,basisB_CV,basisC_CV)
rm(basisA_CV,basisB_CV,basisC_CV)
##################
n <- length(obsMod); k <- ncol(XMat); p<-ncol(M)
Q<-diag(p)
mn=rep(0,max(k,p))
print(n);print(k);print(p)

keepInR <- list(X=XMat, M=M, pCov=100*diag(k))
consts   <- list(n=n, k=k , p=p)
data     <- list(Z=obsMod, Q=Q , mn=mn)
inits    <- list(beta=c(0,0),tau=2 , 
                 delta=rep(0,p),  shape=1)
source("source/nimble.R")
deregisterDistributions("dGammaVector") 
deregisterDistributions("dRdmnormB")

Rmodel<-nimbleModel(code=model_gamma_string, data = data,  
                    constants=consts, inits = inits,calculate = FALSE)
cRmodel <- compileNimble(Rmodel)
confRModel <- configureMCMC(Rmodel, print = TRUE, useConjugacy = FALSE)
confRModel$printSamplers()
confRModel$addMonitors(c("beta","tau","delta","shape"))  
mcmcRmodel <- buildMCMC(confRModel)
cmcmcRmodel <- compileNimble(mcmcRmodel)


niter=2000
pt<-proc.time()
cmcmcRmodel$run(niter)
ptFinal<-proc.time()-pt
ptFinal
chain_output <- as.matrix(cmcmcRmodel$mvSamples)
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(chain_output,
                                       totTime=ptFinal[3], bmseThresh = 0.01),3)

save(summaryMat,chain_output,ptFinal,file="../output/mcmcPrecip_bisquare.RData")

for(i in 1:50){
  pt<-proc.time()
  cmcmcRmodel$run(niter, reset = FALSE)
  ptFinal<-ptFinal+(proc.time()-pt)
  chain_output <- as.matrix(cmcmcRmodel$mvSamples)
  summaryMat<-list()
  summaryMat[[1]]<-round(summaryFunction(chain_output,
                                         totTime=ptFinal[3], bmseThresh = 0.01),3)
  save(summaryMat,chain_output,ptFinal,file="../output/mcmcPrecip_bisquare.RData")
}

