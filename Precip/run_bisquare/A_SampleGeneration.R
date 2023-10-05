########################################################
# A_SampleGeneration_Cluster.R 
# Partitions spatial domain into clusters
########################################################
rm(list=ls())
library(fields);library(classInt)
load("../samples/water_vapor_20190328.RData")
set.seed(12345)

# Process Data 
coords.all<-locs # Coordinates 
range.coords<-apply(coords.all,2,range) # Rescale data to be on the unit grid
coords.all<-coords.all+matrix(-(range.coords[1,]),nrow=nrow(coords.all),ncol=2,byrow = TRUE)
range.coords<-apply(coords.all,2,range);range.coords
coords.all<-coords.all*matrix(1/abs(range.coords[2,]),nrow=nrow(coords.all),ncol=2,byrow = TRUE)
apply(coords.all,2,range)

# Split sample into training and testing datasets
nTot<-length(z)
nCV=floor(0.1*nTot)
nMod=nTot-nCV;
idMod<-sample(1:nTot, nMod)
idCV<-(1:nTot)[-idMod]

# Generate Split into Training and Testing Samples
comboLocation<-as.matrix(coords.all)
rm(coords.all)
gridLocationMod<-comboLocation[idMod,]
gridLocationCV<-comboLocation[idCV,]
XMatFull<-cbind(1,comboLocation[,1])
XMat<-XMatFull[idMod,]
XMatCV<-XMatFull[idCV,]

#Observations
obsFull<-z # Observations
obsMod<-z[idMod]
obsCV<-z[idCV]

# Save Data
save(comboLocation,gridLocationCV, gridLocationMod, 
     XMatFull, XMat, XMatCV, 
     idMod, idCV, 
     obsFull, obsCV, obsMod,
     file="../samples/precipData.RData")

