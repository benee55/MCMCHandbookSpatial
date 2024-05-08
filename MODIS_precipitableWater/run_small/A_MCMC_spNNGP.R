################################################################################################
################################################################################################
rm(list=ls())
setwd("~/MCMCHandbook/PrecipitationExample/run_bisquare_small/")
library(fields) ; library(mvtnorm) ; library(classInt); 
load("../samples/precipDataSmall.RData") # Dataset

library(spNNGP)
n.samples <- 50000
starting <- list(phi = 1/1, sigma.sq = 1, tau.sq = 1)
priors <- list(phi.Unif = c(1 / 3, 1 / 0.001), sigma.sq.IG = c(2, 1), tau.sq.IG = c(2,
                                                                                  1))

cov.model <- "exponential"
obsMod<-log(obsMod)


tuning <- list(phi = 0.05)
ord <- order(gridLocationMod[, 1] + gridLocationMod[, 2])
sim.s <- spNNGP(formula = obsMod ~ 0+XMat, 
                coords = gridLocationMod, 
                starting = starting,
                tuning = tuning,
                priors = priors, 
                cov.model = cov.model, 
                n.samples = n.samples, 
                n.neighbors = 5,
                method = "latent",
                ord = ord, 
                n.report = 5000, 
                fit.rep = TRUE,
                sub.sample = list(start = 100),
                return.neighbor.info = TRUE)

totTime<-sim.s$run.time
 # save(sim.s , totTime , file = "../output/spNNGP_small.RData")
# load("../output/spNNGP_small.RData")
# totTime<-sim.s$run.time
# Prediction 
p.s <- predict(sim.s, X.0 = XMatCV, coords.0 = gridLocationCV, n.omp.threads=1)
predCV<-expectedValue<-exp(apply(p.s$p.y.0,1,mean,na.rm=TRUE))
#Root Mean squared error
rmspe<-sqrt(mean((predCV-obsCV)^2))
save(predCV, rmspe,totTime,file = "../output/spNNGP_small_2.RData")
# save(sim.s,p.s, rmspe,file = "../output/spNNGP_small.RData")
