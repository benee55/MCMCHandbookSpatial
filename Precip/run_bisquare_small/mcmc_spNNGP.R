################################################################################################
################################################################################################
rm(list=ls())

library(fields) ; library(mvtnorm) ; library(classInt); 
setwd("~/Dropbox/GMU_BenLee/HandbookMCMC/Precip/run_bisquare_small/") # Set Directory
load("../samples/precipDataSmall.RData") # Dataset
load("../samples/BisquareBasisFunctions.RData") # Dataset

library(spNNGP)

n.samples <- 10000

starting <- list(phi = 1/1, sigma.sq = 1, tau.sq = 1)

priors <- list(phi.Unif = c(1 / 3, 1 / 0.001), sigma.sq.IG = c(2, 1), tau.sq.IG = c(2,
                                                                                  1))

cov.model <- "exponential"
obsMod<-log(obsMod)


tuning <- list(phi = 0.05)

# gridLocationMod<-gridLocationMod
# obsMod<-obsMod[1:100]
# XMat<-XMat[1:100,]
# gridLocationMod<-gridLocationMod[1:100,]
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
                # n.omp.threads = 1, 
                n.report = 500, 
                fit.rep = TRUE,
                sub.sample = list(start = 1000),
                return.neighbor.info = TRUE)




plot.ts(sim.s$p.beta.samples[,1])
plot.ts(sim.s$p.beta.samples[,2])
plot.ts(sim.s$p.theta.samples[,1])
plot.ts(sim.s$p.theta.samples[,2])
plot.ts(sim.s$p.theta.samples[,3])
sim.s$run.time
