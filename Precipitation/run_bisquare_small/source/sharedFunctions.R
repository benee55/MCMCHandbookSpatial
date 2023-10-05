
################################
plotRF<-function(dat,rangeDat=dat,label="Plot",location,length.out=10,pch=16,cex=1){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1,alpha = 1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=pch,cex=cex,
       main=label)
}

## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}

# Logit
logitTr<-function(x){
  foo<-exp(x)
  return(foo/(1+foo))
}
# Matern Cov Function + Acceptance Rate function
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}
# Calculate Acceptance Rate for MCMC algorithms
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}

################################################################
# Gamma Data
################################################################
cvFunction.gamma<-function(mcmcDat, # Matrix of posterior samples 
                           mBase,   # Moran's Basis Functions
                           XMatCV,  # Design Matrix (test)
                           obsCV,   # Test Data
                           burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
{
  basisMat<-mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # Index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  foo<-apply(cvPred+cvPreXB,1,mean) # Linear Predictor 
  predVal<-exp(foo) # Expected Value 
  cvSummary<-mean((predVal-obsCV)^2)  # Mean squared prediction error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat,predVal)) # MSPE and predicted values
}


# Summarize MCMC output
summaryFunction<-function(mcmcDat,totTime,bmseThresh=0.01){
  
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                    apply(mcmcDat,2,hpd),
                    apply(mcmcDat,2,accRateFunc),
                    bmmat(mcmcDat)[,2],
                    abs(apply(mcmcDat,2,mean))*bmseThresh,
                    apply(mcmcDat,2,ess),
                    apply(mcmcDat,2,ess)/totTime)
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                          "Accept","BMSE",paste(bmseThresh,"x mean"),
                          "ESS","ESS/sec")
  return(summaryMat)
}

# Covariance functions

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })


matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    
    return(result)
  })
