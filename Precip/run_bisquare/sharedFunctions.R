
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

# Rank Selection
MLE_FindRank_Poisson_PICAR<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,MoransOperatorEig){
  CVMSPE<-vector("numeric")
  betaMatList<-list()
  betaMatList[[1]]<-betaMatList[[2]]<-matrix(NA,nrow=length(dimSeq),ncol=3)
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-c(1:dimSeq[jk])
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM])
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM])
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "poisson")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    betaMatList[[1]][jk,]<-rbind(estMean,lowCI,highCI)[,1]
    betaMatList[[2]][jk,]<-rbind(estMean,lowCI,highCI)[,2]
    foo<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs)  
    predCV<-foo
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(list(CVMSPE,betaMatList))
}


MLE_FindRank_Poisson_basis<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,modBasis, cvBasis){
  CVMSPE<-vector("numeric")
  betaMatList<-list()
  betaMatList[[1]]<-betaMatList[[2]]<-matrix(NA,nrow=length(dimSeq),ncol=3)
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-c(1:dimSeq[jk])
    mBase<-(modBasis[,keepM])
    mBaseCV<-(cvBasis[,keepM])
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "poisson")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    betaMatList[[1]][jk,]<-rbind(estMean,lowCI,highCI)[,1]
    betaMatList[[2]][jk,]<-rbind(estMean,lowCI,highCI)[,2]
    foo<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs)  
    predCV<-foo
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(list(CVMSPE,betaMatList))
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


# Cross validation

# Rank Selection
MLE_FindRank_Binary_PICAR<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,MoransOperatorEig){
  CVMSPE<-vector("numeric")
  betaMatList<-list()
  betaMatList[[1]]<-betaMatList[[2]]<-matrix(NA,nrow=length(dimSeq),ncol=3)
  for(jk in 1:length(dimSeq)){
    if(jk%%10==0){print(jk)}
    print(dimSeq[jk])
    keepM<-c(1:dimSeq[jk])
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM])
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM])
    lm1<-glm(obsMod~0+cbind(XMat,as.matrix(mBase)),family = "binomial")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[c(length(lm1$coefficients)-1,length(lm1$coefficients))]
    betaMatList[[1]][jk,]<-rbind(estMean,lowCI,highCI)[,1]
    betaMatList[[2]][jk,]<-rbind(estMean,lowCI,highCI)[,2]
    predCV<-ifelse(logitTr(cbind(XMatCV,as.matrix(mBaseCV))%*%coeffs)>0.5,1,0)
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(CVMSPE)
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
