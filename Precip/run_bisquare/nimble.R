################################################################
#
# Nimble Functions for fitting hierarchical spatial models
#
################################################################
#
# Covariance Functions
#
################################################################

# for a vector of Gamma random variables
dGammaVector <- nimbleFunction(
  run = function(x    = double(1),
                 shape = double(0),
                 linPred = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    
    returnType(double(0))
    shape = shape
    rate = shape/exp(linPred)
    logLike <- sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
    if(log) return(logLike) else return(exp(logLike))
  }
  
)


## Use R's lexical scoping to create R functions to be called
## via nimbleRcall that will find large objects like
## X and pCov in their enclosing environment.
makeRfuns <- function(keepInR) {
  X <- keepInR$X
  pCov <- keepInR$pCov
  chol_pCov <- chol(pCov)
  M <- keepInR$M
  
  RcalculateXB <- function(beta) {
    as.numeric(X %*% beta)
  }
  RcalculateW <- function(delta) {
    as.numeric(M %*% delta)
  }
  dRdmnormB <- function(x, mean, log) {
    nimble::dmnorm_chol(x, mean, cholesky = chol_pCov, prec_param = FALSE, log = log)
  }
  list(RcalculateXB = RcalculateXB, RcalculateW=RcalculateW, dRdmnormB = dRdmnormB)
}

## Create the functions and put them in .GlobalEnv
Rfuns <- makeRfuns(keepInR)
RcalculateXB_internal <- Rfuns$RcalculateXB ## I think these need to be in .GlobalEnv
RcalculateW_internal <- Rfuns$RcalculateW ## I think these need to be in .GlobalEnv
dRdmnormB_internal <- Rfuns$dRdmnormB

## Define the nimbleRcall nimbleFunctions to call the R functions.
RcalculateXB <- nimbleRcall(function(beta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateXB_internal')
RcalculateW <- nimbleRcall(function(delta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateW_internal')
dRdmnormB <- nimbleRcall(function(x = double(1), mean = double(1), log = integer(0, default = 0)) {}, returnType = double(), Rfun = 'dRdmnormB_internal')

######################################################################
# Bisquare for Positive data - Model fit using nimble
######################################################################
model_gamma_string <- nimbleCode({
  # Data Model
  Z[1:n] ~ dGammaVector(shape=shape ,linPred=mu[1:n])
  mu[1:n] <- (XB[1:n]+W[1:n])
  XB[1:n] <- RcalculateXB(beta[1:k])
  W[1:n] <- RcalculateW(delta[1:p])
  
  precMat[1:p,1:p]<- tau * Q[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  
  # Parameter Model
  beta[1:k] ~ dRdmnormB(mean = mn[1:k])
  tau   ~   dgamma(shape=0.5,rate=2000) 
  shape   ~   dinvgamma(shape=0.5,rate=0.5) 
})

