
## glppmGroupedLassoFit <- function(model,initPar=NULL,fisherInformation=TRUE,control=list(),gamma=1,...) {
##   nrPar <- parDim <- dim(getModelMatrix(model))[2]
##   modelFixedPar <- model@fixedCoefficients
##   Omega <- model@Omega

##   ### Figure out which terms there are in the formula
  
## #  for(i in terms) {
## #    notTermsPar <- 
## #    fixedPar <- unique(c(notTermsPar$,modelFixedPar$))
    
    
##     mll <- function(par,...) {
##       tmpPar[-fixedPar$which] <- par
##       computeMinusLogLikelihood(model,tmpPar,...) + par %*% Omega[-fixedPar$which,-fixedPar$which] %*% par
##     }
    
##     dmll <- function(par,...) {
##       tmpPar[-fixedPar$which] <- par
##       computeDMinusLogLikelihood(model,tmpPar,...)[-fixedPar$which] + 2*par %*% Omega[-fixedPar$which,-fixedPar$which]
##     }
    
    
##       if(length(fixedPar) == 0) {
##         model@coefficients <- model@optimResult$par
##       } else {
##         model@coefficients[-fixedPar$which] <- model@optimResult$par
##         model@coefficients[fixedPar$which] <- fixedPar$value
##       }
##     names(model@coefficients) <- dimnames(model@modelMatrix)[[2]]
    
##   mll <- function(par,...) computeMinusLogLikelihood(model,par,...) + par %*% Omega %*% par
##   dmll <- function(par,...) computeDMinusLogLikelihood(model,par,...) + 2*par %*% Omega


## ## } else {

  
##     }
##    }

##   ### The actual minimization
                
##   model@optimResult <- optim(initPar,mll,gr=dmll,method="BFGS",control=control,...)
