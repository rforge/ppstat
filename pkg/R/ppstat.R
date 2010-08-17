###======================================================================
### ppstat - point process statistics
###======================================================================

## bSpline replaces bs for convenience and provides easier computations
## of B-spline bases that we are using. Does not issue a warning when
## evaluated outside of knots. This will in fact be the rule rather
## than the exception.

bSpline <- function(x,knots,...,trunc=0,whichColumns=NULL){
  if(length(knots) <= 4) stop("Need at least 5 knots")
  design <- splineDesign(knots,x,...,outer.ok=TRUE)
  if(!is.null(whichColumns)) design <- design[,whichColumns]
  if(!is.null(trunc)) design[x <= trunc,] <- 0
  return(design)
}

## Truncated exponential

tExp <- function(x,tLevel=1e-4){
  y <- exp(x)
  y[y < tLevel] <- 0
  return(y)
}

## Constant function

const <- function(x,y,c=1){
  c*as.numeric(x<=y)
}

     
## Implementation of main function to call for the estimation of
## generalized linear point process models.

glppm <- function(formula,data,family,pointProcessModel,Delta=0.1,support=c(0,1),initPar=NULL,fisherInformation=TRUE,fixedPar=list(),Omega=NULL,...){

  ## At this point it may be desirable to make certain checks and to parse some
  ## user specified information. For instance a choice of evaluation positions
  ## or specification of intelligent evaluation schemes.
  ## This has at the moment been moved to the class ProcessData

  if(length(support) == 1) support <- c(0,support)
  if(support[2] - support[1] <= 0) stop("support has to be either a positive number or a positive interval")
  if(Delta > support[2] - support[1]) stop("Delta has to be smaller than the length of the support")
  
  (function(fisherPreComp=NULL,...) {
    if(!is.null(fisherPreComp)) warning("Use of precomputations and argument 'fisherPreComp' in 'glppm' is deprecated",call.=FALSE)
  })(...)

  if(missing(pointProcessModel)){
    if(!(class(data)=="ProcessData")) stop("Function 'glppm' needs data to be of class 'ProcessData'")
    pointProcessModel <- new("PointProcessModel",
                             processData=data,
                             formula=formula,
                             family=family,
                             coefficients=initPar,
                             fixedCoefficients=fixedPar,
                             fisherInformation=fisherInformation,
                             Omega=Omega,
                             support=support,
                             Delta=Delta, 
                             call=match.call(),...)
  } else if(!(class(pointProcessModel)=="PointProcessModel")) {
    stop("Function 'glppm' needs a pointProcessModel of class 'PointProcessModel'")
  }
  return(pointProcessModel)
}

glppmFit <- function(model,initPar=NULL,fisherInformation=TRUE,control=list(),...) {
  nrPar <- parDim <- dim(getModelMatrix(model))[2]
  fixedPar <- model@fixedCoefficients
  Omega <- model@Omega

  control <- c(list(maxit=1000),control)
  
  if(length(fixedPar) != 0) {
    nrPar <- parDim -length(fixedPar$which)
    tmpPar <- numeric(parDim)
    tmpPar[fixedPar$which] <- fixedPar$value
  }

  ### Setting up the initial parameters

  if(is.null(initPar)){
    if(length(coefficients(model)) == parDim) {
      if(length(fixedPar) != 0) {
        initPar <- coefficients(model)[-fixedPar$which]
      } else {
        initPar <- coefficients(model)
      }      
    } else {
      initPar <- rep(0,nrPar)
    }
  } else { #initPar is not NULL
    if(length(initPar) != parDim) {
      initPar <- rep(0,nrPar)
      warning("Incorrect length of initial parameter vector. Initial parameters all set to 0.")
    } else {
      if(length(fixedPar) != 0) { 
        initPar <- initPar[-fixedPar$which]
      }
    }
  }

  ### Setting up the objective function to minimize

  if(!model@penalization){
    if(length(fixedPar) == 0) {
      mll <- function(par,...) computeMinusLogLikelihood(model,par,...)
      dmll <- function(par,...) computeDMinusLogLikelihood(model,par,...)
    } else {
      mll <- function(par,...) {
        tmpPar[-fixedPar$which] <- par
        computeMinusLogLikelihood(model,tmpPar,...)
      }
      dmll <- function(par,...) {
        tmpPar[-fixedPar$which] <- par
        computeDMinusLogLikelihood(model,tmpPar,...)[-fixedPar$which]
      }
    }
  } else {
     if(length(fixedPar) == 0) {
      mll <- function(par,...) computeMinusLogLikelihood(model,par,...) + par %*% Omega %*% par
      dmll <- function(par,...) computeDMinusLogLikelihood(model,par,...) + 2*par %*% Omega
    } else {
      mll <- function(par,...) {
        tmpPar[-fixedPar$which] <- par
        computeMinusLogLikelihood(model,tmpPar,...) + par %*% Omega[-fixedPar$which,-fixedPar$which] %*% par
      }
      dmll <- function(par,...) {
        tmpPar[-fixedPar$which] <- par
        computeDMinusLogLikelihood(model,tmpPar,...)[-fixedPar$which] + 2*par %*% Omega[-fixedPar$which,-fixedPar$which]
      }
    }
   }

  ### The actual minimization
                
  model@optimResult <- optim(initPar,mll,gr=dmll,method="BFGS",control=control,...)

  if(length(fixedPar) == 0) {
    model@coefficients <- model@optimResult$par
  } else {
    model@coefficients[-fixedPar$which] <- model@optimResult$par
    model@coefficients[fixedPar$which] <- fixedPar$value
  }
  names(model@coefficients) <- dimnames(model@modelMatrix)[[2]]

  ### Computation of the estimated covariance matrix
  if(fisherInformation) {
    vcovInv <- computeDDMinusLogLikelihood(model)
    if(model@penalization) vcovInv <- vcovInv + 2*Omega ## This requires some more thought ....
    vcov <- matrix(0,nrow=dim(vcovInv)[1],ncol=dim(vcovInv)[2])
    
    if(length(fixedPar) == 0) {
      tmp <- try(solve(vcovInv),silent=TRUE)
      if(class(tmp)=="try-error") {
        cat("Covariance matrix can not be estimated.",tmp[1]," Check convergence status or parameterization\n.")
      } else {
        vcov <- tmp
      }
    } else {
      tmp <- try(solve(vcovInv[-fixedPar$which,-fixedPar$which]),silent=TRUE)
      if(class(tmp)=="try-error") {
        cat("Fisher information singular:\n",tmp[1],"\nCheck convergence status and parameterization.")
      } else {
        vcov[-fixedPar$which,-fixedPar$which] <- tmp
      }
    }
    rownames(vcov) <- names(model@coefficients)
    colnames(vcov) <- names(model@coefficients)
    model@var <- (vcov + t(vcov))/2   ## To assure symmetry
  } else {
    model@var <- matrix(0,length(model@coefficients),length(model@coefficients))
  }
  return(model)
}

### Todo: This should be changed to an S4 method and specialized to the
### different point process classes.

print.summary.glppm <-  function (x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    if(length(x$fixedCoefficients$which) == 1) {
      cat("\nNote that", rownames(x$coefficients)[x$fixedCoefficients$which], "was fixed\n\n")
    } else if(length(x$fixedCoefficients$which) > 1) {
      cat("\nNote that", paste(rownames(x$coefficients)[x$fixedCoefficients$which],collapse=", "), "were fixed\n\n")
    }
    if(x$penalization) cat("\nWarning: Parameters estimated with a quadratic penalty term.\nCorrections of the degrees of freedom is currently not implemented\nand AIC is not corrected either. The estimated standard errors may be\ninaccurate and the univariate test statistics can be misleading.\n\n")
    cat("\nMinus-log-likelihood: ", format(x$mll, digits = max(4, digits + 1)), " on ",x$df, " degrees of freedom",
        "\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), 
        "\n\n", "Number of function evaluations: ", x$iter[1],
        "\n","Number of gradient evaluations: ", x$iter[2],
        "\n", sep = "")
    if(x$convergence!=0) cat("\nWarning: Algorithm did not converge.\n         'optim' convergence status:",x$convergence,"\n")  
    cat("\n")
    invisible(x)
}



### Stub for generalized additive point process models
### with smoothing. Have to figure out where to automatically
### place the knots - any theory here?
### 10/3-2010 This stub is deprecated and will be removed. 

### gappm <- function(formula,data){
###  if(!(class(data)=="PointProcessData")) stop("Function 'gappm' needs data to be of class 'PointProcessData'")
###
### }


## Stub for vector generalized linear point process model.
## First planned implementation is a matter of exploiting variation
## independence of parameters, which allows for the estimation
## to be broken down

### vgappm <- function(formula,data){
###  if(!(class(data)=="PointProcessData")) stop("Function 'gappm' needs data to be of class 'PointProcessData'")
###
### }


### This is the main function to be called to estimate a generalized linear
### point process smooth.

glpps <- function(formula,data,family,pointProcessSmooth,Delta=0.1,support=c(0,1),initSmooth=NULL,Omega=NULL,...){

  if(length(support) == 1) support <- c(0,support)
  if(support[2] - support[1] <= 0) stop("support has to be either a positive number or a positive interval")
  if(Delta > support[2] - support[1]) stop("Delta has to be smaller than the length of the support")
  
  if(missing(pointProcessModel)){
    if(!(class(data)=="ProcessData")) stop("Function 'glppm' needs data to be of class 'ProcessData'")
    pointProcessSmooth <- new("PointProcessSmooth",
                              processData=data,
                              formula=formula,
                              family=family,
                              support=support,
                              g=initSmooth,
                              Omega=Omega,
                              call=match.call(),...)
  } else if(!(class(pointProcessModel)=="PointProcessSmooth")) {
    stop("Function 'glpps' needs a pointProcessModel of class 'PointProcessModel'")
  }
  return(pointProcessModel)
}


glppsFit <- function(model,initSmooth=NULL,control=list(),...)
  {
  }
