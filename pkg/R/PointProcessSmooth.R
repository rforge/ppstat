pointProcessSmooth <- function(
                              formula,
                              data,
                              family,
                              support,
                              Delta,
                              basisPoints,
                              coefficients,
                              fixedCoefficients = list(),
                              fit = TRUE,
                              modelMatrix = TRUE,
                              varMethod = 'Fisher',
                              ...) {
  
  processDataEnv <- new.env(.GlobalEnv)
  processDataEnv$processData <- data
  lockEnvironment(processDataEnv,bindings=TRUE)

  modelMatrixEnv <- new.env(parent=.GlobalEnv)
  modelMatrixEnv$modelMatrix <- Matrix()

  ## TODO: Check if ... has a basisEnv, in which case that environment is used instead.
  basisEnv <- new.env(parent=.GlobalEnv)
  basisEnv$basis <- list()

  if(missing(Omega))
    {
      Omega <- matrix()  ### TODO: Here goes new stuff ...
      penalization <- TRUE
    }

  if(missing(basisPoints))
    {
      if(!(missing(support) & missing(Delta)))
        {
          if(length(support) == 1) support <- c(0,max(support[1],0))
          basisPoints <- seq(support[1],support[2],Delta)
        } else {
          stop("Must specify either 'support' and 'Delta' or 'basisPoints'.")
        }
    } else {
      support = range(basisPoints)
      Delta = min(diff(basisPoints))
    }
      
  delta <- as.numeric(unlist(tapply(getPosition(getContinuousProcess(data)),
                                    getId(getContinuousProcess(data)),
                                    function(x) c(diff(x),0)),use.names=FALSE))

  varMethods <- c('none','Fisher')
  if(!(varMethod %in% varMethods)){
    warning(paste("Method '", varMethod, "' for variance matrix estimation is currently not supported. Using method 'Fisher'.", sep=""))
    varMethod <- 'Fisher'
  }
  

  model <- new("PointProcessSmooth",
               processDataEnv = processDataEnv,
               delta = delta,
               formula = formula,
               family = family,
               call = match.call(),
               support = support,
               basisPoints = basisPoints,
               Delta = Delta,
               modelMatrixEnv = modelMatrixEnv,
               Omega = Omega,
               penalization = penalization,
               varMethod = varMethod,
               basisEnv = basisEnv,
               ...)
  
  if(modelMatrix) {
    model <- computeModelMatrix(model)
  }

  parDim <- dim(getModelMatrix(model))[2]
  
  if(missing(coefficients)){
    coefficients <- rep(0,parDim)
  } else {
    if(length(coefficients) != parDim) {
      coefficients <- rep(0,parDim)
      warning("Incorrect length of initial parameter vector. Initial parameters all set to 0")
    }
  }            

  if(length(fixedCoefficients) != 0) coefficients[fixedCoefficients$which] <- fixedCoefficients$value
  model@fixedCoefficients <- fixedCoefficients
  coefficients(model) <- coefficients
  
  if(fit){
    model <- ppmFit(model,...)
  }

  return(model)
}

## TODO: New summary function for an object of class 'PointProcessSmooth'.

setMethod("summary","PointProcessSmooth",
          function(object,...) {
            callNextMethod()
          }
          )
