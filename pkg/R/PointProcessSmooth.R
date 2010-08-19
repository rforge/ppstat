### TODO: The class PointProcessSmooth is going to change to the class where we do 
### spline based smoothing -- formally only "correct" with the identity phi. 
### Need to implement an automatic knot subset selection scheme, automatic (cross-
### validation or ?) penalty parameter selection.

setClass("PointProcessSmooth",
         contains="PointProcessModel"
         )

setMethod("initialize","PointProcessSmooth",
          function(.Object,
                   processData,
                   formula,
                   family,
                   Delta,
                   support,
                   storeBasis=TRUE,
                   fit=TRUE,
                   coefficients=NULL,
                   fixedCoefficients=list(),
                   modelMatrix=TRUE,
                   fisherInformation=modelMatrix,
                   Omega=NULL,
                   call=NULL,...){
            .Object@processDataEnv <- new.env(.GlobalEnv)
            .Object@processDataEnv$processData <- processData
            lockEnvironment(.Object@processDataEnv,bindings=TRUE)
            .Object@formula <- formula
            .Object@family <- family
            .Object@support <- support
            .Object@Delta <- Delta
            .Object@basisPoints <- seq(support[1],support[2],Delta)
            .Object@storeBasis <- storeBasis
            if(!is.null(Omega)) {
              if(any(Omega != t(Omega))) {
                Omega <- (Omega + t(Omega))/2
                warning("Penalization matrix Omega is not symmetric and is replaced by  (Omega + t(Omega))/2.")
              }
              if(!isTRUE(all.equal(min(eigen(Omega,only.values=TRUE,symmetric=TRUE)$values,0),0))) {
                stop("Penalization matrix Omega is not positive semi-definite.")
              }
              .Object@Omega <- Omega
              .Object@penalization <- TRUE
            } else if(is.null(Omega)) .Object@penalization <- FALSE
           
            .Object@delta <- as.numeric(unlist(tapply(getPosition(getContinuousProcess(getProcessData(.Object))),
                                                      getId(getContinuousProcess(getProcessData(.Object))),
                                                      function(x) c(diff(x),0)),use.names=FALSE))
            
            
            if(modelMatrix) {
              .Object <- computeModelMatrix(.Object)

              parDim <- dim(.Object@modelMatrix)[2]

              if(is.null(coefficients)){
                coefficients <- rep(0,parDim)
              } else {
                if(length(coefficients) != parDim) {
                  coefficients <- rep(0,parDim)
                  warning("Incorrect dimension of initial parameter vector. Initial parameters all set to 0")
                }
              }

              if(length(fixedCoefficients) != 0) coefficients[fixedCoefficients$which] <- fixedCoefficients$value
              .Object@fixedCoefficients <- fixedCoefficients
              .Object@coefficients <- coefficients
            
              if(fit){
                .Object <- glppmFit(.Object,coefficients,fisherInformation=fisherInformation,...)
              }
            }
            
            if(!is.null(call)) .Object@call <- call
            return(.Object)

          } 
        )
