setClass("PointProcessModel",
         representation(
                        modelMatrixEnv = "environment",  ### The modelMatrix as a 'Matrix' is in this environment. Locked after computation.
                        modelMatrixCol = "numeric",      ### The 'active' columns. Initially NULL, set in update and used in getModelMatrix.
                        coefficients = "numeric",
                        fixedCoefficients = "list",
                        Omega = "matrix",
                        penalization = "logical",
                        var = "matrix",
                        optimResult = "list",
                        basisEnv = "environment",      ### Evaluations of basis functions in support at Delta-grid values are in 
                                                       ### the list 'basis' in this environment.
                        basisPoints = "numeric"        ### The 'basisPoints' in contains the evaluation points for the basis functions.
                        ),
         validity = function(object) {
           if(isTRUE(object@penalization) && !isTRUE(all.equal(min(eigen(object@Omega,only.values=TRUE,symmetric=TRUE)$values,0),0)))
             stop("Penalization matrix 'Omega' is not positive semi-definite.")
           if(isTRUE(object@support[2] - object@support[1] <= 0))
             stop("Variable 'support' has to be either a positive number or a positive interval.")
           if(isTRUE(object@Delta > object@support[2] - object@support[1]))
             stop("Variable 'Delta' has to be smaller than the length of the support")
           return(TRUE)
         },
         contains="PointProcess"
         )


pointProcessModel <- function(
                              formula,
                              data,
                              family,
                              support,
                              Delta,
                              basisPoints,
                              Omega,
                              coefficients,
                              fixedCoefficients = list(),
                              fit = TRUE,
                              modelMatrix = TRUE,
                              fisherInformation = modelMatrix
                              ,...) {
  
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
      Omega <- matrix()
      penalization <- FALSE
    } else {
      if(any(Omega != t(Omega)))
        {
          Omega <- (Omega + t(Omega))/2
          warning("Penalization matrix Omega is not symmetric and is replaced by (Omega + t(Omega))/2.")
        }
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

  .Object <- new("PointProcessModel",
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
                 basisEnv = basisEnv,
                 ...)
                 
  if(modelMatrix) {
    .Object <- computeModelMatrix(.Object)
  }

  parDim <- dim(getModelMatrix(.Object))[2]
  
  if(missing(coefficients)){
    coefficients <- rep(0,parDim)
  } else {
    if(length(coefficients) != parDim) {
      coefficients <- rep(0,parDim)
      warning("Incorrect length of initial parameter vector. Initial parameters all set to 0")
    }
  }            

  if(length(fixedCoefficients) != 0) coefficients[fixedCoefficients$which] <- fixedCoefficients$value
  .Object@fixedCoefficients <- fixedCoefficients
  .Object@coefficients <- coefficients
            
  if(fit){
    .Object <- glppmFit(.Object,coefficients,fisherInformation=fisherInformation,...)
  }

  return(.Object)
}

pointProcessSmooth <- function(
                              formula,
                              data,
                              family,
                              support,
                              Delta,
                              basisPoints,
                              Omega,
                              coefficients,
                              fixedCoefficients = list(),
                              fit = TRUE,
                              modelMatrix = TRUE,
                              fisherInformation = modelMatrix
                              ,...) {
  
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
      Omega <- matrix()
      penalization <- FALSE
    } else {
      if(any(Omega != t(Omega)))
        {
          Omega <- (Omega + t(Omega))/2
          warning("Penalization matrix Omega is not symmetric and is replaced by (Omega + t(Omega))/2.")
        }
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

  .Object <- new("PointProcessModel",
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
                 basisEnv = basisEnv,
                 ...)
                 
  if(modelMatrix) {
    .Object <- computeModelMatrix(.Object)
  }

  parDim <- dim(getModelMatrix(.Object))[2]
  
  if(missing(coefficients)){
    coefficients <- rep(0,parDim)
  } else {
    if(length(coefficients) != parDim) {
      coefficients <- rep(0,parDim)
      warning("Incorrect length of initial parameter vector. Initial parameters all set to 0")
    }
  }            

  if(length(fixedCoefficients) != 0) coefficients[fixedCoefficients$which] <- fixedCoefficients$value
  .Object@fixedCoefficients <- fixedCoefficients
  .Object@coefficients <- coefficients
            
  if(fit){
    .Object <- glppmFit(.Object,coefficients,fisherInformation=fisherInformation,...)
  }

  return(.Object)
}

## TODO: Implement the stress-release type of model with
## intensity acccumulating according to time but "released"
## according to mark values.

## TODO: Implement the use of interaction terms.

computeBasis <- function(term,x,basisEnv,varLabels) {

  ## Basis evaluations for 'term' are computed if
  ## not already computed, locked and available in 'basis' in the
  ## 'basisEnv' environment.
    
  form <- as.formula(paste("~",term,"-1"))
  variables <- all.vars(form)
  if(missing(varLabels)) varLabels = variables
  
  if(all(variables %in% varLabels)) {
    if(environmentIsLocked(basisEnv)) return(TRUE)
    if(length(variables) > 1){
      stop(paste("Interactions of two or more variables in '",term,
                 "' is currently not supported.",sep=""))
    } else {
      x <- as.data.frame(x)
      colnames(x) <- variables
      basisEnv$basis[[term]] <- model.matrix(form,x)
      return(TRUE)
    }
  } else {
    return(FALSE)
  }
}

  

setMethod("subset","PointProcessModel",
          function(x,...){
            pointProcessModel(formula = x@formula,
                              data = subset(getProcessData(x),...),
                              family = x@family,
                              Delta = x@Delta,
                              support = x@support,
                              basisPoints = x@basisPoints,
                              Omega = x@Omega,
                              coefficients = x@coefficients,
                              fixedCoefficients = x@fixedCoefficients,
                              basisEnv = x@basisEnv
                              )
          }
          )



setMethod("getModelMatrix","PointProcessModel",
          function(object,...){
            if(length(object@modelMatrixCol) == 0) {
              return(object@modelMatrixEnv$modelMatrix)
            } else {
              return(object@modelMatrixEnv$modelMatrix[,object@modelMatrixCol])
            }
          }
          )


setMethod("coefficients","PointProcessModel",
          function(object,...){
            return(object@coefficients)
          }
          )

setReplaceMethod("coefficients","PointProcessModel",
                 function(.Object,value){
                   .Object@coefficients <- value
                   return(.Object)
                 }
                 )

setMethod("formula","PointProcessModel",
          function(x,...){
            return(x@formula)
          }
          )

setReplaceMethod("formula","PointProcessModel",
                 function(.Object,value){
                   .Object@formula <- value
                   .Object@basisEnv <- new.env(parent=.GlobalEnv)
                   basisEnv$basis <- list()
                   return(.Object)
                 }
                 )

setMethod("vcov","PointProcessModel",
          function(object,...){
            return(object@var)
          }
          )

setMethod("computeModelMatrix","PointProcessModel",
          function(object,evaluationPositions=NULL,...){

            ## The 'object' of class PointProcessModel contains the data
            ## as an object of class ProcessData and the formula for the
            ## model specification. The 'evalPositions' below corresponding to
            ## the model matrix rows are either given by the 'evaluationPositions'
            ## argument or extracted from the the ProcessData object (default).

            if(is.null(evaluationPositions)) {
              evalPositions <- tapply(getPosition(getContinuousProcess(getProcessData(object))),
                                      getId(getContinuousProcess(getProcessData(object))),list)
            } else {
              evalPositions <- evaluationPositions
            }
            
            ## The observed points ('positions') for the marked point process,
            ## the corresponding 'id' labels and 'marks' are extracted.
            
            markedPointProcess <- getMarkedPointProcess(getProcessData(object))
            positions <- getPosition(markedPointProcess)

            id <- factor(getId(markedPointProcess))
            idLevels <- levels(id)
            
            marks <- getMarkType(markedPointProcess)
            markLevels <- levels(marks)
            
            ## The formula object is extracted and decomposed into terms.
            ## Each term label (in 'termLabels') is processed below,
            ## and the corresponding columns in the model matrix are computed.
            
            mt <- terms(object@formula) 
            termLabels <- attr(mt,"term.labels")
            notMarkTerms <- character()
            
            ## The points where the basis functions are evaluated are extracted
            ## and the list of model matrices ('design') is set up, which holds model
            ## matrices for the different terms. 
            
            design <- list()
            
            ## Model matrix computations for the terms involving the marks.
            ## Terms involving 'id' and the continuous process components below.
            ## Those terms are collected in 'notMarkTerms' in the loop below.
            
            for(term in termLabels) {
              if(computeBasis(term,object@basisPoints,object@basisEnv,markLevels))
                {
                ## The call to computeBasis is invoked for its side effect
                ## of computing the basis evaluations if that is not already
                ## done (which is checked by checking if the environment 'basisEnv' 
                ## is locked). The function also checks if 'term' holds the correct
                ## number of variables from markLevels.

                  designList <- list()

                ## Model matrix computed for each value of 'id' and stored in
                ## 'designList'. 

                ## Central loop over 'idLevels' and computations of the model
                ## matrix in the C function 'computeModelMatrix'. Result is
                ## converted to a sparse matrix, bound together in one matrix
                ## below and stored in the list 'design'.
                
                ## TODO: Can the C level computation return a sparse matrix
                ## directly?

                for(i in idLevels) {
                  mark <- all.vars(parse(text=term))
                  posi <- sort(positions[marks == mark & id == i])
                  ## Is this the right place to order the observed positions?
                  designList[[i]] <- Matrix(.Call("computeModelMatrix",
                                                  evalPositions[[i]],
                                                  object@basisEnv$basis[[term]],
                                                  object@Delta,
                                                  posi,
                                                  PACKAGE="ppstat"),sparse=TRUE)
                }
                design[[term]] <- do.call("rBind",designList)
                colnames(design[[term]]) <- colnames(object@basisEnv$basis[[term]])

              } else {
                ## The term does not involve marks
                notMarkTerms <- c(notMarkTerms,term) 
              }
            }
            
            ## Is there an intercept in the model? If so, add the intercept
            ## explicitly to the vector 'notMarkType'.

            if(attr(mt,"intercept") == 1){
              notMarkTerms <- paste(c(notMarkTerms,"1"),collapse="+")
            } else if(length(notMarkTerms) > 0) {
              notMarkTerms <-  paste(paste(notMarkTerms,collapse="+"),"-1")
            }

            ## Model matrix computations for terms involving 'id' and
            ## continuous time process components.

            if(length(notMarkTerms) > 0){
              form <-  as.formula(paste("~",notMarkTerms))
              variables <- all.vars(form)

              if(all(variables %in% c("id",colnames(getValue(getContinuousProcess(getProcessData(object))))))) {
                values <- data.frame(id=getId(getContinuousProcess(getProcessData(object))))
                notIdVariables <- variables[variables != "id"]

                if(length(notIdVariables) > 0) {
                  tmp <- as.matrix(getValue(getContinuousProcess(getProcessData(object)))[,notIdVariables,drop=FALSE])
                  rownames(tmp) <- rownames(values)
                  values <- cbind(values,tmp)
                }

                tmp <- model.matrix(form,values)
                X0 <- Matrix(tmp,dimnames=dimnames(tmp),sparse=TRUE)
              } else {
                stop(paste("Use of non existing variable(s) in:", form))
              }
            }

            if(environmentIsLocked(object@modelMatrixEnv)) object@modelMatrixEnv <- new.env(parent=.GlobalEnv)
            object@modelMatrixEnv$modelMatrix <- cBind(X0,do.call("cBind",design))
            lockEnvironment(object@modelMatrixEnv,binding=TRUE)
            lockEnvironment(object@basisEnv,binding=TRUE)            
            return(object)
          }
          )

            
setMethod("computeLinearPredictor","PointProcessModel",
          function(object,coefficients=NULL,...){
            if(is.null(coefficients)) {
              coefficients <- coefficients(object)
            } 
                                     
            eta =  as.numeric(getModelMatrix(object) %*% coefficients)
            return(eta)
          }
          )

setMethod("predict","PointProcessModel",
          function(object,...) {
            eta <- computeLinearPredictor(object,...)
            return(object@family@phi(eta))
          }
          )

setMethod("computeDMinusLogLikelihood","PointProcessModel",
          function(object,coefficients=NULL,...){
            eta <- computeLinearPredictor(object,coefficients,...)
             if(attr(terms(object@formula),"response") != 0) {
              response <- all.vars(object@formula,unique=FALSE)[attr(terms(object@formula),"response")]
            } else stop("no response variable specified")

            if(object@family@link == "log") {

              dmll <- as.vector(t(exp(eta)*object@delta)%*%getModelMatrix(object)) -
                colSums(getModelMatrix(object)[getMarkTypePosition(getProcessData(object),response),])

            } else {
              
              etaP <- eta[getMarkTypePosition(getProcessData(object),response)]
              mmP <- getModelMatrix(object)[getMarkTypePosition(getProcessData(object),response),]

              dmll <-  as.vector(t(object@family@Dphi(eta)*object@delta)%*%getModelMatrix(object)) -
                as.vector(t(object@family@Dphi(etaP)/object@family@phi(etaP))%*%mmP)
              
            }
            
            return(dmll)
            
          }
          )



setMethod("computeDDMinusLogLikelihood","PointProcessModel",
          function(object,coefficients=NULL,...){
            eta <- computeLinearPredictor(object,coefficients,...)
             if(attr(terms(object@formula),"response") != 0) {
              response <- all.vars(object@formula,unique=FALSE)[attr(terms(object@formula),"response")]
            } else stop("no response variable specified")

             if(object@family@link == "log"){

               ddmll <-  as(crossprod(getModelMatrix(object),exp(eta)*object@delta*getModelMatrix(object)),"matrix")

             } else if(object@family@link == "identity"){

               etaP <- eta[getMarkTypePosition(getProcessData(object),response)]
               mmP <- getModelMatrix(object)[getMarkTypePosition(getProcessData(object),response),]

               ddmll <-  as(crossprod(mmP,1/object@family@phi(etaP)^2*mmP),"matrix")

             } else {

               etaP <- eta[getMarkTypePosition(getProcessData(object),response)]
               mmP <- getModelMatrix(object)[getMarkTypePosition(getProcessData(object),response),]

               ddmll <-  as(crossprod(getModelMatrix(object),object@family@D2phi(eta)*object@delta*getModelMatrix(object)),"matrix") -
                 as(crossprod(mmP,(object@family@D2phi(etaP)*object@family@phi(etaP) - object@family@Dphi(etaP)^2)/object@family@phi(etaP)^2*mmP),"matrix")

             }
            
            return(ddmll)
            
          }
          )


setMethod("update","PointProcessModel",
          function(object,...){
            .local <- function(object,formula,fisherInformation=TRUE,warmStart=TRUE,fixedPar=list(),...){

              updatedFormula <- update(formula(object),formula)
              formula(object) <- updatedFormula

              updatedTerms <- attr(terms(updatedFormula),"term.labels")
              if(attr(terms(updatedFormula),"intercept")==1) updatedTerms <- c(updatedTerms,"Intercept")

              ## TODO: The use of grep might be improved to get a better update function
              
              updatedModelMatrixColumns <- sapply(updatedTerms,grep,colnames(getModelMatrix(object)),fixed=TRUE,simplify=FALSE)
                            
              if(any(sapply(updatedModelMatrixColumns,function(x) length(x) == 0))) {
                object <- computeModelMatrix(object)
                coefficients(object) <- rep(0,dim(getModelMatrix(object))[2])
              } else {
                col <- sort(unlist(updatedModelMatrixColumns))
                if(any(attr(terms(object@formula),"order") >= 2)) {
                  warning("Original model specification includes interaction terms. The updated model may not be as desired.")
                }
                object@modelMatrixCol <- col
                if(warmStart) object@coefficients <- coefficients(object)[col]
              }

              if(length(fixedPar) != 0) object@coefficients[fixedPar$which] <- fixedPar$value
              object@fixedCoefficients <- fixedPar
                          
              return(glppmFit(object,initPar=coefficients(object),fisherInformation=fisherInformation,...))
            }
            object@call <- match.call()
            .local(object,...)
          }
          )

setMethod("getLinearFilter",c(object="PointProcessModel",se="logical"),
          function(object,se=FALSE,nr=NULL...){
            
            markLevels <- levels(getMarkType(getMarkedPointProcess(getProcessData(object))))
            termLabels <- attr(terms(formula(object)),"term.labels")

            linearFilter <- list()
            if(isTRUE(se)) linearFilterSE <- list()
            if(is.null(nr)) nr <- dim(object@basisEnv$basis[[term]])[1]
            
            for(term in termLabels){
              if(computeBasis(term,object@basisPoints,object@basisEnv,markLevels))
                {
                  i <- seq_len(nr)*floor(dim(object@basisEnv$basis[[term]])[1]/nr)
                  design <- object@basisEnv$basis[[term]][i,]
                  linearFilter[[term]] <- design %*% coefficients(object)[dimnames(design)[[2]]]
                  if(isTRUE(se))
                    linearFilterSE[[term]] <- sqrt(rowSums(design %*% vcov(object)[dimnames(design)[[2]],dimnames(design)[[2]]] * design))
                }
            }
            
            lockEnvironment(object@basisEnv,bindings=TRUE)
            if(isTRUE(se)) {
              return(list(linearFilter=cbind(data.frame(x=object@basisPoints[i]), as.data.frame(linearFilter)),se=linearFilterSE))
            } else {
              return(cbind(data.frame(x=object@basisPoints[i]), as.data.frame(linearFilter)))
            }
              
          }
          )

setMethod("plot",signature(x="PointProcessModel",y="missing"),
          function(x,y,trans=NULL,alpha=0.05,...){
            linearFilter <- getLinearFilter(x,se=TRUE,nr=1000)
            moltenFilter <- melt(linearFilter$linearFilter,id.vars="x")
            q <- qnorm(1-alpha/2)
            plotData <- cbind(moltenFilter,data.frame(cf.lower = moltenFilter$value-q*unlist(linearFilter$se),
                                                      cf.upper = moltenFilter$value+q*unlist(linearFilter$se)))
                              
            if(!is.null(trans)) plotData[,c("value","cf.lower","cf.upper")] <- do.call(trans,list(plotData[,c("value","cf.lower","cf.upper")]))

            linearFilterPlot <- ggplot(data=plotData,aes(x=x,y=value)) +
              facet_grid(.~variable) +
                geom_ribbon(aes(min=cf.lower,max=cf.upper),fill=alpha("blue",0.2)) +
                  scale_x_continuous("position") +
                    scale_y_continuous("") + geom_line()

            return(linearFilterPlot)
          }
          )

setMethod("print","PointProcessModel",
          function(x,digits= max(3, getOption("digits") - 3), ...){
            cat("\nCall:\n", deparse(x@call), "\n\n", sep = "")
            if (length(coefficients(x))) {
              cat("Coefficients:\n")
              print.default(format(coefficients(x), digits = digits), print.gap = 2, 
                            quote = FALSE)
            }
            else cat("No coefficients\n")
            cat("\n")
            invisible(x)
          }
          )

setMethod("show","PointProcessModel",
          function(object) print(x=object)
          )


### TODO: This should return an S4 object instead with appropriate view method.

setMethod("summary","PointProcessModel",
          function(object,...) {
            result <- list()

            result$df <- length(object@coefficients) - length(object@fixedCoefficients$which)
            result$call <- object@call
            result$mll <- object@optimResult$value
            result$iter <- object@optimResult$counts
            result$convergence <- object@optimResult$convergence
            result$aic <- 2*(result$mll +  result$df)

            # Some way of summarizing 'residuals'

            se <- sqrt(diag(object@var))
            z <- object@coefficients/se
            q <- 2*(1-pnorm(abs(z)))
            z[se==0] <- NA
            q[se==0] <- NA
            se[se==0] <- NA
            
            result$coefficients <- matrix(c(object@coefficients,se,z,q),ncol=4)
            zapNames <- sapply(strsplit(names(object@coefficients),""),
                               function(x) {
                                 end <- x[seq(max(length(x)-5,1),length(x))]
                                 if(length(x) > 18) end <- c("...",x[seq(max(length(x)-2,1),length(x))])
                                 paste(c(x[1:min(max(length(x)-6,1),12)],end),sep="",collapse="")
                               })
                             
            
            dimnames(result$coefficients) <- list(zapNames,c("Estimate","Std. Error", "z value", "Pr(> |z|)"))

            result$fixedCoefficients <- object@fixedCoefficients
            result$penalization <- object@penalization
            class(result) <- c("summary.glppm")
            return(result)

          }
          )



            
