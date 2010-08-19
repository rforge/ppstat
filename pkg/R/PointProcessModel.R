setClass("PointProcessModel",
         representation(
                        modelMatrixEnv = "environment",  ### The modelMatrix as a 'Matrix' is in this environment. Locked after computation.
                        modelMatrixCol = "numeric",      ### The 'active' columns. Set in update, used in getModelMatrix and reset in computeModelMAtrix
                        coefficients = "numeric",
                        fixedCoefficients = "list",
                        Omega = "matrix",
                        penalization = "logical",
                        var = "matrix",
                        varMethod = "character",            ### Which method is used to compute the estimate of the variance. 'pointProcessModel' has default 'Fisher'.
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
                              varMethod = "Fisher",
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

  model <- new("PointProcessModel",
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
  model@coefficients <- coefficients
            
  if(fit){
    model <- glppmFit(model,...)
  }

  return(model)
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
                              varMethod = "Fisher",
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

  model <- new("PointProcessModel",
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
  model@coefficients <- coefficients
  
  if(fit){
    model <- glppmFit(model,...)
  }

  return(model)
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
          function(model,...){
            if(length(model@modelMatrixCol) == 0) {
              return(model@modelMatrixEnv$modelMatrix)
            } else {
              return(model@modelMatrixEnv$modelMatrix[,model@modelMatrixCol])
            }
          }
          )

## TODO: Implement setReplaceMethod for getModelMatrix. 


setMethod("coefficients","PointProcessModel",
          function(object,...){
            return(object@coefficients)
          }
          )

setReplaceMethod("coefficients",c(model="PointProcessModel",value="numeric"),
                 function(model,value){
                   model@coefficients <- value
                   return(model)
                 }
                 )

setMethod("vcov","PointProcessModel",
          function(object,...){
            attr(object@var,"method") <- object@varMethod
            return(object@var)
          }
          )

setMethod("computeModelMatrix","PointProcessModel",
          function(model,evaluationPositions=NULL,...){

            ## Resetting the selected columns
            
            model@modelMatrixCol <- numeric()

            ## The 'model' of class PointProcessModel contains the data
            ## as an object of class ProcessData and the formula for the
            ## model specification. The 'evalPositions' below corresponding to
            ## the model matrix rows are either given by the 'evaluationPositions'
            ## argument or extracted from the the ProcessData object (default).

            if(is.null(evaluationPositions)) {
              evalPositions <- tapply(getPosition(getContinuousProcess(getProcessData(model))),
                                      getId(getContinuousProcess(getProcessData(model))),list)
            } else {
              evalPositions <- evaluationPositions
            }
            
            ## The observed points ('positions') for the marked point process,
            ## the corresponding 'id' labels and 'marks' are extracted.
            
            markedPointProcess <- getMarkedPointProcess(getProcessData(model))
            positions <- getPosition(markedPointProcess)

            id <- factor(getId(markedPointProcess))
            idLevels <- levels(id)
            
            marks <- getMarkType(markedPointProcess)
            markLevels <- levels(marks)
            
            ## The formula object is extracted and decomposed into terms.
            ## Each term label (in 'termLabels') is processed below,
            ## and the corresponding columns in the model matrix are computed.
            
            mt <- terms(model@formula) 
            termLabels <- attr(mt,"term.labels")
            notMarkTerms <- character()
            
            ## The points where the basis functions are evaluated are extracted
            ## and the list of model matrices ('design') is set up, which holds model
            ## matrices for the different terms. 'assign' will be an attribute to  
            ## the model matrix of length equal to the number of columns, and for
            ## each column pointing to the term number. 
            
            design <- list()
            assign <- numeric()
            
            ## Model matrix computations for the terms involving the marks.
            ## Terms involving 'id' and the continuous process components below.
            ## Those terms are collected in 'notMarkTerms' in the loop below.
            
            for(term in termLabels) {
              if(computeBasis(term,model@basisPoints,model@basisEnv,markLevels))
                {
                  termNr <- which(term == termLabels)
                  assign <- c(assign,rep(termNr,dim(model@basisEnv$basis[[term]])[2])) 
                  
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
                                                  model@basisEnv$basis[[term]],
                                                  model@Delta,
                                                  posi,
                                                  PACKAGE="ppstat"),sparse=TRUE)
                }
                design[[term]] <- do.call("rBind",designList)
                colnames(design[[term]]) <- colnames(model@basisEnv$basis[[term]])

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

              if(all(variables %in% c("id",colnames(getValue(getContinuousProcess(getProcessData(model))))))) {
                values <- data.frame(id=getId(getContinuousProcess(getProcessData(model))))
                notIdVariables <- variables[variables != "id"]

                if(length(notIdVariables) > 0) {
                  tmp <- as.matrix(getValue(getContinuousProcess(getProcessData(model)))[,notIdVariables,drop=FALSE])
                  rownames(tmp) <- rownames(values)
                  values <- cbind(values,tmp)
                }

                tmp <- model.matrix(form,values)
                termPos <- c(0,sapply(attr(terms(form),"term.labels"),function(t) which(t == termLabels),USE.NAMES=FALSE))
                assign <- c(termPos[attr(tmp,"assign")+1],assign)
                X0 <- Matrix(tmp,dimnames=dimnames(tmp),sparse=TRUE)
              } else {
                stop(paste("Use of non existing variable(s) in:", form))
              }
            }

            if(environmentIsLocked(model@modelMatrixEnv)) model@modelMatrixEnv <- new.env(parent=.GlobalEnv)
            model@modelMatrixEnv$modelMatrix <- cBind(X0,do.call("cBind",design))
            attr(model@modelMatrixEnv$modelMatrix,"assign") <- assign
            lockEnvironment(model@modelMatrixEnv,binding=TRUE)
            lockEnvironment(model@basisEnv,binding=TRUE)            
            return(model)
          }
          )

            
setMethod("computeLinearPredictor","PointProcessModel",
          function(model,coefficients=NULL,...){
            if(is.null(coefficients)) {
              coefficients <- coefficients(model)
            } 
                                     
            eta =  as.numeric(getModelMatrix(model) %*% coefficients)
            return(eta)
          }
          )

setMethod("predict","PointProcessModel",
          function(object,...) {
            eta <- computeLinearPredictor(object,...)
            return(model@family@phi(eta))
          }
          )

setMethod("computeDMinusLogLikelihood","PointProcessModel",
          function(model,coefficients=NULL,...){
            eta <- computeLinearPredictor(model,coefficients,...)
             if(attr(terms(model@formula),"response") != 0) {
              response <- all.vars(model@formula,unique=FALSE)[attr(terms(model@formula),"response")]
            } else stop("no response variable specified")

            if(model@family@link == "log") {

              dmll <- as.vector(t(exp(eta)*model@delta)%*%getModelMatrix(model)) -
                colSums(getModelMatrix(model)[getMarkTypePosition(getProcessData(model),response),])

            } else {
              
              etaP <- eta[getMarkTypePosition(getProcessData(model),response)]
              mmP <- getModelMatrix(model)[getMarkTypePosition(getProcessData(model),response),]

              dmll <-  as.vector(t(model@family@Dphi(eta)*model@delta)%*%getModelMatrix(model)) -
                as.vector(t(model@family@Dphi(etaP)/model@family@phi(etaP))%*%mmP)
              
            }
            
            return(dmll)
            
          }
          )



setMethod("computeDDMinusLogLikelihood","PointProcessModel",
          function(model,coefficients=NULL,...){
            eta <- computeLinearPredictor(model,coefficients,...)
             if(attr(terms(model@formula),"response") != 0) {
              response <- all.vars(model@formula,unique=FALSE)[attr(terms(model@formula),"response")]
            } else stop("no response variable specified")

             if(model@family@link == "log"){

               ddmll <-  as(crossprod(getModelMatrix(model),exp(eta)*model@delta*getModelMatrix(model)),"matrix")

             } else if(model@family@link == "identity"){

               etaP <- eta[getMarkTypePosition(getProcessData(model),response)]
               mmP <- getModelMatrix(model)[getMarkTypePosition(getProcessData(model),response),]

               ddmll <-  as(crossprod(mmP,1/model@family@phi(etaP)^2*mmP),"matrix")

             } else {

               etaP <- eta[getMarkTypePosition(getProcessData(model),response)]
               mmP <- getModelMatrix(model)[getMarkTypePosition(getProcessData(model),response),]

               ddmll <-  as(crossprod(getModelMatrix(model),model@family@D2phi(eta)*model@delta*getModelMatrix(model)),"matrix") -
                 as(crossprod(mmP,(model@family@D2phi(etaP)*model@family@phi(etaP) - model@family@Dphi(etaP)^2)/model@family@phi(etaP)^2*mmP),"matrix")

             }
            
            return(ddmll)
            
          }
          )


setMethod("update","PointProcessModel",
          function(object,...){
            .local <- function(model,formula,warmStart = TRUE,fixedCoefficients = list(),...){

              termLabels <- attr(terms(formula(model)),"term.labels")
              updatedFormula <- update(formula(model),formula)
              updatedTermLabels <- attr(terms(updatedFormula),"term.labels")
              if(attr(terms(formula(model)),"intercept")==1) termLabels <- c("Intercept",termLabels)
              if(attr(terms(updatedFormula),"intercept")==1) updatedTermLabel <- c("Intercept",updatedTermLabels)
              
              if(all(updatedTermLabels %in% termLabels)) {
                col <- which(termLabels[attr(getModelMatrix(tmpPPM),"assign")+1] %in% updatedTermLabels)
                if(any(attr(terms(formula(model)),"order") >= 2)) {
                  warning(paste(c(object@call,"Original model formula includes interaction terms. Check that updated model is as expected."),collapse="\n"),call.=FALSE)
                }
                formula(model) <- updatedFormula
                model@modelMatrixCol <- col
                if(warmStart) coefficients(model) <- coefficients(model)[col] 
              } else {
                formula(model) <- updatedFormula
                model <- computeModelMatrix(model)
                coefficients(model) <- rep(0,dim(getModelMatrix(model))[2])
              }

              if(length(fixedCoefficients) != 0) model@coefficients[fixedCoefficients$which] <- fixedCoefficients$value
              model@fixedCoefficients <- fixedCoefficients
                          
              return(glppmFit(model,...))
            }
            object@call <- match.call()
            .local(object,...)
          }
          )

setMethod("getLinearFilter",c(model="PointProcessModel",se="logical"),
          function(model,se=FALSE,nr=NULL...){
            
            markLevels <- levels(getMarkType(getMarkedPointProcess(getProcessData(model))))
            termLabels <- attr(terms(formula(model)),"term.labels")

            linearFilter <- list()
            if(isTRUE(se)) linearFilterSE <- list()
            if(is.null(nr)) nr <- dim(model@basisEnv$basis[[term]])[1]
            
            for(term in termLabels){
              if(computeBasis(term,model@basisPoints,model@basisEnv,markLevels))
                {
                  i <- seq_len(nr)*floor(dim(model@basisEnv$basis[[term]])[1]/nr)
                  design <- model@basisEnv$basis[[term]][i,]
                  linearFilter[[term]] <- design %*% coefficients(model)[dimnames(design)[[2]]]
                  if(isTRUE(se))
                    linearFilterSE[[term]] <- sqrt(rowSums(design %*% vcov(model)[dimnames(design)[[2]],dimnames(design)[[2]]] * design))
                }
            }
            
            lockEnvironment(model@basisEnv,bindings=TRUE)
            if(isTRUE(se)) {
              return(list(linearFilter=cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)),se=linearFilterSE))
            } else {
              return(cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)))
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
          function(model) print(x=model)
          )


### TODO: This should return an S4 object instead with appropriate view method.
### TODO: Implementation of standard model diagnostics.

setMethod("summary","PointProcessModel",
          function(object,...) {
            result <- list()

            result$df <- length(object@coefficients) - length(object@fixedCoefficients$which)
            result$call <- object@call
            result$mll <- object@optimResult$value
            result$iter <- object@optimResult$counts
            result$convergence <- object@optimResult$convergence
            result$aic <- 2*(result$mll +  result$df)

            ## TODO: Implement some way of summarizing 'residuals'

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

setMethod("computeVar","PointProcessModel",
          function(model,...){
            if(attr(vcov(model),"method") == "Fisher") {
              vcovInv <- computeDDMinusLogLikelihood(model)
              if(model@penalization) vcovInv <- vcovInv + 2*Omega ## TODO: This requires some more thought ....
              vcov <- matrix(0,nrow=dim(vcovInv)[1],ncol=dim(vcovInv)[2])
              
              if(length(model@fixedCoefficients) == 0) {
                tmp <- try(solve(vcovInv),silent=TRUE)
                if(class(tmp)=="try-error") {
                  cat("Fisher information singular:\n",tmp[1]," Check convergence status or parameterization\n.")
                } else {
                  vcov <- tmp
                }
              } else {
                tmp <- try(solve(vcovInv[-model@fixedCoefficients$which,-model@fixedCoefficients$which]),silent=TRUE)
                if(class(tmp)=="try-error") {
                  cat("Fisher information singular:\n",tmp[1],"\nCheck convergence status and parameterization.")
                } else {
                  vcov[-model@fixedCoefficients$which,-model@fixedCoefficients$which] <- tmp
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
          )
                   
setMethod("glppmFit","PointProcessModel",
          function(model,control=list(),...) {
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
            
            if(length(coefficients(model)) == parDim) {
              if(length(fixedPar) != 0) {
                initPar <- coefficients(model)[-fixedPar$which]
              } else {
                initPar <- coefficients(model)
              }    
            } else {
              initPar <- rep(0,nrPar)
              warning("Length of initial parameter vector worng. Initial parameters all set to 0.")
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
            names(model@coefficients) <- dimnames(getModelMatrix(model))[[2]]

### Computation of the estimated covariance matrix

            model <- computeVar(model)
            return(model)
          }
          )            
