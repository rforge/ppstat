setClass("PointProcessModel",
         representation(
                        ## Evaluations of basis functions in support at Delta-grid values are in 
                        ## the list 'basis' in this environment.
                        basisEnv = "environment",
                        
                        ## The 'basisPoints' contains the evaluation points for the basis functions.
                        basisPoints = "numeric",
                        
                        coefficients = "numeric",
                        fixedCoefficients = "list",
                        
                        ## The 'active' columns. Set in update, used in getModelMatrix and reset in computeModelMAtrix
                        modelMatrixCol = "numeric",
                        
                        ## The modelMatrix as a 'Matrix' is in this environment. Locked after computation.
                        modelMatrixEnv = "environment",

                        Omega = "matrix",
                        penalization = "logical",
                        var = "matrix",
                        
                        ## Which method is used to compute the estimate of the variance. 'pointProcessModel' has default 'Fisher'.
                        varMethod = "character"      
                        ),
         validity = function(object) {
           if(isTRUE(object@penalization) && !isTRUE(all.equal(min(eigen(object@Omega,only.values=TRUE,symmetric=TRUE)$values,0),0)))
             stop("Penalization matrix 'Omega' is not positive semi-definite.")
           if(isTRUE(object@support[2] - object@support[1] <= 0))
             stop("Variable 'support' has to be an interval.")
           if(isTRUE(object@Delta > object@support[2] - object@support[1]))
             stop("Variable 'Delta' has to be smaller than the length of the support.")
           return(TRUE)
         },
         contains="PointProcess"
         )


## TODO: Implement the stress-release type of model with
## intensity acccumulating according to time but "released"
## according to mark values.


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
                              modelMatrix = TRUE,
                              fit = modelMatrix,
                              varMethod = 'Fisher',
                              basisEnv,
                              ...) {
  
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
               delta = delta,
               formula = formula,
               family = family,
               call = match.call(),
               support = support,
               basisPoints = basisPoints,
               Delta = Delta,
               Omega = Omega,
               penalization = penalization,
               varMethod = varMethod)

  setProcessData(model) <- data

  if(missing(basisEnv))  {
    setBasis(model) <- list()
  } else {
    setBasis(model) <- basisEnv$basis
  }
  
  if(modelMatrix) {
    model <- computeModelMatrix(model)
  } else {
    setModelMatrix(model) <- Matrix()
  }

  parDim <- dim(getModelMatrix(model))[2]
  
  if(missing(coefficients)){
    coefficients <- rep(0,parDim)
  } else {
    if(length(coefficients) != parDim && parDim != 0) {
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

## TODO: Implement the use of interaction terms.

setMethod("computeBasis",c(model="PointProcessModel",form="terms"),
          function(model,form,...) {

            ## Basis evaluations for 'form' are computed if
            ## not already computed, locked and available in 'model@basisEnv$basis'.
            
            if(attr(form,"response")==1) form <- terms(formula(form)[-2])
            term <- attr(form,"term.labels")
            if(environmentIsLocked(model@basisEnv))
              {
                if(length(getBasis(model,term)) > 0)
                  {
                    return(TRUE)
                  }  else {
                    stop(paste("Term '",term,"' does not have precomputed basis values even though the basis environment is locked.",sep=""))
                  }
              }

            variables <- all.vars(form[[2]])
            form <- update(form,~.-1)
            
            if(length(variables) > 1){
              stop(paste("Basis computations with two or more variables in '",term,
                         "' is currently not supported.",sep=""))
            } else {
              x <- as.data.frame(model@basisPoints)
              colnames(x) <- variables
              model@basisEnv$basis[[term]] <- model.matrix(form,x)
              return(TRUE)
            }
            
            return(FALSE)
          }
          )

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

setMethod("computeLinearPredictor","PointProcessModel",
          function(model,coefficients=NULL,...){
            if(is.null(coefficients)) {
              coefficients <- coefficients(model)
            }                                    
            eta =  as.numeric(getModelMatrix(model) %*% coefficients)
            return(eta)
          }
          )

setMethod("computeModelMatrix","PointProcessModel",
          function(model,evaluationPositions=NULL,...){

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
            continuousProcess <- getContinuousProcess(getProcessData(model))
            DcontinuousVar <- paste("d.",colnames(getValue(continuousProcess)),sep="")

            id <- factor(getId(markedPointProcess))
            idLevels <- levels(id)
            
            marks <- getMarkType(markedPointProcess)
            markLevels <- levels(marks)
            
            ## The formula object is extracted and decomposed into terms.
            ## Each term label (in 'termLabels') is processed below,
            ## and the corresponding columns in the model matrix are computed.
            
            mt <- terms(model@formula)
            termLabels <- attr(mt,"term.labels")
            nrTerms <- length(termLabels)
            notFilterTerms <- character()
            
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
            
            for(i in 1:nrTerms) {
              term <- termLabels[i]

              variable <- all.vars(mt[i][[3]])
              
              ## The model matrix computed separately for terms
              ## involving marks, terms being linear filters of
              ## continuous processes and terms not being linear
              ## filters. The two former is done by a loop over for
              ## each value of 'id' and the result is stored in
              ## 'designList'.
              
              ## TODO: Can the C level computation return a sparse matrix
              ## directly?
              
              if(all(variable %in% markLevels))
                {
                  ## The call to computeBasis is invoked for its side effect
                  ## of computing the basis evaluations if that is not already
                  ## done (which is checked by checking if the environment 'basisEnv' 
                  ## is locked).
                  computeBasis(model,mt[i])
                  
                  assign <- c(assign,rep(i,dim(getBasis(model,term))[2]))
                  designList <- list()

                  ## Central loop over 'idLevels' and computations
                  ## of the model matrix in the C function
                  ## 'computePointProcessFilterMatrix'. Result is
                  ## converted to a sparse matrix, bound together in one
                  ## matrix below and stored in the list 'design'.
                  
                  for(i in idLevels) {
                    posi <- sort(positions[marks == variable & id == i])
                    ## Is this the right place to order the observed positions?
                    designList[[i]] <- Matrix(.Call("computePointProcessFilterMatrix",
                                                    evalPositions[[i]],
                                                    getBasis(model,term),
                                                    model@Delta,
                                                    posi,
                                                    PACKAGE="ppstat"),sparse=TRUE)
                  }
                  design[[term]] <- do.call("rBind",designList)
                  colnames(design[[term]]) <- colnames(getBasis(model,term))
                } else if(all(variable %in% DcontinuousVar)) {
                  computeBasis(model,mt[i])
                  assign <- c(assign,rep(i,dim(getBasis(model,term))[2]))
                  designList <- list()    

                  ## Central loop over 'idLevels' and computations of
                  ## the model matrix in the C function
                  ## 'computeContinuousProcessFilterMatrix'. Result is
                  ## converted to a sparse matrix, bound together in one
                  ## matrix below and stored in the list 'design'.
                  
                  for(i in idLevels) {
                    values <- getValue(continuousProcess)[id == i,variable== DcontinuousVar]
                    ## Is this the right place to order the observed positions?
                    designList[[i]] <- Matrix(.Call("computeContinuousProcessFilterMatrix",
                                                    evalPositions[[i]],
                                                    getBasis(model,term),
                                                    model@Delta,
                                                    values,
                                                    PACKAGE="ppstat"),sparse=TRUE)
                  }
                  design[[term]] <- do.call("rBind",designList)
                  colnames(design[[term]]) <- colnames(getBasis(model,term))
                } else {
                  ## The term does not involve filters
                  notFilterTerms <- c(notFilterTerms,term) 
                }
            }            

            ## Is there an intercept in the model? If so, add the intercept
            ## explicitly to the vector 'notMarkType'.
              
            if(attr(mt,"intercept") == 1){
              notFilterTerms <- paste(c(notFilterTerms,"1"),collapse="+")
              } else if(length(notFilterTerms) > 0) {
                notFilterTerms <-  paste(paste(notFilterTerms,collapse="+"),"-1")
              }
              
              ## Model matrix computations for terms involving 'id' and
              ## continuous time process components.
              
              if(length(notFilterTerms) > 0){
                form <-  as.formula(paste("~",notFilterTerms))
                variables <- all.vars(form)
                
                if(all(variables %in% c("id",colnames(getValue(continuousProcess))))) {
                  values <- data.frame(id=getId(continuousProcess))
                  notIdVariables <- variables[variables != "id"]
                  
                  if(length(notIdVariables) > 0) {
                    tmp <- as.matrix(getValue(continuousProcess)[,notIdVariables,drop=FALSE])
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
            
            
            
            modelMatrix <- cBind(X0,do.call("cBind",design))
            attr(modelMatrix,"assign") <- assign
            attr(modelMatrix,"formula") <- formula(model)
            setModelMatrix(model) <- modelMatrix
            lockEnvironment(model@modelMatrixEnv,binding=TRUE)
            lockEnvironment(model@basisEnv,binding=TRUE)            
            return(model)
          }
          )

setMethod("computeVar","PointProcessModel",
          function(model,...){
            if(attr(vcov(model),"method") == "none"){
              model@var <- matrix(0,length(coefficients(model)),length(coefficients(model)))
            } 
            if(attr(vcov(model),"method") == "Fisher") {
              vcovInv <- computeDDMinusLogLikelihood(model)
              if(model@penalization) vcovInv <- vcovInv + 2*model@Omega ## TODO: This requires some more thought ....
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
            }
            
            return(model)
          }
          )

setMethod("getLinearFilter",c(model="PointProcessModel",se="logical"),
          function(model,se=FALSE,nr=NULL...){
            filterTerms <- which(all.vars(formula(model)[[3]]) %in%
            c(levels(getMarkType(getMarkedPointProcess(getProcessData(model)))),
              paste("d.",colnames(getValue(getContinuousProcess(getProcessData(model)))),sep="")))
            
            mt <- terms(formula(model))

            linearFilter <- list()
            if(isTRUE(se)) linearFilterSE <- list()
            if(is.null(nr)) nr <- length(model@basisPoints)
            
            for(j in filterTerms){
              computeBasis(model,mt[j])
              term <- attr(mt[j],"term.labels")
              i <- seq_len(nr)*floor(dim(getBasis(model,term))[1]/nr)
              design <- getBasis(model,term)[i,,drop=FALSE]
              varName <- paste(all.vars(parse(text=term)),collapse=".")
              linearFilter[[varName]] <- design %*% coefficients(model)[dimnames(design)[[2]]]
              if(isTRUE(se))
                linearFilterSE[[varName]] <- sqrt(rowSums(design %*% vcov(model)[dimnames(design)[[2]],dimnames(design)[[2]]] * design))
            }
          
            lockEnvironment(model@basisEnv,bindings=TRUE)
            if(isTRUE(se)) {
              return(list(linearFilter=cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)),se=linearFilterSE))
            } else {
              return(cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)))
            }
              
          }
          )

setMethod("getBasis",c(model="PointProcessModel",term="ANY"),
          function(model,term,...){
            if(missing(term)) return(model@basisEnv$basis)
            return(model@basisEnv$basis[[term]])
          }
          )

setMethod("getModelMatrix",c(model="PointProcessModel",col="ANY"),
          function(model,col,...){
            if(missing(col)) col <- model@modelMatrixCol
            if(length(col) == 0) {
              return(model@modelMatrixEnv$modelMatrix)
            } else {
              modelMatrix <- model@modelMatrixEnv$modelMatrix[,col,drop=FALSE]
              attr(modelMatrix,"assign") <- attr(model@modelMatrixEnv$modelMatrix,"assign")[col]
              attr(modelMatrix,"formula") <- attr(model@modelMatrixEnv$modelMatrix,"formula")
              return(modelMatrix)
            }
          }
          )

setMethod("getModelMatrixEnv","PointProcessModel",
          function(model,...){
            return(list(modelMatrixEnv=model@modelMatrixEnv,modelMatrixCol=model@modelMatrixCol))
          }
          )

setReplaceMethod("setBasis",c(model="PointProcessModel", term="character", value="numeric"),
                 function(model,term,value){
                   if(environmentIsLocked(model@basisEnv))
                     model@basisEnv <- new.env(parent=.GlobalEnv)
                   
                   model@basisEnv$basis[[term]] <- value
                   
                   return(model)
                 }
                 )

setReplaceMethod("setBasis",c(model="PointProcessModel", term="ANY", value="list"),
                 function(model,term,value){
                   if(environmentIsLocked(model@basisEnv))
                     model@basisEnv <- new.env(parent=.GlobalEnv)
                   
                   model@basisEnv$basis <- value
                   
                   return(model)
                 }
                 )

setReplaceMethod("setModelMatrix",c(model="PointProcessModel",value="Matrix"),
          function(model,value){
            model@modelMatrixEnv <- new.env(parent=.GlobalEnv)
            model@modelMatrixEnv$modelMatrix <- value
            model@modelMatrixCol <- numeric()
            return(model)
          }
          )

setReplaceMethod("setModelMatrixEnv",c(model="PointProcessModel",value="list"),
                 function(model,value){
                   if(all(c("modelMatrixEnv","modelMatrixCol") %in% names(value))){
                     model@modelMatrixEnv <- value$modelMatrixEnv
                     model@modelMatrixCol <- value$modelMatrixCol
                   } else {
                     error("Right hand side of the assignment needs to be a list with two entries named 'modelMatrixEnv' and 'modelMatrixCol'")
                   }
                   return(model)
                 }
                 )

setMethod("predict","PointProcessModel",
          function(object,...) {
            eta <- computeLinearPredictor(object,...)
            return(object@family@phi(eta))
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

setMethod("ppmFit","PointProcessModel",
          function(model,control=list(),...) {
            nrPar <- parDim <- dim(getModelMatrix(model))[2]
            fixedPar <- model@fixedCoefficients
            Omega <- model@Omega

            ## If the model is a submodel we temporarily change the
            ## full model matrix to the submodel matrix in this function
            ## to avoid repeated subsetting of the full matrix. 
            
            if(length(model@modelMatrixCol) > 0) {
              modelMatrixEnv <- getModelMatrixEnv(model)
              setModelMatrix(model) <- getModelMatrix(model)
            }
              
            control <- c(list(maxit=1000),control)
            
            if(length(fixedPar) != 0) {
              nrPar <- parDim -length(fixedPar$which)
              tmpPar <- numeric(parDim)
              tmpPar[fixedPar$which] <- fixedPar$value
            }
            
            ## Setting up the initial parameters
            
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
            

            ## Setting up the objective function to minimize

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

            ## The actual minimization
            
            model@optimResult <- optim(initPar,mll,gr=dmll,method="BFGS",control=control,...)

            if(length(fixedPar) == 0) {
              model@coefficients <- model@optimResult$par
            } else {
              model@coefficients[-fixedPar$which] <- model@optimResult$par
              model@coefficients[fixedPar$which] <- fixedPar$value
            }
            names(model@coefficients) <- dimnames(getModelMatrix(model))[[2]]

            ## Resetting the full model matrix if required

            if(exists("modelMatrixEnv")) {
              setModelMatrixEnv(model) <- modelMatrixEnv
            }
            
            ## Computation of the estimated covariance matrix
            
            model <- computeVar(model)
            return(model)
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

setMethod("subset","PointProcessModel",
          function(x,...){
            pointProcessModel(formula = formula(x),
                              data = subset(getProcessData(x),...),
                              family = x@family,
                              Delta = x@Delta,
                              support = x@support,
                              basisPoints = x@basisPoints,
                              Omega = x@Omega,
                              coefficients = coefficients(x),
                              fixedCoefficients = x@fixedCoefficients,
                              basisEnv = x@basisEnv
                              )
          }
          )

### TODO: Summary should return an S4 object instead with appropriate view method.
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

setMethod("vcov","PointProcessModel",
          function(object,...){
            attr(object@var,"method") <- object@varMethod
            return(object@var)
          }
          )

setMethod("update","PointProcessModel",
          function(object,...){
            .local <- function(model,formula = .~.,warmStart = TRUE,fixedCoefficients = list(),...){

              modelFormula <- formula(model)
              updatedFormula <- update(modelFormula,formula)
              updatedTermLabels <- attr(terms(updatedFormula),"term.labels")
              superTermLabels <- attr(terms(attr(getModelMatrix(model,numeric()),"formula")),"term.labels")
              
              if(attr(terms(updatedFormula),"intercept")==1) updatedTermLabels <- c("Intercept",updatedTermLabels)
              if(attr(terms(attr(getModelMatrix(model,numeric()),"formula")),"intercept")==1) superTermLabels <- c("Intercept",superTermLabels)

              superTermLabels <- superTermLabels[unique(attr(getModelMatrix(model,numeric()),"assign"))+1]
              
              formula(model) <- updatedFormula
              
              if(length(model@modelMatrixCol) == 0){
                tmpCoef <- coefficients(model)
              } else {
                tmpCoef <- rep(0,dim(getModelMatrix(model,numeric()))[2])
                tmpCoef[model@modelMatrixCol] <- coefficients(model)
              }

              if(all(updatedTermLabels %in% superTermLabels)) {
                if("Intercept" %in%  superTermLabels) {
                  col <- which(superTermLabels[attr(getModelMatrix(model,numeric()),"assign")+1] %in% updatedTermLabels)
                } else {
                  col <- which(superTermLabels[attr(getModelMatrix(model,numeric()),"assign")] %in% updatedTermLabels)
                }
                ## TODO: The following checks whether there are interactions in the model. There is a general problem with
                ## choices of contrasts in submodels that should be dealt with. Removing e.g. the intercept gives a peculiar
                ## model it seems ... 
                if(any(attr(terms(modelFormula),"order") >= 2)) {
                  warning(paste(c(object@call,"Original model formula includes interaction terms. Check that updated model is as expected."),collapse="\n"),call.=FALSE)
                }
                if(length(col) == length(tmpCoef)) {
                  model@modelMatrixCol <- numeric()
                } else {
                  model@modelMatrixCol <- col
                }
                if(warmStart) {
                  coefficients(model) <- tmpCoef[col]
                } else {
                  coefficients(model) <- rep(0,length(col))
                }
              } else {
                model <- computeModelMatrix(model)
                coefficients(model) <- rep(0,dim(getModelMatrix(model)[2]))
              }
              
              if(length(fixedCoefficients) != 0) model@coefficients[fixedCoefficients$which] <- fixedCoefficients$value
              model@fixedCoefficients <- fixedCoefficients
                          
              return(ppmFit(model,...))
            }
            object@call <- match.call()
            .local(object,...)
          }
          )
