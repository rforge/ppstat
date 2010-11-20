pointProcessModel <- function(
                              formula,
                              data,
                              family,
                              support,
                              N = 200,
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
      if(!(missing(support)))
        {
          if(length(support) == 1)
            support <- c(0,max(support[1],0))
          
          if(missing(Delta))
            Delta <- (support[2] - support[1])/N
          
          basisPoints <- sort(unique(c(0, seq(support[1], support[2], Delta))))
        } else {
          stop("Must specify either 'support' or 'basisPoints'.")
        }
    } else {
      basisPoints <- sort(unique(c(0,basisPoints)))
      support = range(basisPoints)
      Delta = min(diff(basisPoints))
    }
      
  delta <- as.numeric(unlist(tapply(getPosition(data),
                                    getId(data),
                                    function(x) c(diff(x),0)), use.names=FALSE))

  model <- new("PointProcessModel",
               delta = delta,
               family = family,
               formula = .~.,
               call = match.call(),
               processData = data,
               support = support,
               basisPoints = basisPoints,
               Delta = Delta,
               Omega = Omega,
               penalization = penalization,
               varMethod = varMethod)

  formula(model) <- formula
  
  if(!anticipating(model)) {
    model@basisPoints <- model@basisPoints[model@basisPoints >= 0]
    model@support[1] <- max(0,model@support[1])
  }
    
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
    coefficients <- rep(.Machine$double.eps, parDim)
    names(coefficients) <- dimnames(getModelMatrix(model))[[2]]
  } else {
    if(length(coefficients) != parDim && parDim != 0) {
      coefficients <- rep(.Machine$double.eps, parDim)
      names(coefficients) <- dimnames(getModelMatrix(model))[[2]]
      warning(paste("Incorrect length of initial parameter vector. Initial parameters all set to", .Machine$double.eps))
    }
  }


  if(length(fixedCoefficients) != 0)
    coefficients[fixedCoefficients$which] <- fixedCoefficients$value

  model@fixedCoefficients <- fixedCoefficients
  coefficients(model) <- coefficients
            
  if(fit)
    model <- ppmFit(model,...)
  

  return(model)
}

## TODO: Implement the use of interaction terms.

setMethod("computeBasis", c(model = "PointProcessModel", form = "ANY"),
          function(model, form, ...) {
            
            if(class(form)[1] != "terms")
              stop("The 'form' argument must be of S3-class 'terms'")
            
            ## Basis evaluations for 'form' are computed if
            ## not already computed, locked and available in 'model@basisEnv$basis'.
            
            if(attr(form, "response") == 1) form <- delete.response(form)
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
          function(model, coefficients = NULL, ...){
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            eta <- computeLinearPredictor(model, coefficients, ...)
            
            if(model@family@link == "log") {

              dmll <- as.vector(t(exp(eta)*model@delta)%*%getModelMatrix(model)) -
                colSums(getModelMatrix(model)[getPointPointer(processData(model), response(model)), , drop = FALSE])

            } else {
              
              etaP <- eta[getPointPointer(processData(model), response(model))]
              mmP <- getModelMatrix(model)[getPointPointer(processData(model), response(model)), , drop = FALSE]

              dmll <-  as.vector(t(model@family@Dphi(eta)*model@delta) %*% getModelMatrix(model)) -
                as.vector(t(model@family@Dphi(etaP)/model@family@phi(etaP)) %*% mmP)
              
            }
            
            return(dmll)
            
          }
          )

setMethod("computeDDMinusLogLikelihood", "PointProcessModel",
          function(model, coefficients = NULL, ...){
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")

            eta <- computeLinearPredictor(model, coefficients, ...)
            
            if(model@family@link == "log"){
              
              ddmll <-  as(crossprod(getModelMatrix(model), exp(eta)*model@delta*getModelMatrix(model)), "matrix")
              
            } else if(model@family@link == "identity"){
              
              etaP <- eta[getPointPointer(processData(model), response(model))]
              mmP <- getModelMatrix(model)[getPointPointer(processData(model), response(model)), , drop = FALSE]
              
              ddmll <-  as(crossprod(mmP, 1/model@family@phi(etaP)^2*mmP), "matrix")
              
            } else {
              
              etaP <- eta[getPointPointer(processData(model), response(model))]
              mmP <- getModelMatrix(model)[getPointPointer(processData(model), response(model)), , drop = FALSE]
              
              ddmll <-  as(crossprod(getModelMatrix(model), model@family@D2phi(eta)*model@delta*getModelMatrix(model)), "matrix") -
                as(crossprod(mmP, (model@family@D2phi(etaP)*model@family@phi(etaP) - model@family@Dphi(etaP)^2)/model@family@phi(etaP)^2*mmP), "matrix")
              
            }
            
            return(ddmll)
            
          }
          )


setMethod("computeWeights", "PointProcessModel",
          function(model, coefficients = NULL, ...) {
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            
            eta <- computeLinearPredictor(model, coefficients, ...)
          

            if(model@family@link == "log"){

              w <-  exp(eta)*model@delta

             } else if(model@family@link == "identity"){

               points <- getPointPointer(processData(model), response(model))
               etaP <- eta[points]
               w <- rep(0, length(eta))
               w[points] <- 1/model@family@phi(etaP)^2

             } else {

               points <- getPointPointer(processData(model), response(model))
               etaP <- eta[points]
               w <- model@family@D2phi(eta)*model@delta
               w[points] <- w[points] - (model@family@D2phi(etaP)*model@family@phi(etaP) - model@family@Dphi(etaP)^2)/model@family@phi(etaP)^2
             
             }
            
            return(w)
          }
          )

setMethod("computeWorkingResponse", "PointProcessModel",
          function(model, coefficients = NULL, ...){
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            eta <- computeLinearPredictor(model, coefficients, ...)
            points <- getPointPointer(processData(model), response(model))
            w <- computeWeights(model)
            w[w < .Machine$double.eps] <- 1
            
            if(model@family@link == "log") {

              wr <- exp(eta)*model@delta
              wr[points] <- wr[points] - 1
              wr <- eta - wr/w
              
            } else {
              
              etaP <- eta[points]

              wr <-  model@family@Dphi(eta)*model@delta
              wr[points] <- wr[points] - model@family@Dphi(etaP)/model@family@phi(etaP)
              wr <- eta - wr/w
              
            }
            
            return(wr)
            
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
          function(model, evaluationPositions = NULL, ...){

            ## The 'model' of class PointProcessModel contains the data
            ## as an object of class MarkedPointProcess and the formula for the
            ## model specification. The 'evalPositions' below corresponding to
            ## the model matrix rows are either given by the 'evaluationPositions'
            ## argument or extracted from the the MarkedPointProcess object (default).

            if(is.null(evaluationPositions)) {
              evalPositions <- tapply(getPosition(processData(model)),
                                      getId(processData(model)), list)
            } else {
              evalPositions <- evaluationPositions
            }

            ## Checks if the model is allowed to be anticipating and sets
            ## the 'zero' accordingly.

            if(anticipating(model)) {
              zero <- which(model@basisPoints == 0)
            } else {
              zero <- 0
            }
              
            
            ## The observed points ('positions') for the marked point process,
            ## the corresponding 'id' labels and 'marks' are extracted.
            
            processData <- processData(model)
            positions <- getPointPosition(processData)
            DcontinuousVar <- paste(colnames(getValue(processData)),".d",sep="")

            id <- factor(getPointId(processData))
            idLevels <- levels(id)
            
            marks <- getMarkType(processData)
            markLevels <- levels(marks)
            
            ## The formula object is extracted and decomposed into terms.
            ## Each term label (in 'termLabels') is processed below,
            ## and the corresponding columns in the model matrix are computed.
            
            mt <- delete.response(terms(formula(model)))
            termLabels <- attr(mt,"term.labels")
            nrTerms <- length(termLabels)
            notFilterTerms <- numeric()
            
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

              variable <- all.vars(mt[i])
              
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
                  computeBasis(model, mt[i])
                  
                  assign <- c(assign, rep(i, dim(getBasis(model,term))[2]))
                  designList <- list()

                  ## Central loop over 'idLevels' and computations of
                  ## the model matrix in the C function
                  ## 'computeFilterMatrix'. Result is converted to a
                  ## sparse matrix, bound together in one matrix below
                  ## and stored in the list 'design'.
                  
                  for(i in idLevels) {
                    posi <- positions[marks == variable & id == i]
                    ## posi is sorted for a valid data object. This is
                    ## assumed in the following computation.
                    designList[[i]] <- Matrix(.Call(computeFilterMatrix,
                                                    evalPositions[[i]],
                                                    getBasis(model, term),
                                                    model@Delta,
                                                    posi,
                                                    zero,
                                                    'p'), sparse=TRUE)
                  }
                  design[[term]] <- do.call("rBind",designList)
                  colnames(design[[term]]) <- colnames(getBasis(model,term))
                } else if(all(variable %in% c(DcontinuousVar,"d.position","d.time"))) {
                  computeBasis(model, mt[i])
                  assign <- c(assign, rep(i, dim(getBasis(model,term))[2]))
                  designList <- list()    

                  ## Central loop over 'idLevels' and computations of
                  ## the model matrix in the C function
                  ## 'computeFilterMatrix'. Result is
                  ## converted to a sparse matrix, bound together in one
                  ## matrix below and stored in the list 'design'.
                  if(variable %in% c("d.position","d.time")) {
                    values <- getPosition(processData)
                  } else {
                    values <- getValue(processData)[ , variable == DcontinuousVar]
                  }
                  
                  for(i in idLevels) {
                    valuesi <- values[getId(processData) == i]
                    
                    designList[[i]] <- Matrix(.Call(computeFilterMatrix,
                                                    evalPositions[[i]],
                                                    getBasis(model,term),
                                                    model@Delta,
                                                    valuesi,
                                                    zero,
                                                    'c'), sparse=TRUE)
                  }
                  design[[term]] <- do.call("rBind", designList)
                  colnames(design[[term]]) <- colnames(getBasis(model,term))
                } else {
                  ## The term does not involve filters
                  notFilterTerms <- c(notFilterTerms,i) 
                }
            }            

            ## Model matrix computations for terms involving 'id',
            ## 'position/time', unit variables and non-filtered
            ## continuous time process components.
              
            if(length(notFilterTerms) > 0 || attr(mt, "intercept") == 1){
              form <-  mt[notFilterTerms]
              attr(form, "intercept") <-  attr(mt, "intercept")
              variables <- all.vars(form)
              
              if(all(variables %in% c(processData@idVar, processData@positionVar, colnames(getValue(processData)), colnames(getUnitData(processData))))) {
                ## TODO: Investigate if this build of the model matrix can exploit
                ## sparse matrices more directly, if it is necessary to create the
                ## values object, or if we could refer to a suitable environment?
                
                values <- list()

                if(processData@idVar %in% variables || attr(mt, "intercept") == 1) {
                  values[[1]] <- data.frame(getId(processData))
                  names(values[[1]]) <- processData@idVar
                }
                
                if(processData@positionVar %in% variables) {
                  values[[2]] <- data.frame(getPosition(processData))
                  names(values[[2]]) <-  processData@positionVar
                }

                unitVariables <- colnames(getUnitData(processData)) %in% variables
                if(any(unitVariables)) {
                  values[[3]] <- getUnitData(processData)[as.numeric(getId(processData)),
                                                          unitVariables,
                                                          drop = FALSE]
                }
                
                otherVariables <-  colnames(getValue(processData)) %in% variables 
                if(any(otherVariables)) {
                  values[[4]] <- as.data.frame(as.matrix(getValue(processData)[ , otherVariables, drop=FALSE]))
                  rownames(values[[4]]) <- NULL
                }
                
                values <- do.call("cbind", values[!sapply(values, is.null)])
                ## This sparse.model.matrix solution avoids the dense
                ## model matrix, but does not produce the assign
                ## attribute ...
                ## X0 <- sparse.model.matrix(form, values)

                tmp <- model.matrix(form, values)
                assign <- c(c(0,notFilterTerms)[attr(tmp,"assign")+1],assign)
                X0 <- Matrix(tmp, dimnames = dimnames(tmp), sparse=TRUE)
                assign <- c(c(0, notFilterTerms)[attr(X0, "assign") + 1], assign)
              } else {
                stop(paste("Use of non existing variable(s) in:", form))
              }
            }
            
            
            if(exists("X0")) {
              modelMatrix <- cBind(X0, do.call("cBind",design))
            } else {
              modelMatrix <- do.call("cBind",design)
            }
            attr(modelMatrix, "assign") <- assign
            form <- formula(model)
            attr(form, "filterTerms") <- which(!(seq(along=termLabels) %in% notFilterTerms))
            formula(model) <- form
            attr(modelMatrix, "formula") <- form
            setModelMatrix(model) <- modelMatrix
            lockEnvironment(model@modelMatrixEnv, binding=TRUE)
            lockEnvironment(model@basisEnv, binding=TRUE)            
            return(model)
          }
          )

setMethod("computeVar", "PointProcessModel",
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

setMethod("getLinearFilter", c(model = "PointProcessModel"),
          function(model, se = FALSE, nr, ...){
            mt <- delete.response(terms(formula(model)))
            ## TODO: filterTerms attribute first computed in compute model
            ## matrix. Should be computed globally prior to this!
            filterTerms <- attr(formula(model), "filterTerms")
            
            linearFilter <- list()
            design <- list()
            if(isTRUE(se))
              linearFilterSE <- list()
            
            if(missing(nr))
              nr <- length(model@basisPoints)

            for(j in filterTerms){
              computeBasis(model, mt[j])
              term <- attr(mt[j], "term.labels")
              NR <- dim(getBasis(model, term))[1]
              i <- seq_len(min(nr,NR))*floor(max(NR/nr,1))
              varName <- paste(all.vars(parse(text = term)), collapse = ".")
              design[[varName]] <- cbind(design[[varName]], getBasis(model, term)[i, ,drop=FALSE])}

            for(j in seq(along=design)){
              linearFilter[[j]] <- design[[j]] %*% coefficients(model)[dimnames(design[[j]])[[2]]]
              if(isTRUE(se))
                linearFilterSE[[j]] <- sqrt(rowSums(design[[j]] %*% vcov(model)[dimnames(design[[j]])[[2]], dimnames(design[[j]])[[2]]] * design[[j]]))
            }

            names(linearFilter) <- names(design)
          
            lockEnvironment(model@basisEnv, bindings=TRUE)
            if(isTRUE(se)) {
              return(list(linearFilter = cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)), se = linearFilterSE))
            } else {
              return(cbind(data.frame(x = model@basisPoints[i]), as.data.frame(linearFilter)))
            }
              
          }
          )

setMethod("getBasis", c(model = "PointProcessModel", term = "ANY"),
          function(model, term,...){
            if(missing(term)) return(model@basisEnv$basis)
            return(model@basisEnv$basis[[term]])
          }
          )

setMethod("getModelMatrix", c(model = "PointProcessModel", col = "ANY"),
          function(model, col,...){
            if(missing(col)) col <- model@modelMatrixCol
            if(length(col) == 0) {
              return(model@modelMatrixEnv$modelMatrix)
            } else {
              modelMatrix <- model@modelMatrixEnv$modelMatrix[ , col, drop=FALSE]
              attr(modelMatrix, "assign") <- attr(model@modelMatrixEnv$modelMatrix, "assign")[col]
              attr(modelMatrix, "formula") <- attr(model@modelMatrixEnv$modelMatrix, "formula")
              return(modelMatrix)
            }
          }
          )

setMethod("getModelMatrixEnv", "PointProcessModel",
          function(model,...){
            return(list(modelMatrixEnv = model@modelMatrixEnv, modelMatrixCol = model@modelMatrixCol))
          }
          )

setReplaceMethod("setBasis", c(model = "PointProcessModel", term = "character", value = "numeric"),
                 function(model, term, value){
                   if(environmentIsLocked(model@basisEnv))
                     model@basisEnv <- new.env(parent=.GlobalEnv)
                   
                   model@basisEnv$basis[[term]] <- value
                   
                   return(model)
                 }
                 )

setReplaceMethod("setBasis", c(model = "PointProcessModel", term = "ANY", value = "list"),
                 function(model, term, value){
                   if(environmentIsLocked(model@basisEnv))
                     model@basisEnv <- new.env(parent=.GlobalEnv)
                   
                   model@basisEnv$basis <- value
                   
                   return(model)
                 }
                 )

setReplaceMethod("setModelMatrix", c(model = "PointProcessModel", value = "Matrix"),
          function(model, value){
            model@modelMatrixEnv <- new.env(parent=.GlobalEnv)
            model@modelMatrixEnv$modelMatrix <- value
            model@modelMatrixCol <- numeric()
            return(model)
          }
          )

setReplaceMethod("setModelMatrixEnv", c(model = "PointProcessModel", value = "list"),
                 function(model,value){
                   if(all(c("modelMatrixEnv","modelMatrixCol") %in% names(value))){
                     model@modelMatrixEnv <- value$modelMatrixEnv
                     model@modelMatrixCol <- value$modelMatrixCol
                   } else {
                     stop("Right hand side of the assignment needs to be a list with two entries named 'modelMatrixEnv' and 'modelMatrixCol'")
                   }
                   return(model)
                 }
                 )

setMethod("predict","PointProcessModel",
          function(object, ...) {
            eta <- computeLinearPredictor(object,...)
            return(object@family@phi(eta))
          }
          )

setMethod("termPlot","PointProcessModel",
          function(model, alpha = 0.05, layer = geom_line(), trans = NULL, ...) {
            if(alpha <= 0 || alpha > 1)
              stop("The 'alpha' level must be in (0,1]")

            if(length(attr(formula(model), "filterTerms")) == 0){
              print("No filter function terms to plot")
              return(invisible())
            }
            
            if(alpha == 1) {
              se <- FALSE
            } else {
              se <- TRUE
              q <- qnorm(1-alpha/2)
            }
            
            linearFilter <- getLinearFilter(model, se = se, nr = 400)
            if(se) {
              moltenFilter <- melt(linearFilter$linearFilter, id.vars = "x")
              plotData <- cbind(moltenFilter,
                                data.frame(cf.lower = moltenFilter$value-q*unlist(linearFilter$se),
                                           cf.upper = moltenFilter$value+q*unlist(linearFilter$se)))
              if(!is.null(trans))
                plotData[, c("value", "cf.lower", "cf.upper")] <- do.call(trans,list(plotData[, c("value", "cf.lower", "cf.upper")]))
            } else {
              plotData <- melt(linearFilter, id.vars = "x")
               if(!is.null(trans))
                 plotData$value <- do.call(trans, plotData$value)
            }
           

            linearFilterPlot <- ggplot(data = plotData, aes(x = x, y = value)) +
              facet_grid(variable ~ ., scales = "free_y") +
                scale_x_continuous("position") +
                  scale_y_continuous("") + layer

            if(se)
              linearFilterPlot <- linearFilterPlot + geom_ribbon(aes(min = cf.lower, max = cf.upper), fill = alpha("blue", 0.2))

            return(linearFilterPlot)
          }
          )

setMethod("ppmFit", "PointProcessModel",
          function(model, control = list() , ...) {
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
            
            if(!("maxit" %in% names(control)))
              control <- c(list(maxit = 1000), control)
            
            
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
              initPar <- rep(.Machine$double.eps,nrPar)             
            warning("Length of initial parameter vector worng. Initial parameters all set to 0.")
          }
            

            ## Setting up the objective function to minimize

            if(!model@penalization){
              if(length(fixedPar) == 0) {
                mll <- function(par,...) computeMinusLogLikelihood(model, par, ...)
                dmll <- function(par,...) computeDMinusLogLikelihood(model, par, ...)
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

            args <- list(...)
            args[["par"]] <- initPar
            args[["fn"]] <- mll
            args[["gr"]] <- dmll
            args[["control"]] <- control

            if(family(model)@link == "identity") {
              method <- "L-BFGS-B"
              if(attr(terms(formula(model)),"intercept")==1) {
                lower = c(sqrt(.Machine$double.eps),rep(0,dim(getModelMatrix(model))[2]-1))
              } else {
                lower = sqrt(.Machine$double.eps)
              }
            } else {
              method <- "BFGS"
              lower <- -Inf
            }
            
            if(!("method" %in% names(args))) args[["method"]] <- method
            if(!("lower" %in% names(args))) args[["lower"]] <- lower
                        
            model@optimResult <- do.call("optim", args)
            
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

setMethod("IWLS", "PointProcessModel",
          function(model, control, ...) {

            value <- computeMinusLogLikelihood(model)
            reltol <- sqrt(.Machine$double.eps)
            maxit <- 100
            trace <- 1
            
            i <- 1
            while(i < maxit) {
              w <- computeWeights(model)
              z <- computeWorkingResponse(model)
              ## TODO: check when lm.fit.sparse becomes exported.
              ## This is relying on an algorithm in the MatrixModels package still
              ## under development. 
              coefficients(model) <- MatrixModels:::lm.fit.sparse(getModelMatrix(model), z, w)
              val <- computeMinusLogLikelihood(model)
              i <- i + 1
              if(trace > 0)
                cat("Value: ", val, "\n")
              if(value < reltol * (abs(value) + reltol) * value + val)
                break
              value <- val
            }
            
            return(model)
          }
          )

          

setMethod("print","PointProcessModel",
          function(x, digits = max(3, getOption("digits") - 3), ...){
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

setMethod("show", "PointProcessModel",
          function(object) print(x=object)
          )

setMethod("subset", "PointProcessModel",
          function(x, ...){
            pointProcessModel(formula = formula(x),
                              data = subset(processData(x), ...),
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

setMethod("summary", "PointProcessModel",
          function(object,...) {
            result <- list()

            result$df <- length(object@coefficients) - length(object@fixedCoefficients$which)
            result$call <- object@call
            result$mll <- object@optimResult$value
            result$iter <- object@optimResult$counts
            result$convergence <- object@optimResult$convergence
            result$aic <- getInformation(object)

            ## TODO: Implement some way of summarizing 'residuals'

            se <- sqrt(diag(object@var))
            z <- object@coefficients/se
            q <- 2*(1-pnorm(abs(z)))
            z[se==0] <- NA
            q[se==0] <- NA
            se[se==0] <- NA
            
            result$coefficients <- matrix(c(object@coefficients, se, z, q), ncol=4)
            zapNames <- sapply(names(object@coefficients),
                               function(y) {
                                 x <- strsplit(y, "")[[1]]
                                 if(length(x) > 18){
                                   end <- c("...", x[(length(x)-2):length(x)])
                                   return(paste(c(x[1:12], end), sep = "", collapse = ""))
                                 }
                                   return(y)
                               })
                             
            
            dimnames(result$coefficients) <- list(zapNames,c("Estimate","Std. Error", "z value", "Pr(> |z|)"))

            result$fixedCoefficients <- object@fixedCoefficients
            result$penalization <- object@penalization
            class(result) <- c("summary.ppm")
            return(result)

          }
          )

setMethod("vcov", "PointProcessModel",
          function(object,...){
            attr(object@var,"method") <- object@varMethod
            return(object@var)
          }
          )

setMethod("update", "PointProcessModel",
          function(object, formula = .~., warmStart = TRUE, fixedCoefficients = list(), ...){

            
            modelFormula <- formula(object)
            updatedFormula <- update(modelFormula, formula)
            updatedTermLabels <- attr(terms(updatedFormula), "term.labels")
            superTermLabels <- attr(terms(attr(getModelMatrix(object, numeric()), "formula")), "term.labels")
            attr(updatedFormula, "filterTerms") <- which(updatedTermLabels %in% attr(terms(modelFormula), "term.labels")[attr(modelFormula, "filterTerms")])
            
            if(attr(terms(updatedFormula), "intercept") == 1) updatedTermLabels <- c("Intercept", updatedTermLabels)
            if(attr(terms(attr(getModelMatrix(object,numeric()), "formula")), "intercept") == 1) superTermLabels <- c("Intercept", superTermLabels)
           
            formula(object) <- updatedFormula
            
            if(length(object@modelMatrixCol) == 0){
              tmpCoef <- coefficients(object)
            } else {
              tmpCoef <- rep(.Machine$double.eps ,dim(getModelMatrix(object, numeric()))[2])
              tmpCoef[object@modelMatrixCol] <- coefficients(object)
            }

            if(all(updatedTermLabels %in% superTermLabels)) {
              if("Intercept" %in%  superTermLabels) {
                col <- which(superTermLabels[attr(getModelMatrix(object,numeric()),"assign")+1] %in% updatedTermLabels)
              } else {
                col <- which(superTermLabels[attr(getModelMatrix(object,numeric()),"assign")] %in% updatedTermLabels)
              }
              ## TODO: The following checks whether there are interactions in the model. There is a general problem with
              ## choices of contrasts in submodels that should be dealt with. Removing e.g. the intercept gives a peculiar
              ## model it seems ... 
              if(any(attr(terms(modelFormula), "order") >= 2)) {
                warning(paste(c(object@call,"Original model formula includes interaction terms. Check that updated model is as expected."),collapse="\n"), call. = FALSE)
              }
              if(length(col) == length(tmpCoef)) {
                object@modelMatrixCol <- numeric()
              } else {
                object@modelMatrixCol <- col
              }
              if(warmStart) {
                coefficients(object) <- tmpCoef[col]
              } else {
                coefficients(object) <- rep(0,length(col))
              }
            } else {
              object <- computeModelMatrix(object)
              coefficients(object) <- rep(.Machine$double.eps, dim(getModelMatrix(object))[2])
            }
            
            if(length(fixedCoefficients) != 0) {
              object@coefficients[fixedCoefficients$which] <- fixedCoefficients$value
              object@fixedCoefficients <- fixedCoefficients
            } else if(length(object@fixedCoefficients) != 0 && exists("col", inherits = FALSE)) {
              fixedCoefficients <-  list() 
              fixedCoefficients$which <- which(col %in% object@fixedCoefficients$which)
              fixedCoefficients$value <- object@fixedCoefficients$value[object@fixedCoefficients$which %in% col]
              object@fixedCoefficients <- fixedCoefficients
              object@coefficients[fixedCoefficients$which] <- fixedCoefficients$value
            } else {
              object@fixedCoefficients <- list()
            }   
            

            call <- as.list(object@call)
            call$formula <- formula(object)
            object@call <- as.call(call)
            
            return(ppmFit(object, ...))
          }
          )

setMethod("getInformation", "PointProcessModel",
          function(model, k = 2, ...) {
            df <- length(model@coefficients) - length(model@fixedCoefficients$which)
            mll <- computeMinusLogLikelihood(model)
            return(2*mll + k*df)
          }
          )

setMethod("stepInformation", "PointProcessModel",
          function(model, direction = "both", trace = 1, steps = 1000, warmStart = TRUE, k = 2, ...) {
            if(trace == 1) 
              cat(" Step \t  AIC \t\t Direction\n")
            
            AIC <- list()
            models <- list()
            models[["current"]] <- model
            models$current@varMethod <- 'none'
            step <- 1
            initialForm <- terms(formula(model))
            initialTerms <- attr(initialForm, "term.labels")
            currentTerms <- initialTerms
            termsToAdd <- character()
            termsToRemove <- currentTerms
            minDirection <- "start"
            while(step < steps) {
              AIC$current <- getInformation(models$current, k, ...)

              if(trace == 1)
                cat(" ", step, "\t", AIC$current, "\t ", minDirection, "\n")

              if(trace == 2) {
                cat("Step:", step, "\tAIC:", AIC$current, "\n")
                cat("Model:", deparse(formula(models$current)), "\n")
                cat("Direction:", minDirection, "\n\n")
                
              }

              if(direction == "both") {
                if(length(termsToAdd) > 0) {
                  AIC$forward <- numeric(length(termsToAdd))
                  names(AIC$forward) <- termsToAdd
                  models$forward <- list()
                  for(term in termsToAdd) {
                    models$forward[[term]] <- update(models$current,
                                                     as.formula(paste(".~. +", term)),
                                                     warmStart = warmStart)
                    AIC$forward[term] <- getInformation(models$forward[[term]], k, ...)
                  }
                }
              }
              
              if(direction %in% c("both", "backward")) {
                if(length(termsToRemove) > 1 || (length(termsToRemove) == 1 & (attr(initialForm, "intercept") == 1))) {
                  AIC$backward <- numeric(length(termsToRemove))
                  names(AIC$backward) <- termsToRemove
                  models$backward <- list()
                  for(term in termsToRemove) {
                    models$backward[[term]] <- update(models$current,
                                                      as.formula(paste(".~. -", term)),
                                                      warmStart = warmStart)
                    AIC$backward[term] <- getInformation(models$backward[[term]], k, ...)
                  }
                }
              }

              minDirection <- names(which.min(sapply(AIC, min)))

              if(minDirection == "current")              
                break()
              
              if(minDirection == "forward") {
                termIndex <- which.min(AIC$forward)
                models$current <- models$forward[[termIndex]]
                termsToRemove <- currentTerms
                currentTerms <- c(currentTerms, termsToAdd[termIndex]) 
                termsToAdd <- termsToAdd[-termIndex]
              }
              
              if(minDirection == "backward") {
                termIndex <- which.min(AIC$backward)
                models$current <-  models$backward[[termIndex]]
                termsToAdd <- initialTerms[!(initialTerms %in% currentTerms)]
                currentTerms <- termsToRemove[-termIndex]
                termsToRemove <- currentTerms
              }
              
              step <- step + 1
            }
            models$current@varMethod <- model@varMethod
            model <- ppmFit(models$current)
            return(invisible(model))
          }
          )
