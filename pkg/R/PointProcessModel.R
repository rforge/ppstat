pointProcessModel <- function(
                              formula,
                              data,
                              family,
                              support = 1,
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
      if(length(support) == 1)
        support <- c(0, max(support[1], 0))
          
      if(missing(Delta))
        Delta <- (support[2] - support[1])/N
          
      basisPoints <- sort(unique(c(0, seq(support[1], support[2], Delta))))
    } else {
      basisPoints <- sort(unique(c(0, basisPoints)))
      support = range(basisPoints)
      Delta = min(diff(basisPoints))
    }
      
  delta <- as.numeric(unlist(tapply(getPosition(data),
                                    getId(data),
                                    function(x) c(0, diff(x))),
                             use.names=FALSE)
                      )
                      
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

  ## Parsing the formula: Either the model is multivariate with the
  ## formula a list of formulae, one for each model, or the response
  ## is a vector -- assuming that the right hand side of each model is 
  ## the same as a starting point. Alternatively, the formula is just a
  ## single formula. 
  
  if(is.list(formula)) {
    ## Extract all terms used in any of the formulae to create a
    ## super formula.
    superFormula <- reformulate(unlist(lapply(formula, function(f) attr(terms(f), "term.labels"))))
    formula(model) <- superFormula
    response <- NULL
  } else {
    response <- attr(terms(formula), "variables")[[2]]
    if(length(response) > 2 && response[[1]] == as.symbol("c")) {
      response <- response[seq.int(2, length(response))]
      formula(model) <- update(formula, as.formula(paste(deparse(response[[1]]), "~ .")))
      call <- as.list(model@call)
      call$formula <- formula(model)
      model@call <- as.call(call)
    } else {
      formula(model) <- formula
      response <- NULL
    }
  }
  
  if(!anticipating(model)) {
    model@basisPoints <- model@basisPoints[model@basisPoints >= 0]
    model@support[1] <- max(0, model@support[1])
  }
  
  if(missing(basisEnv))  {
    setBasis(model) <- list()
    model <- computeBasis(model)
  } else {
    setBasis(model) <- basisEnv$basis
  }

  if(modelMatrix) {
    model <- computeModelMatrix(model)
  } else {
    model <- updateModelMatrix(model,
                               Matrix(),
                               assign = numeric(),
                               form = formula(~0)
                               )
  }
  
  if(missing(coefficients))
    coefficients <- .Machine$double.eps
  
  if(length(fixedCoefficients) != 0)
    coefficients[fixedCoefficients$which] <- fixedCoefficients$value
  
  model@fixedCoefficients <- fixedCoefficients
  coefficients(model) <- coefficients
  
  if(!is.null(response)) {
    models <- list()
    models[[1]] <- model
    for(i in seq.int(2, length(response))) {
      models[[i]] <- update(model, as.formula(paste(deparse(response[[i]]), "~ .")), fit = FALSE, ...)
    }
    model <- new("MultivariatePointProcess")
    setModels(model) <- models
  }
  
  if(is.list(formula)) {
    models <- list()
    for(i in seq_along(formula)) {
      models[[i]] <- update(model, formula[[i]], fit = FALSE, ...)
    }
    model <- new("MultivariatePointProcess")
    setModels(model) <- models 
  }

  if(fit) {
    model <- ppmFit(model, ...)
  } else {
    model <- computeVar(model, method = "none")
##  TODO: correct to work with multivariate models
##    model@optimResult <- list(value = computeMinusLogLikelihood(model),
##                              counts = c(0, 0),
##                              convergence = NA)
  }
  
  return(model)
}

## TODO: Implement the use of interaction terms.

setMethod("computeBasis", "PointProcessModel",
          function(model, ...) {
       
            ## Basis evaluations are computed if not
            ## already computed, locked and available in
            ## 'model@basisEnv$basis'.
            
            if(!bindingIsLocked("basis", model@basisEnv)) {
              processData <- processData(model)
              DcontinuousVar <- paste(colnames(getValue(processData)), ".d", sep="")
              DpositionVar <- paste(processData@positionVar, ".d", sep="")
              markLevels <- levels(getMarkType(processData, drop = FALSE))
              
              ## The formula object is extracted and decomposed into terms.
              ## Each term label (in 'termLabels') is processed below,
              ## and the corresponding basis functions are computed.
              
              mt <- delete.response(terms(formula(model)))
              termLabels <- attr(mt, "term.labels")
              filterTerms <- numeric()
              for(i in seq_along(termLabels)) {
                term <- termLabels[i]
                variable <- all.vars(mt[i])
                
                if(all(variable %in% c(markLevels, DcontinuousVar, DpositionVar))) {                  
                  if(length(variable) > 1){
                    stop(paste("Basis computations with two or more variables in '", term,
                               "' is currently not supported.", sep=""))
                  } else {
                    filterTerms <- c(filterTerms, i)
                    form <- update(mt[i], ~ . -1)
                    x <- as.data.frame(model@basisPoints)
                    colnames(x) <- variable
                    model@basisEnv$basis[[term]] <- model.matrix(form, x)
                  }
                }
              }
              setFilterTerms(model) <- filterTerms
              lockBinding("basis", model@basisEnv)
            }
            return(model)
          }
          )

setMethod("coefficients", "PointProcessModel",
          function(object, ...){
            return(object@coefficients)
          }
          )

setReplaceMethod("coefficients", c(model = "PointProcessModel", value = "numeric"),
                 function(model, value){
                   nc <- length(model@coefficients)
                   d <- dim(getModelMatrix(model))[2]
                   if (d == nc) {
                     model@coefficients[] <- value
                   } else {
                     nv <- length(value)
                     if (nv == 1 && d > 1) {
                       value <- rep(value, d)
                     } else if(nv != d) {
                       value <- rep(.Machine$double.eps, d)
                       warning(paste("Incorrect length of parameter vector. Parameters all set to", .Machine$double.eps))
                     }
                     if(d > 0) 
                       names(value) <- dimnames(getModelMatrix(model))[[2]]
                     model@coefficients <- value
                   }           
          
                   return(model)
                 }
                 )

setMethod("computeDMinusLogLikelihood", "PointProcessModel",
          function(model, coefficients = NULL, ...){
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            eta <- computeLinearPredictor(model, coefficients, ...)
            
            if(model@family@link == "log") {

              dmll <- as.vector(t(exp(eta)*model@delta) %*% getModelMatrix(model)) -
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
          function(model, coefficients = NULL, method = c("poisson", "empirical"), ...) {
            if(isTRUE(response(model) == ""))
              stop("No response variable specified.")
            
            eta <- computeLinearPredictor(model, coefficients, ...)
          
            if(model@family@link == "log"){

              w <-  exp(eta)*model@delta

             } else if(method[1] == "empirical") {
               if(model@family@link == "identity"){

                 points <- getPointPointer(processData(model), response(model))
                 etaP <- eta[points]
                 w <- rep(0, length(eta))
                 w[points] <- 1/etaP^2
                 
               } else {
                 
                 points <- getPointPointer(processData(model), response(model))
                 etaP <- eta[points]
                 w <- model@family@D2phi(eta)*model@delta
                 w[points] <- w[points] - (model@family@D2phi(etaP)*model@family@phi(etaP) - model@family@Dphi(etaP)^2)/model@family@phi(etaP)^2
                 
               }
                              
             } else if(method[1] == "poisson") {

               w <- (model@family@Dphi(eta)^2*model@delta)/model@family@phi(eta)
               
             } else {
               stop(paste("Method", method[1], "is not a valid weight method."))
             }

            ## Negative and small weights are set equal to 0. 
            
            w[w < .Machine$double.eps] <- 0
            
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
            ## Small weights should not play a role. For numeric
            ## stability, they are put equal to 1 here.
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
            
            return(list(weights = w,
                        workingResponse = wr
                        )
                   )
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

setMethod("computeModelMatrix", "PointProcessModel",
          function(model, evaluationPositions = NULL, ...){

            ## Making sure that the basis evaluations are computed.
            ## The call to computeBasis is invoked for its side effect
            ## of computing the basis evaluations if that is not already
            ## done (which is checked by checking if the 'basis' symbol
            ## is locked in the environment 'basisEnv').

            model <- computeBasis(model)

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
              zero <- which(model@basisPoints == 0) - 1
            } else {
              zero <- 0
            }
              
            ## The observed points ('positions') for the marked point process,
            ## the corresponding 'id' labels and 'marks' are extracted.
            
            processData <- processData(model)
            positions <- getPointPosition(processData)
            DcontinuousVar <- paste(colnames(getValue(processData)), ".d", sep="")
            DpositionVar <- paste(processData@positionVar, "d.", sep = "")
            
            id <- factor(getPointId(processData))
            idLevels <- levels(id)
            
            marks <- getMarkType(processData, drop = FALSE)
            markLevels <- levels(marks)
            
            ## The formula object is extracted and decomposed into terms.
            ## Each term label (in 'termLabels') is processed below,
            ## and the corresponding columns in the model matrix are computed.
            
            mt <- delete.response(terms(formula(model)))
            termLabels <- attr(mt, "term.labels")
            filterTerms <- getFilterTerms(model)
            
            ## The points where the basis functions are evaluated are extracted
            ## and the list of model matrices ('design') is set up, which holds model
            ## matrices for the different terms. 'assign' will be an attribute to  
            ## the model matrix of length equal to the number of columns, and for
            ## each column pointing to the term number. 
                        
            ## Model matrix computations for the terms involving filters:

            design <- lapplyParallel(filterTerms,
                                     function(i, ...) {
                                       term <- termLabels[i]
                                       variable <- all.vars(mt[i])
                                       
                                       ## The model matrix computed separately for terms
                                       ## involving marks and terms being linear filters of
                                       ## continuous processes and is done by a loop over 
                                       ## each value of 'id' whose result is stored in
                                       ## 'designList'.
                                       
                                       ## TODO: Can the C level computation return a sparse matrix
                                       ## directly?
                                       
                                       if(all(variable %in% markLevels))
                                         {
                                           
                                           designList <- list()

                                           ## Central loop over 'idLevels' and computations of
                                           ## the model matrix in the C function
                                           ## 'computeFilterMatrix'. Result is converted to a
                                           ## sparse matrix, bound together in one matrix below
                                           ## and stored in the variable 'localDesign'.
                                           
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
                                           localDesign <- do.call("rBind", designList)
                                           colnames(localDesign) <- colnames(getBasis(model, term))
                                         } else if(all(variable %in% c(DcontinuousVar, DpositionVar))) {

                                           designList <- list()    

                                           ## Central loop over 'idLevels' and computations of
                                           ## the model matrix in the C function
                                           ## 'computeFilterMatrix'. Result is converted to a
                                           ## sparse matrix, bound together in one matrix below
                                           ## and stored in the variable 'localDesign'.
                                           
                                           if(variable == DpositionVar) {
                                             values <- getPosition(processData)
                                           } else {
                                             values <- getNumerics(processData)[ , variable == DcontinuousVar]
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
                                           localDesign <- do.call("rBind", designList)
                                           colnames(localDesign) <- colnames(getBasis(model,term))
                                         } 
                                       localDesign ## The return value
                                     },
                                     ## With multicore backend the
                                     ## jobs are spawned sequentially
                                     ## and not prescheduled. This is
                                     ## expected to be a sensible
                                     ## strategy here because
                                     ## different types of filter
                                     ## terms may involve highly
                                     ## different computation times,
                                     ## and the total number of terms
                                     ## is generally relatively small.
                                     mc.preschedule = FALSE 
                                     ) ## End lapplyParalle


            
            assign <- unlist(lapply(filterTerms,
                                    function(i) {
                                      rep(i, dim(getBasis(model, termLabels[i]))[2])
                                    }
                                    )
                             )
            names(design) <- termLabels[filterTerms]
            
            ## Terms that do not involve filters
            if(length(filterTerms) > 0) {
              notFilterTerms <- seq_along(termLabels)[-filterTerms]
            } else {
              notFilterTerms <- seq_along(termLabels)
            }

            ## Model matrix computations for terms involving 'id',
            ## 'position/time', unit variables and non-filtered
            ## continuous time process components.
              
            if(length(notFilterTerms) > 0){
              form <-  mt[notFilterTerms]
              attr(form, "intercept") <-  attr(mt, "intercept")
              variables <- all.vars(form)
              otherVariables <- c(colNames(processData, type = "numeric"),
                                  colNames(processData, type = "factor"),
                                  colNames(processData, type = "unit"))
              
              if(all(variables %in% c(processData@idVar,
                                      processData@positionVar,
                                      otherVariables))) {
                
                ## TODO: Investigate if this build of the model matrix can exploit
                ## sparse matrices more directly, if it is necessary to create the
                ## values object, or if we could refer to a suitable environment?
                
                values <- vector("list", 2)

                if(processData@idVar %in% variables) {
                  values[[1]] <- getId(processData)
                  names(values)[1] <- processData@idVar
                }
                
                if(processData@positionVar %in% variables) {
                  values[[2]] <- getPosition(processData)
                  names(values)[2] <-  processData@positionVar
                }
               
                otherVariables <-  otherVariables[otherVariables %in% variables]
                if(length(otherVariables) > 0) {
                  values <- c(values, getColumns(processData, otherVariables, drop = FALSE))
                  names(values)[-c(1, 2)] <- otherVariables
                }
                
                values <- values[!sapply(values, is.null)]

                ## This sparse.model.matrix solution avoids the dense
                ## model matrix, but does not produce the assign
                ## attribute ...
                ## modelMatrix0 <- sparse.model.matrix(form, values)

                tmp <- model.matrix(form, values)
                assign <- c(c(0, notFilterTerms)[attr(tmp, "assign") + 1], assign)
                modelMatrix0 <- Matrix(tmp, dimnames = dimnames(tmp), sparse=TRUE)
                assign <- c(c(0, notFilterTerms)[attr(modelMatrix0, "assign") + 1], assign)
              } else {
                stop(paste("Use of non existing variable(s) in:", form))
              }
            } else if (attr(mt, "intercept") == 1) {
              modelMatrix0 <- Matrix(rep(1, dim(processData)[1]), sparse=TRUE)
              colnames(modelMatrix0) <- "(Intercept)"
              assign <- c(0, assign)
            } 
            
            if(exists("modelMatrix0")) {
              modelMatrix <- cBind(modelMatrix0, do.call("cBind", design))
            } else {
              modelMatrix <- do.call("cBind", design)
            }
            form <- formula(model)
            attr(form, "filterTerms") <- which(!(seq(along=termLabels) %in% notFilterTerms))
            model <- updateModelMatrix(model, modelMatrix, assign, form)
            lockEnvironment(model@modelMatrixEnv, binding = TRUE)
            return(model)
          }
          )

setMethod("computeVar", "PointProcessModel",
          function(model, method = attr(vcov(model), "method"), ...){
            switch(method,
                   subset = {
                     varMatrix <- vcov(model)
                     i <- which(rownames(varMatrix) %in% names(coefficients(model)))
                     model@var <- varMatrix[i, i, drop = FALSE]
                   },
                   none = {
                     model@var <- matrix(0, length(coefficients(model)),
                                         length(coefficients(model)))
                   },
                   lsSandwich = {
                     X <- ppstat:::getModelMatrix(model)
                     points <- getPointPointer(processData(model),
                                                      ppstat:::response(model))
                     vcov <- matrix(0, nrow = dim(X)[2], ncol = dim(X)[2])
                     fixed <- model@fixedCoefficients
                     if(length(fixed) == 0) {
                       Iinv <- crossprod(sqrt(model@delta) * X)
                       K <- crossprod(X[points, ])
                     } else {
                       Iinv <- crossprod(sqrt(model@delta) * X[, -fixed$which])
                       K <- crossprod(X[points, -fixed$which])
                     }                       
                     Iinv <- try(solve(Iinv), silent = TRUE)
                     if(class(Iinv) == "try-error") 
                       message("Model matrix not of full column rank:\n", Iinv[1], " Check parameterization.")
                     vcov <- as(Iinv %*% K %*% Iinv, "matrix") ### Sandwich estimator
                     rownames(vcov) <- names(model@coefficients)
                     colnames(vcov) <- names(model@coefficients)
                     model@var <- vcov
                   },
                   Fisher = {
                     vcovInv <- computeDDMinusLogLikelihood(model)
                     if(model@penalization)
                       vcovInv <- vcovInv + 2*model@Omega
                     ## TODO: This requires some more thought ....
                     vcov <- matrix(0, nrow = dim(vcovInv)[1],
                                    ncol=dim(vcovInv)[2])
              
                     if(length(model@fixedCoefficients) == 0) {
                       tmp <- try(solve(vcovInv), silent = TRUE)
                       if(class(tmp) == "try-error") {
                         message("Fisher information singular:\n", tmp[1], " Check convergence status or parameterization.")
                       } else {
                         vcov <- tmp
                       }
                     } else {
                       tmp <- try(solve(vcovInv[-model@fixedCoefficients$which,
                                                -model@fixedCoefficients$which]),
                                  silent=TRUE)
                       if(class(tmp) == "try-error") {
                         message("Fisher information singular:\n", tmp[1], " Check convergence status and parameterization.")
                       } else {
                         vcov[-model@fixedCoefficients$which,
                              -model@fixedCoefficients$which] <- tmp
                       }
                     }
                     
                     rownames(vcov) <- names(model@coefficients)
                     colnames(vcov) <- names(model@coefficients)
                     model@var <- (vcov + t(vcov))/2   ## To ensure symmetry
                   }
                   )
            return(model)
          }
          )

setMethod("getFilterTerms", "PointProcessModel",
          function(model, ...) {
            return(model@filterTerms)
          }
          )

setMethod("getLinearFilter", "PointProcessModel",
          function(model, se = FALSE, nr, ...){
            mt <- delete.response(terms(formula(model)))
            model <- computeBasis(model)
            filterTerms <- getFilterTerms(model)
            
            linearFilter <- list()
            design <- list()
            if(isTRUE(se))
              linearFilterSE <- list()
            
            if(missing(nr)) {
              nr <- length(model@basisPoints)
            } else {
              nr <- min(nr, length(model@basisPoints))
            }

            for(j in filterTerms){
              term <- attr(mt[j], "term.labels")
              NR <- dim(getBasis(model, term))[1]
              i <- round(seq_len(min(nr, NR))*NR/nr)
              varName <- paste(all.vars(parse(text = term)), collapse = ".")
              design[[varName]] <- cbind(design[[varName]], getBasis(model, term)[i, ,drop=FALSE])}

            for(j in seq_along(design)){
              linearFilter[[j]] <- design[[j]] %*% coefficients(model)[dimnames(design[[j]])[[2]]]
              if(isTRUE(se))
                linearFilterSE[[j]] <- sqrt(rowSums(design[[j]] %*% vcov(model)[dimnames(design[[j]])[[2]], dimnames(design[[j]])[[2]]] * design[[j]]))
            }

            names(linearFilter) <- names(design)
          
            if(isTRUE(se)) {
              return(list(linearFilter = cbind(data.frame(x=model@basisPoints[i]), as.data.frame(linearFilter)), se = linearFilterSE))
            } else {
              return(cbind(data.frame(x = model@basisPoints[i]), as.data.frame(linearFilter)))
            }
              
          }
          )

setMethod("getAssign", "PointProcessModel",
          function(model, col, ...){
            if(missing(col))
              col <- model@modelMatrixCol
            if(length(col) == 0) {
              assign <- model@modelMatrixEnv$assign
            } else {
              assign <- model@modelMatrixEnv$assign[col]
            }
            return(assign)
          }
          )

setMethod("getBasis", c(model = "PointProcessModel", term = "ANY"),
          function(model, term,...){
            if(missing(term))
              return(model@basisEnv$basis)
            return(model@basisEnv$basis[[term]])
          }
          )

setMethod("getModelMatrix", c(model = "PointProcessModel", col = "ANY"),
          function(model, col,...){
            if(missing(col))
              col <- model@modelMatrixCol
            if(length(col) == 0) {
              modelMatrix <- model@modelMatrixEnv$modelMatrix
            } else {
              modelMatrix <- model@modelMatrixEnv$modelMatrix[, col, drop = FALSE]
            }
            return(modelMatrix)
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
                     model@basisEnv <- new.env(parent = emptyenv())
                   
                   model@basisEnv$basis[[term]] <- value
                   lockEnvironment(model@basisEnv)
                   
                   return(model)
                 }
                 )

setReplaceMethod("setBasis", c(model = "PointProcessModel", term = "ANY", value = "list"),
                 function(model, term, value){
                   if(environmentIsLocked(model@basisEnv))
                     model@basisEnv <- new.env(parent = emptyenv())
                   
                   model@basisEnv$basis <- value
                   lockEnvironment(model@basisEnv)

                   return(model)
                 }
                 )

setReplaceMethod("setModelMatrixEnv", c(model = "PointProcessModel", value = "list"),
                 function(model, value){
                   if(all(c("modelMatrixEnv", "modelMatrixCol") %in% names(value))){
                     model@modelMatrixEnv <- value$modelMatrixEnv
                     model@modelMatrixCol <- value$modelMatrixCol
                   } else {
                     stop("Right hand side of the assignment needs to be a list with two entries named 'modelMatrixEnv' and 'modelMatrixCol'.")
                   }
                   return(model)
                 }
                 )

setMethod("predict", "PointProcessModel",
          function(object, ...) {
            eta <- computeLinearPredictor(object, ...)
            return(object@family@phi(eta))
          }
          )

setMethod("getTermPlotData", "PointProcessModel",
          function(model, alpha = 0.05, trans = NULL, ...) {
            if(alpha <= 0 || alpha > 1)
              stop("The 'alpha' level must be in (0,1]")
             
            if(isTRUE(all.equal(alpha, 1))) {
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

            return(plotData)
          }
          )
          

setMethod("termPlot", "PointProcessModel",
          function(model, alpha = 0.05, layer = geom_line(), trans = NULL, ...) {
            if(length(getFilterTerms(model)) == 0){
              print("No filter function terms to plot")
              return(invisible())
            }

            plotData <- getTermPlotData(model = model, alpha = alpha, trans = trans, ...)
                      
            linearFilterPlot <- ggplot(data = plotData, aes(x = x, y = value)) +
              facet_grid(variable ~ ., scales = "free_y") +
                scale_x_continuous("position") +
                  scale_y_continuous("") + layer
            
            if(!isTRUE(all.equal(alpha, 1)))
              linearFilterPlot <- linearFilterPlot +
                geom_ribbon(aes(min = cf.lower, max = cf.upper),
                            fill = "blue", alpha = 0.2)

            return(linearFilterPlot)
          }
          )

setMethod("ppmFit", "PointProcessModel",
          function(model, control = list(), optim = 'optim', ...) {
            model <- switch(optim,
                            optim = optimFit(model = model, control = control, ...),
                            iwls = iwlsFit(model = model, control = control, ...),
                            ls = lsFit(model = model, control = control, ...),
                            poisson = glmFit(model = model, control = control, ...),
                            glmnet = glmnetFit(model = model, control = control, ...),
                            stop(paste("No optimization method '", optim,
                                       "' available.", sep = ""), call. = FALSE)
                            )
            ## Computation of the estimated covariance matrix
            computeVar(model)
          }
          )
              
setMethod("optimFit", "PointProcessModel",
          function(model, control = list() , ...) {
            nrPar <- parDim <- dim(getModelMatrix(model))[2]
            fixedPar <- model@fixedCoefficients
            Omega <- model@Omega

            ## If the model is a submodel we temporarily change the
            ## full model matrix to the submodel matrix in this function
            ## to avoid repeated subsetting of the full matrix. 
            
            if(length(model@modelMatrixCol) > 0) {
              modelMatrixEnv <- getModelMatrixEnv(model)
              model <- updateModelMatrix(model)
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
              initPar <- rep(.Machine$double.eps, nrPar)             
            warning("Length of initial parameter vector wrong. Initial parameters all set to 0.")
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
              if(attr(terms(formula(model)), "intercept")==1) {
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
            
            return(model)
          }
          )            

setMethod("grpLassoFit", "PointProcessModel",
          function(model, control, ...) {
            
            w <- computeWeights(model)
            z <- computeWorkingResponse(model)
            
          }
          )

setMethod("iwlsFit", "PointProcessModel",
          function(model, control = list(), ...) {

            value <- computeMinusLogLikelihood(model)
            reltol <- sqrt(.Machine$double.eps)
            maxit <- 1000  ## Currently hardcoded!
            trace <- 1
            
            i <- 1
            while(i < maxit) {
              ## Computation of weights and working response returned in a list
              ## with entries 'weights' and 'workingReponse'.
              wr <- computeWorkingResponse(model)
              ## TODO: check if lm.fit.sparse is exported. The
              ## following computation relies on an algorithm in the
              ## MatrixModels package still under development and
              ## currently not exported.
              coefficients(model) <- MatrixModels:::lm.fit.sparse(getModelMatrix(model),
                                                                  wr$workingResponse,
                                                                  wr$weights)
              val <- computeMinusLogLikelihood(model)
              i <- i + 1
              if(trace > 0)
                cat("Value: ", val, "\n")
              if(val < value && value < reltol * (abs(value) + reltol) + val)
                break
              value <- val
            }

            if(i < maxit) {
              convergence <- 0
            } else {
              convergence <- 1
            }

            model@optimResult <- list(value = value,
                                      counts = c(i, 0),
                                      convergence = convergence
                                      )
            return(model)
          }
          )

setMethod("glmnetFit", "PointProcessModel",
          function(model, control = list(), refit = FALSE, ...) {
            hasGlmnet <- require("glmnet")
            if(!hasGlmnet) {
              stop("Package 'glmnet' is not installed.")
            } else {
              offset <- log(model@delta)
              weights <- rep(1, length(offset))
              weights[offset == -Inf] <- 0
              offset[offset == -Inf] <- 0
              y <- rep(0, length(offset))
              y[getPointPointer(processData(model), ppstat:::response(model))] <- 1
              ## TODO: How do we automatically compute a lambda sequence
              ##              lambda <- exp(seq(-4, -15, -0.1))

              glmnetFit <- glmnet(x = ppstat:::getModelMatrix(model),
                                  y = y,
                                  family = "poisson",
                                  offset = offset,
                                  weights = weights,
                                  standardize = FALSE,
                                  alpha = 1)
              browser()
              coef <- coefficients(glmnetFit)
              mll <- sapply(seq(1, dim(coef)[2]),
                            function(i) {
                              computeMinusLogLikelihood(testModel, coef[ ,i])
                            }
                            )
              
              coefficients(model) <- coefficients(glmnetFit, "lambda.min")[, 1]
              nzCoef <- tapply(coefficients(model),
                               ppstat:::getAssign(model),
                               function(x) any(x != 0)
                               )
              
              nzCoef <- which(nzCoef[names(nzCoef) != "0"])
              
              zc <- which(coefficients(model) == 0)
              model@fixedCoefficients <- list(which = zc,
                                              values = rep(0, length(zc)))
              model <- update(model, nzCoef, fit = refit)
              model <- ppstat:::computeVar(model)
              
              optimResult <- list(value = computeMinusLogLikelihood(model),
                                  counts = c(glmnetFit$npasses, 0),
                                  convergence = 0
                                  )
              model@optimResult <- optimResult
            }
            return(model)
          }
          )

setMethod("lsFit", "PointProcessModel",
          function(model, control = list(), ...) {

            if(model@family@link == "identity") {
              if(model@varMethod != "none")
                model@varMethod <- "lsSandwich"
              fixedPar <- model@fixedCoefficients
              w <- sqrt(model@delta)
              X <- w * ppstat:::getModelMatrix(model)
              y <- rep(0, dim(X)[1])
              if (model@penalization) {
                if (length(fixedPar) == 0) {
                  Omega <- model@Omega
                }
                else {
                  Omega <- model@Omega[-fixedPar$which, -fixedPar$which]
                }
                OmegaSVD <- svd(Omega, 0)
                L <- Matrix(sqrt(OmegaSVD$d) * OmegaSVD$v)
                y <- c(y, rep(0, dim(Omega)[1]))
                X <- rBind(X, L)
              }
              points <- getPointPointer(processData(model), ppstat:::response(model))
              y[points] <- 1/w[points]
              if (length(fixedPar) == 0) {
                model@coefficients <- MatrixModels:::lm.fit.sparse(X, 
                                                                   y)
              }
              else {
                model@coefficients[-fixedPar$which] <- MatrixModels:::lm.fit.sparse(X[, 
                                                                                      -fixedPar$which], y)
                model@coefficients[fixedPar$which] <- fixedPar$value
              }
              names(model@coefficients) <- dimnames(X)[[2]]
              model@optimResult <- list(value = computeMinusLogLikelihood(model),
                                        counts = c(1, 0),
                                        convergence = 0
                                        )
            } else {
              stop("Use of 'lsFit' only supported with the identity link function.")
            }
            return(model)
          }
          )

setMethod("glmFit", "PointProcessModel",
          function(model, control = list(), ...) {
            offset <- log(model@delta)
            weights <- rep(1, length(offset))
            weights[offset == -Inf] <- 0
            offset[offset == -Inf] <- 0
            y <- rep(0, length(offset))
            y[getPointPointer(processData(model), ppstat:::response(model))] <- 1
            glmFit <- glm.fit(x = as(ppstat:::getModelMatrix(model), "matrix"),
                              y = y,
                              family = poisson(),
                              offset = offset,
                              weights = weights,
                              control = control)
            
            coefficients(model) <- glmFit$coefficients
            
            if(isTRUE(glmFit$converged)) {
                convergence <- 0
              } else {
                convergence <- glmFit$converged
              }
            
            optimResult <- list(value = computeMinusLogLikelihood(model),
                                counts = c(glmFit$iter, 0),
                                convergence = convergence
                                )
            model@optimResult <- optimResult
            
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
            z[se == 0] <- NA
            q[se == 0] <- NA
            se[se == 0] <- NA
            
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
            attr(object@var, "method") <- object@varMethod
            return(object@var)
          }
          )

setMethod("update", "PointProcessModel",
          function(object, formula = .~., warmStart = TRUE, fixedCoefficients = list(), fit = TRUE, ...){
            
            modelFormula <- formula(object)
            if(class(formula) != "formula") {
              selectTerms <- as.numeric(formula)
              selectTerms <- selectTerms[selectTerms > 0]
              formula <- formula(terms(modelFormula)[selectTerms])
            }
            superFormula <- terms(object@modelMatrixEnv$formula)
            updatedFormula <- update(modelFormula, formula)
            updatedTermLabels <- attr(terms(updatedFormula), "term.labels")
            superTermLabels <- attr(superFormula, "term.labels")
            superFilterTerms <- superTermLabels[attr(object@modelMatrixEnv$formula,
                                                     "filterTerms")]
            setFilterTerms(object) <- which(updatedTermLabels %in%
                                            superFilterTerms)
            
            if(attr(terms(updatedFormula), "intercept") == 1)
              updatedTermLabels <- c("Intercept", updatedTermLabels)
            if(attr(superFormula, "intercept") == 1)
              superTermLabels <- c("Intercept", superTermLabels)
           
            formula(object) <- updatedFormula
            
            if(length(object@modelMatrixCol) == 0){
              tmpCoef <- coefficients(object)
            } else {
              tmpCoef <- rep(.Machine$double.eps, dim(getModelMatrix(object, numeric()))[2])
              tmpCoef[object@modelMatrixCol] <- coefficients(object)
            }

            if(all(updatedTermLabels %in% superTermLabels)) {
              if("Intercept" %in% superTermLabels) {
                col <- which(superTermLabels[getAssign(object, numeric()) + 1] %in% updatedTermLabels)
              } else {
                col <- which(superTermLabels[getAssign(object, numeric())] %in% updatedTermLabels)
              }
              ## TODO: The following checks whether there are
              ## interactions in the model. There is a general problem
              ## with choices of contrasts in submodels that should be
              ## dealt with. Removing e.g. the intercept gives a
              ## peculiar model it seems ...
              if(any(attr(terms(modelFormula), "order") >= 2)) {
                warning(paste(c(object@call, "Original model formula includes interaction terms. Check that updated model is as expected."),
                              collapse="\n"), call. = FALSE)
              }
              if(length(col) == length(tmpCoef)) {
                object@modelMatrixCol <- numeric()
              } else {
                object@modelMatrixCol <- col
              }
              if(warmStart) {
                coefficients(object) <- tmpCoef[col]
              } else {
                coefficients(object) <- rep(0, length(col))
              }
            } else {
              unlockBinding("basis", object@basisEnv)
              object <- computeModelMatrix(object)
              coefficients(object) <- rep(.Machine$double.eps,
                                          dim(getModelMatrix(object))[2])
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

            if(fit) {
              object <- ppmFit(object, ...)
            } else {
              object <- computeVar(object, method = "subset")
            }
            
            return(object)
          }
          )

setMethod("getInformation", "PointProcessModel",
          function(model, k = 2, ...) {
            df <- length(model@coefficients) - length(model@fixedCoefficients$which)
            mll <- computeMinusLogLikelihood(model)
            return(2*mll + k*df)
          }
          )

setReplaceMethod("setFilterTerms", "PointProcessModel",
                 function(model, value) {
                   model@filterTerms <- value
                   return(model)
                 }
                 )

## TODO: Is it possible to implement the following using parallelization?

setMethod("stepInformation", "PointProcessModel",
          function(model, scope = ~ 0, direction = "both", trace = 1,
                   steps = 1000, warmStart = TRUE, k = 2, ...) {
            if(trace == 1) 
              cat(" Step \t  AIC \t\t Direction\n")
            
            AIC <- list()
            models <- list()
            models[["current"]] <- model
            models$current@varMethod <- 'none'
            step <- 1
            termsToAdd <- attr(terms(scope), "term.labels")
            initialForm <- terms(formula(model))
            initialTerms <- attr(initialForm, "term.labels")
            currentTerms <- initialTerms
            termsToRemove <- currentTerms
            minDirection <- "start"
            while(step < steps) {
              AIC$current <- getInformation(models$current, k, ...)

              if(trace == 1)
                cat(" ", step, "\t", AIC$current, "\t ", minDirection, "\n")

              if(trace == 2) {
                cat("Step:", step, "\tAIC:", AIC$current,"\n")
                cat("Model:", deparse(formula(models$current)),"\n")
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


setMethod("updateModelMatrix", "PointProcessModel",
          function(model, modelMatrix = getModelMatrix(model), assign = getAssign(model), form){
            force(modelMatrix)
            model@modelMatrixEnv <- new.env(parent = emptyenv())
            model@modelMatrixEnv$modelMatrix <- modelMatrix
            model@modelMatrixEnv$assign <- assign
            if(!missing(form))
              model@modelMatrixEnv$formula <- form
            model@modelMatrixCol <- numeric()
            return(model)
          }
          )
