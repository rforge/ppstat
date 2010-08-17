### TODO: The current implementation below is a preliminary attempt to make the
### general gradient algorithm work. This should change into a class
### relying on reproducing kernels.

setClass("PointProcessKernel",
         representation(
                        modelMatrix = "Matrix",
                        recurrenceMatrix = "list", ### This is a list of sparse matrices of class "Matrix"
###                        integrationMatrix = "list"
                        labels = "list",          ### a list with two entries holding the term labels and idLevels
                                                  ### used in the namings of the entries in the recurrenceMatrix
                        g = "matrix",            ### A matrix of g-evalutions at 0,Delta,2Delta,...
                        coefficients = "numeric",
                        Omega = "matrix",
                        penalization = "logical",
                        var = "matrix",
                        optimResult="list"
                        ),
         contains="PointProcess")

setMethod("initialize","PointProcessKernel",
          function(.Object,
                   processData,
                   formula,
                   family,
                   support,
                   Delta,
                   fit=TRUE,
                   recurrenceMatrix = TRUE,
                   g = NULL,
                   Omega=NULL,
                   call=NULL,...){
            .Object@processDataEnv <- new.env(.GlobalEnv)
            .Object@processDataEnv$processData <- processData
            lockEnvironment(.Object@processDataEnv,bindings=TRUE)
            .Object@formula <- formula
            .Object@family <- family
            .Object@support <- support
            .Object@Delta <- Delta
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
            
            .Object@delta <- as.numeric(unlist(tapply(getPosition(getContinuousProcess(processData)),
                                                      getId(getContinuousProcess(processData)),
                                                      function(x) c(diff(x),0)),use.names=FALSE))
            
            if(recurrenceMatrix) {
              .Object <- computeRecurrenceMatrix(.Object)
            }

            if(is.null(g)) {
              if(recurrenceMatrix) {
                ##        dimensions <- lapply(.Object@integrationMatrix,function(x) dim(x)[2])
                ##                 nrTerms <- max(dimensions)
                g <- matrix(0,nrow=1+supportBound/Delta,ncol=length(.Object@labels$termLabels))
                colnames(g) <- .Object@labels$termLabels 
              } else {
                g <- 0
              }
            } else {
              g <- g
            }
              
            if(fit){
              .Object <- glppsFit(.Object,...)
            }
          
          
          if(!is.null(call)) .Object@call <- call
          return(.Object)
          }
          )

setMethod("coefficients","PointProcessKernel",
          function(object,...){
            return(list(object@g,object@coefficients))
          }
          )

setMethod("computeRecurrenceMatrix","PointProcessKernel",
          function(object,evaluationPositions=NULL,...){

            fR <- function(s) s*as.numeric(s >= object@support[1] & s <= object@support[2])

            design <- list()
            designList <- list()
            
            if(is.null(evaluationPositions)) {
              evalPositions <- tapply(getPosition(getContinuousProcess(getProcessData(object))),
                                      getId(getContinuousProcess(getProcessData(object))),list)
            } else {
              evalPositions <- evaluationPositions
            }
            
            markedPointProcess <- getMarkedPointProcess(getProcessData(object))
            positions <- getPosition(markedPointProcess)

            id <- factor(getId(markedPointProcess))
            idLevels <- levels(id)

            marks <- getMarkType(markedPointProcess)
            markLevels <- levels(marks)

            mt <- terms(object@formula) 
            termLabels <- attr(mt,"term.labels")
            notMarkTerms <- character()
            
            for(term in termLabels){

              form <- as.formula(paste("~",term,"-1"))
              variables <- all.vars(form)
              mark <- markLevels[markLevels %in% variables]
              
              if(length(mark) > 1) {
                stop(paste("Interaction of two mark types in '",
                           term,"' is currently not implemented.",sep=""))
              } else if(length(mark) == 1) {
                if(mark != term) {
                  stop(paste("Use of '",term,"' is not implemented",sep=""))
                } else {

                  form <- as.formula(paste("~fR(",mark,")-1"))
                  posi <- positions[marks==mark]
                  idi <- id[marks==mark]
                  idList <- factor()
                  last <- as.list(rep(0,length(idLevels)))
                  names(last) <- idLevels
                  for(j in seq(along=posi)){
                    afterPosi <- evalPositions[[idi[j]]] > posi[j]
                    val <- list(evalPositions[[idi[j]]][afterPosi] - posi[j])
                    names(val) <- mark
                    tmp <- model.matrix(form,val)
                    if(!(idi[j] %in% idLevels[idList])) {
                      designList[[paste(term,idi[j],sep="")]] <-  Matrix(0,nrow=length(evalPositions[[idi[j]]]),
                                                                         ncol=length(which(idi == idi[j]))*dim(tmp)[2],
                                                                         sparse=TRUE)
                      idList <- c(idList,idi[j])
                    }
                    designList[[paste(term,idi[j],sep="")]][afterPosi,
                                                            seq(last[[idi[j]]]+1,
                                                                last[[idi[j]]]+dim(tmp)[2])] <- tmp
                    last[[idi[j]]] <- last[[idi[j]]]+dim(tmp)[2]                   
                  }
                } 
              } else if(length(mark) == 0) notMarkTerms <- c(notMarkTerms,term)
            }

##           for(k in idLevels) {
##               labels <- paste(termLabels,k,sep="")
##               labels <- labels[labels %in% names(designList)]
##               dimensions <- sapply(designList[[labels]],function(x) dim(x)[2])
##               integration[[k]] <- bdiag(dimensions,function(d) rep(1,d))
##               colnames(integration[[k]]) <- termLabels[labels %in% names(designList)]
##               design[[k]] <- do.call(cBind, designList[[labels]])
##               }

            object@recurrenceMatrix <- designList
            object@labels <- list(termLabels = termLabels,idLevels=idLevels)
###            object@integrationMatrix <- integration

### Model matrix computations for the id and continuous time process
### components
            
            if(attr(mt,"intercept") == 1){
              notMarkTerms <- paste(c(notMarkTerms,"1"),collapse="+")
            } else if(length(notMarkTerms) >0) {
              notMarkTerms <-  paste(paste(notMarkTerms,collapse="+"),"-1")
            }

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
                object@modelMatrix <- Matrix(tmp,dimnames=dimnames(tmp))
              } else {
                stop(paste("Use of non existing variable(s) in:", form))
              }
            }

            object@coefficients <- rep(0,dim(object@modelMatrix)[2])
            
            return(object)
          }
          )
            
setMethod("computeLinearPredictor","PointProcessKernel",
          function(object,coefficients=NULL,...){
            if(is.null(coefficients)) {
              g <- coefficients(object)[[1]]
              coefficients <- coefficients(object)[[2]]
            } else {
              g <- coefficients[[1]]
              coefficients <- coefficients(object)[[2]]              
            }

            predictors <- list()
            for(id in object@labels$idLevels) {
              tmp <- list()
              for(term in object@labels$termLabels) {
                label <- paste(term,id,sep="")
                if(any(label == names(object@recurrenceMatrix))){
                  tmp[[term]] <- object@recurrenceMatrix[[label]]
                  lookup <- ceiling(tmp[[term]]@x/object@Delta)
                  tmp[[term]]@x <- g[lookup,term]
                  tmp[[term]] <- rowSums(tmp[[term]])    ### Is there a performance gain by returning a sparse result?
                                                         ### A sparse result creates problems with computations of row sums below?!
                }
              }
              if(length(tmp)==0) stop("No terms for id ",id," in computation of linear predictor.")
              predictors[[id]] <- do.call(cbind,tmp) 
              predictors[[id]] <- predictors[[id]] %*% rep(1,dim(predictors[[id]])[2])
            }

            eta =  do.call(c,predictors) + as.numeric(object@modelMatrix %*% coefficients)

            return(eta)
          }
          )


setMethod("computeDMinusLogLikelihood","PointProcessKernel",
          function(object,coefficients=NULL,...){
            eta <- computeLinearPredictor(object,coefficients=NULL,...)
            if(attr(terms(object@formula),"response") != 0) {
              response <- all.vars(object@formula,unique=FALSE)[attr(terms(object@formula),"response")]
            } else stop("no response variable specified")
            
            seqDelta <- object@Delta*seq(ceiling(object@support[1]/object@Delta),ceiling(object@support[2]/object@Delta))
            dmll <- NULL ### This is bad, should be fixed!
            
            if(object@family@link == "log") {
              mP <- getMarkTypePosition(getProcessData(object),response)
              weights <- exp(eta)*object@delta
              ones <-  rep(1,length(mP))
              
              for(r in seq(along=seqDelta)){
                intKernels <- list()
                for(id in object@labels$idLevels) {
                  tmp <- list()
                  for(term in object@labels$termLabels) {
                    label <- paste(term,id,sep="")
                    if(any(label == names(object@recurrenceMatrix))){
                      tmp[[term]] <- object@recurrenceMatrix[[label]]
                      tmp[[term]]@x <- R(tmp[[term]]@x,seqDelta[r])
                      tmp[[term]] <- rowSums(tmp[[term]])   
                    }
                  }
                  if(length(tmp)==0) stop("No terms for id ",id," in computation of linear predictor.")
                  intKernels[[id]] <- do.call(cbind,tmp) 
                }
                intKernels <- do.call(rbind,intKernels)
                if(is.null(dmll)) dmll <- matrix(0,nrow=length(seqDelta),ncol=dim(intKernels)[2]) ### Stupid if-sentence
                dmll[r,] <-  weights %*% intKernels - ones %*% intKernels[mP,]
            }

            } else {
              
              etaP <- eta[getMarkTypePosition(getProcessData(object),response)]
              mmP <- object@modelMatrix[getMarkTypePosition(getProcessData(object),response),]

              dmll <-  colSums((object@family@Dphi(eta)*object@delta)*object@modelMatrix) -
                colSums(object@family@Dphi(etaP)/object@family@phi(etaP)*mmP)
              
            }
            
            return(dmll)
            
          }
          )
