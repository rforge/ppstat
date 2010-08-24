setClass("PointProcess",
         representation(
                        call = "call",

                        ## The equidistant spacing between the g-evaluations.
                        Delta = "numeric",      

                        delta = "numeric",
                        family = "Family",
                        formula = "formula",

                        ## Results from call to 'optim' goes here.
                        optimResult = "list",
                        
                        ## The processData as ProcessData is in this environment. Locked after initialization.
                        processDataEnv = "environment",
                        
                        ## A vector c(a,b) with the support, [a,b], of the g-functions.
                        support = "numeric",     

                        "VIRTUAL")
         )

setMethod("computeMinusLogLikelihood","PointProcess",
          function(model,coefficients=NULL,...){
            eta <- computeLinearPredictor(model,coefficients,...)
             if(attr(terms(model@formula),"response") != 0) {
              response <- all.vars(model@formula,unique=FALSE)[attr(terms(model@formula),"response")]
            } else stop("no response variable specified")
            

            if(model@family@link == "log"){
              mll <-  sum(exp(eta)*model@delta) -
                sum(eta[getMarkTypePosition(getProcessData(model),response)]) 
            } else {
              mll <-  sum(model@family@phi(eta)*model@delta) -
                sum(log(model@family@phi(eta[getMarkTypePosition(getProcessData(model),response)]))) 
            }
            
            return(mll)
            
          }
          )

setMethod("formula","PointProcess",
          function(x,...){
            return(x@formula)
          }
          )

setReplaceMethod("formula",c(model="PointProcess",value="formula"),
                 function(model,value){
                   model@formula <- value
                   model@basisEnv <- new.env(parent=.GlobalEnv)
                   model@basisEnv$basis <- list()
                   return(model)
                 }
                 )

setMethod("getProcessData","PointProcess",
          function(model,...){
            return(model@processDataEnv$processData)
          }
          )

setReplaceMethod("setProcessData",c(model="PointProcess",value="ProcessData"),
          function(model,value){
            processDataEnv <- new.env(parent=.GlobalEnv)
            processDataEnv$processData <- value
            lockEnvironment(processDataEnv,bindings=TRUE)
            model@processDataEnv <- processDataEnv
            return(model)
          }
          )




