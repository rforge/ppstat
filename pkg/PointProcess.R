setClass("PointProcess",
         representation(
                        processData = "ProcessData",
                        delta = "numeric",
                        formula = "formula",
                        family = "Family",
                        call = "call",
                        support = "numeric",     ### A vector c(a,b) with the support, [a,b], of the g-functions
                        Delta = "numeric",       ### The equidistant spacing between the g-evaluations
                        "VIRTUAL")
         )

setMethod("formula","PointProcess",
          function(x,...){
            return(x@formula)
          }
          )

setReplaceMethod("formula","PointProcess",
          function(.Object,value){
                        
            if(is.formula(value)) {
              .Object@formula <- value
            } else stop("formula needs to be a 'formula'")

            return(.Object)
          }
          )

setMethod("computeMinusLogLikelihood","PointProcess",
          function(object,coefficients=NULL,...){
            eta <- computeLinearPredictor(object,coefficients,...)
             if(attr(terms(object@formula),"response") != 0) {
              response <- all.vars(object@formula,unique=FALSE)[attr(terms(object@formula),"response")]
            } else stop("no response variable specified")
            

            if(object@family@link == "log"){
              mll <-  sum(exp(eta)*object@delta) -
                sum(eta[getMarkTypePosition(object@processData,response)]) 
            } else {
              mll <-  sum(object@family@phi(eta)*object@delta) -
                sum(log(object@family@phi(eta[getMarkTypePosition(object@processData,response)]))) 
            }
            
            return(mll)
            
          }
          )


