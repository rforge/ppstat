setMethod("anticipating", "PointProcess",
          function(model, ...) {
            antipmodels = c("Gibbs")
            if(family(model@family) %in% antipmodels)
              return(TRUE)

            return(FALSE)
          }
          )

setMethod("computeMinusLogLikelihood", "PointProcess",
          function(model, coefficients = NULL, ...) {
            eta <- computeLinearPredictor(model, coefficients, ...)
             if(attr(terms(model@formula), "response") != 0) {
              response <- all.vars(model@formula, unique = FALSE)[attr(terms(model@formula), "response")]
            } else {
              stop("No response variable specified.")
            }
            
            if(model@family@link == "log"){
              mll <-  sum(exp(eta)*model@delta) -
                sum(eta[getPointPointer(processData(model), response)]) 
            } else {
              mll <-  sum(model@family@phi(eta)*model@delta) -
                sum(log(model@family@phi(eta[getPointPointer(processData(model), response)]))) 
            }
            
            return(mll)
          }
          )

setMethod("formula", "PointProcess",
          function(x,...){
            return(x@formula)
          }
          )

setReplaceMethod("formula", c(model = "PointProcess", value = "formula"),
                 function(model, value){
                   model@formula <- value
                   return(model)
                 }
                 )

setMethod("processData", "PointProcess",
          function(model, ...){
            return(model@processData)
          }
          )

setReplaceMethod("processData", c(model = "PointProcess", value = "MarkedPointProcess"),
          function(model, value){
            model@processData <- value
            return(model)
          }
          )




