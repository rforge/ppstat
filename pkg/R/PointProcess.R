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
             if(isTRUE(response(model) == ""))
               stop("No response variable specified.")
            
            if(model@family@link == "log"){
              mll <-  sum(exp(eta)*model@delta) -
                sum(eta[getPointPointer(processData(model), response(model))]) 
            } else {
              mll <-  sum(model@family@phi(eta)*model@delta) -
                sum(log(model@family@phi(eta[getPointPointer(processData(model), response(model))]))) 
            }
            
            return(mll)
          }
          )

setMethod("family", "PointProcess",
          function(object,...) {
            return(object@family)
          }
          )

setMethod("formula", "PointProcess",
          function(x, ...){
            return(x@formula)
          }
          )


setMethod("response", "PointProcess",
          function(model, ...){
            return(model@response)
          }
          )

setReplaceMethod("formula", c(model = "PointProcess", value = "formula"),
                 function(model, value){
                   if(attr(terms(value), "response") != 0) {
                     ## TODO: Is the response always in position 2 in this list/call?
                     response <- attr(terms(value), "variables")[[2]]
                     model@response <- all.vars(response)                  
                   } else {
                     model@response <- ""
                   }
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




