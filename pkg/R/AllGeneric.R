### Generics for ppstat

setGeneric("computeLinearPredictor",function(object,...) standardGeneric("computeLinearPredictor"))
setGeneric("computeMinusLogLikelihood",function(object,coefficients=NULL,...) standardGeneric("computeMinusLogLikelihood"))
setGeneric("computeDMinusLogLikelihood",function(object,coefficients=NULL,...) standardGeneric("computeDMinusLogLikelihood"))
setGeneric("computeDDMinusLogLikelihood",function(object,coefficients=NULL,...) standardGeneric("computeDDMinusLogLikelihood"))
setGeneric("computeModelMatrix",function(object,evaluationPositions=NULL,...) standardGeneric("computeModelMatrix"))
setGeneric("computeRecurrenceMatrix",function(object,evaluationPositions=NULL,A=1000,...) standardGeneric("computeRecurrenceMatrix"))
setGeneric("getModelMatrix",function(object,...) standardGeneric("getModelMatrix"))
setGeneric("getLinearFilter",function(object,se,nr,...) standardGeneric("getLinearFilter"))
setGeneric("getProcessData",function(object,...) standardGeneric("getProcessData"))
setGeneric("coefficients<-",function(.Object,value) standardGeneric("coefficients<-"))
setGeneric("formula<-",function(.Object,value) standardGeneric("formula<-"))
