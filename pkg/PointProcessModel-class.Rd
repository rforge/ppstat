\name{PointProcessModel-class}
\docType{class}
\alias{PointProcessModel-class}
\alias{coefficients,PointProcessModel-method}
\alias{coefficients<-,PointProcessModel-method}
\alias{coefficients<-}
\alias{computeDDMinusLogLikelihood,PointProcessModel-method}
\alias{computeDDMinusLogLikelihood}
\alias{computeDMinusLogLikelihood,PointProcessModel-method}
\alias{computeDMinusLogLikelihood}
\alias{computeLinearPredictor,PointProcessModel-method}
\alias{computeLinearPredictor}
\alias{computeMinusLogLikelihood,PointProcessModel-method}
\alias{computeMinusLogLikelihood}
\alias{computeModelMatrix,PointProcessModel-method}
\alias{computeModelMatrix}
\alias{formula,PointProcessModel-method}
\alias{formula<-,PointProcessModel-method}
\alias{formula<-}
\alias{vcov<-}
\alias{vcov<-,PointProcessModel-method}
\alias{vcov,PointProcessModel-method}
\alias{getModelMatrix,PointProcessModel-method}
\alias{getModelMatrix}
\alias{getLinearFilter,numeric,PointProcessModel-method}
\alias{getLinearFilter}
\alias{initialize,PointProcessModel-method}
\alias{setCoefficients,PointProcessModel-method}
\alias{setFormula,PointProcessModel-method}
\alias{setResponse,PointProcessModel-method}
\alias{show,PointProcessModel-method}
\alias{subset,PointProcessModel-method}

\title{Class "PointProcessModel"}
\description{
  
  An object of class \code{PointProcessModel} is returned from
  \code{glppm} and contains the fitted generalized linear point process
  model.
  
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PointProcessModel", processData, formula, family, coefficients, modelMatrix)}.

The object holds the data used for fitting the model and a formula
specification of the model. When the object is the output of
\code{glppm} the object contains the fitted model parameters (via MLE) and 
the estimated covariance matrix for the parameters (the inverse of the
observed Fisher information).

The model matrix is computed from the data based on the formula when the
object is created. Each row in the matrix corresponds to one observation
point in the \code{position} slot for the continuous process object 
in \code{processData}. 

}
\section{Slots}{
	 \describe{
    \item{\code{processData}:}{Object of class \code{"ProcessData"} -
      the model data.}
    \item{\code{modelMatrix}:}{Object of class \code{"Matrix"} - the
      model matrix.}
    \item{\code{delta}:}{Object of class \code{"numeric"} - the
      interdistance between the observation (time) points.}
    \item{\code{formula}:}{Object of class \code{"formula"} - the model
      formula specifying the model.}
    \item{\code{family}:}{Object of class \code{"Family"} - the family
      object specifying the specific interpretation of the formula
      specification and the link function used.}
    \item{\code{fixedCoefficients}:}{List with two entries, \code{which}
      and \code{value}.}
     \item{\code{penalization}:}{Logical - should penalization be used
       for estimation.}
    \item{\code{Omega}:}{Object of class \code{"matrix"} - the
      positive definite penalization matrix.}
    \item{\code{optimResult}:}{List with results from numerical
      optimization (\code{optim}).}
    \item{\code{coefficients}:}{Object of class \code{"numeric"} - the
      estimated model parameters.}
    \item{\code{var}:}{Object of class \code{"numeric"} - the estimated
      covariance matrix of the estimated model parameters.}
    \item{\code{call}:}{Object of class \code{"character"} - the call to
      glppm that created the point process model.}
  }
}
\section{Methods}{
  \describe{
    \item{coefficients}{\code{signature(object = "PointProcessModel")}:
      Returns the estimated model paramters.}
    \item{coefficients<-}{\code{signature(.Object =
	"PointProcessModel")}: Sets the model parameters.}
    \item{computeDMinusLogLikelihood}{\code{signature(object =
	"PointProcessModel")}: Computes the first derivate of the
      minus-log-likelihood function.}
    \item{computeLinearPredictor}{\code{signature(object =
	"PointProcessModel")}: Computes the linear predictor.}
    \item{computeMinusLogLikelihood}{\code{signature(object =
	"PointProcessModel")}: Computes the minus-log-likelihood function.}
  \item{computeDMinusLogLikelihood}{\code{signature(object =
	"PointProcessModel")}: Computes the derivative minus-log-likelihood function.}
  \item{computeDDMinusLogLikelihood}{\code{signature(object =
	"PointProcessModel")}: Computes the hessian of the minus-log-likelihood function.}
    \item{computeModelMatrix}{\code{signature(object =
	"PointProcessModel")}: Computes the model matrix.}
    \item{formula}{\code{signature(x = "PointProcessModel")}: Returns
      the model formula.}
    \item{formula<-}{\code{signature(.Object = "PointProcessModel")}:
      Sets the model formula.}
    \item{getModelMatrix}{\code{signature(object =
	"PointProcessModel")}: Returns the model matrix.}
    \item{initialize}{\code{signature(.Object = "PointProcessModel")}.}
    \item{setCoefficients}{\code{signature(object =
	"PointProcessModel")}. }
    \item{setFormula}{\code{signature(object = "PointProcessModel")}: ... }
    \item{setResponse}{\code{signature(object = "PointProcessModel")}: ... }
    \item{subset}{\code{signature(x = "PointProcessModel")}: Creates a
      new point process model based on a subset of the data by passing 
the argument(s) to the subset method for the process data, see
\code{\linkS4class{ProcessData}}. The \code{family} and \code{formula}
are retained in the new model. }
    \item{update}{\code{signature(object = "PointProcessModel")}: See \code{\link{update}} }
    \item{vcov}{\code{signature(object = "PointProcessModel")}: Returns the variance-covariance matrix estimated as the
      inverse Fisher information}
    \item{vcov<-}{\code{signature(.Object = "PointProcessModel")}: Sets
      the variance-covariance matrix}
     
	 }
}

\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
	\code{\linkS4class{ProcessData}}, \code{\link{formula}}, \code{\link{glppm}} 
}
\examples{
showClass("PointProcessModel")
}

\keyword{classes}
