\name{PointProcessModel-class}
\docType{class}
\alias{PointProcessModel-class}
\alias{PointProcessSmooth-class}
\alias{coefficients,PointProcessModel-method}
\alias{coefficients<-,PointProcessModel-method}
\alias{computeDDMinusLogLikelihood,PointProcessModel-method}
\alias{computeDDMinusLogLikelihood}
\alias{computeDMinusLogLikelihood,PointProcessModel-method}
\alias{computeDMinusLogLikelihood}
\alias{computeLinearPredictor,PointProcessModel-method}
\alias{computeLinearPredictor}
\alias{computeModelMatrix,PointProcessModel-method}
\alias{computeModelMatrix}
\alias{getLinearFilter,PointProcessModel-method}
\alias{getLinearFilter}
\alias{getModelMatrix,PointProcessModel,ANY-method}
\alias{getModelMatrix}
\alias{setModelMatrix<-,PointProcessModel,ANY-method}
\alias{setModelMatrix<-}
\alias{plot,PointProcessModel-method}
\alias{predict,PointProcessModel-method}
\alias{show,PointProcessModel-method}
\alias{subset,PointProcessModel-method}
\alias{vcov,PointProcessModel-method}

\title{Class "PointProcessModel"}
\description{
  
  An object of class \code{"PointProcessModel"} is returned from the
  constructor \code{pointProcessModel}. It contains the fitted
  generalized linear point process model. The class inherits class
  \code{"PointProcess"}. The constructor \code{pointProcessSmooth}
  returns an object of class \code{"PointProcessSmooth"}, which 
  inherits class \code{"PointProcessModel"}.
  
}
\section{Objects from the Class}{
  
  Objects can be created by calls of the form
  \code{pointProcessModel(formula, data, family, support)} where
  \code{data} is an object of class \code{MarkedPointProcess}.
  
  The model object contains the data used for fitting the model and a
  formula specification of the model. When the object is the result of
  calling \code{pointProcessModel} the object contains the fitted model
  parameters (via MLE) and the estimated covariance matrix for the
  parameters.
  
  A model matrix is by default computed from the data based on the
  formula when the object is created using
  \code{pointProcessModel}. Each row in the matrix corresponds to one
  grid point in the data set.

}
\section{Slots}{
  \describe{
    \item{\code{basisEnv}:}{an \code{"environment"} containing a list
      with basis function evaluations.}

    \item{\code{basisPoints}:}{a \code{"numeric"} containing the
    evaluation points for the basis functions.}
    
    \item{\code{coefficients}:}{a \code{"numeric"} vector containing the
      estimated model parameters.} 
    
    \item{\code{fixedCoefficients}:}{list with two entries, \code{which}
      and \code{value}.}

    \item{\code{modelMatrixCol}:}{a \code{numeric} vector.}

    \item{\code{modelMatrixEnv}:}{an environment that contains an
      object of class \code{"Matrix"} - the model matrix.}

    \item{\code{Omega}:}{object of class \code{"matrix"} - the
      positive semi-definite penalization matrix.}
    
    \item{\code{optimResult}:}{list containing the results from numerical
      optimization (\code{optim}).}
    
    \item{\code{penalization}:}{a \code{"logical"} specifying if
      penalization is used for estimation.}
     
    \item{\code{var}:}{a \code{"matrix"}. The estimated
      variance matrix of the model parameter estimator.}

    \item{\code{varMethod}:}{a \code{"character"} specifying the method
      used for estimating the variance matrix.}
    
  }
}
\section{Methods}{
  \describe{
    
    \item{coefficients}{\code{signature(object = "PointProcessModel")}:
      Returns the estimated model parameters.}

    \item{coefficients<-}{\code{signature(.Object =
	"PointProcessModel")}: Sets the model parameters.}

    \item{computeBasis}{\code{signature(model = "PointProcessModel",
	term = "character", varLabels = "character")}: Computes the
	basis function evaluations.}
    
    \item{computeDDMinusLogLikelihood}{\code{signature(model =
	"PointProcessModel")}: Computes the hessian of the
      minus-log-likelihood function.}
    
    \item{computeDMinusLogLikelihood}{\code{signature(model =
	"PointProcessModel")}: Computes the first derivate of the
      minus-log-likelihood function.}

    \item{computeLinearPredictor}{\code{signature(model =
	"PointProcessModel")}: Computes the linear predictor.}
    
    \item{computeModelMatrix}{\code{signature(model =
	"PointProcessModel")}: Computes the model matrix.}

    \item{getModelMatrix}{\code{signature(model =
	"PointProcessModel")}: Returns the model matrix.}

    \item{getLinearFilter}{\code{signature(model =
	"PointProcessModel")}\code{(model, se = FALSE, nr, ...)}: Returns the values
	of the linear filter functions as a \code{data.frame}. If
	\code{se = TRUE}, the functions returns a list with the linear
	filter functions in the first component and the estimated
	standard errors in the second component.}

     \item{setModelMatrix<-}{\code{signature(model =
	 "PointProcessModel",value = "Matrix")}: Sets the model matrix.}

      \item{predict}{\code{signature(object =
	 "PointProcessModel")}: Returns the intensity. This is the
	 linear predictor, as computed by \code{computeLinearPredictor},
	 transformed by \eqn{\phi}{phi}.}
     
     \item{termPlot}{\code{signature(model =
	 "PointProcessModel")}: Plots the (estimated) linear filter
       functions.}
     
     \item{show}{\code{signature(object =
	 "PointProcessModel")}: Prints the model object.}

      \item{summary}{\code{signature(object =
	 "PointProcessModel")}: Summarizes the model object. See
	 \code{\link[=summary,PointProcessModel-method]{summary}}.}
     
    \item{subset}{\code{signature(x = "PointProcessModel")}: Creates a
      new point process model based on a subset of the data by passing
      the argument(s) to the subset method for the process data, see
      \code{\linkS4class{MarkedPointProcess}}. The \code{family} and
      \code{formula} are retained in the new model. }

    \item{update}{\code{signature(model = "PointProcessModel")}: See
      \code{\link{update}}.}
    
    \item{vcov}{\code{signature(model = "PointProcessModel")}: Returns
      the estimated variance matrix of the parameter estimator.}
     
	 }
}

\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{\code{\link{pointProcessModel}},
  \code{\link{pointProcessSmooth}}, \code{\linkS4class{PointProcess}},
  \code{\linkS4class{ProcessData}},
  \code{\link[=summary,PointProcessModel-method]{summary}}, \code{\link{formula}}.
}
\examples{
showClass("PointProcessModel")
}

\keyword{classes}
