\name{update-methods}
\alias{update,PointProcessModel-method}
\alias{update}
\title{update}
\description{
  
  The function \command{update} will update and re-fit a
  \code{PointProcessModel}.
  
}
\usage{

\S4method{update}{PointProcessModel}(object,formula,fisherInformation=TRUE,warmStart=TRUE,fixedPar=list(),...)

}
\arguments{
  \item{object}{The \code{PointProcessModel} object to be updated.}
  \item{formula}{Changes to the formula - see 'update.formula' for
    details.}
  \item{fisherInformation}{logical. Should the observed Fisher
  information and its inverse be computed. Default is \code{TRUE}, but
  the computations can be time consuming.}
  \item{warmStart}{logical. Should the optimization use the previously
    estimated values of the parameters in the larger model. Default is
    \code{TRUE}.}
  \item{fixedPar}{a list describing parameters that should be fixed in
  the form of a list containing two entries named \code{which} and
  \code{values} that give the position in the parameter vector and the
  value, respectively, of the parameters that should be fixed.}
  \item{...}{Additional arguments to the call}
}

\value{
  
  An object of class \code{\linkS4class{PointProcessModel}} containing the
  estimated parameters.
  
}

\details{

  This function can be used to update and refit an object of class
  \code{\linkS4class{PointProcessModel}} as an alternative to making a
  new object by calling the function \command{glppm} directly.

  If possible, the updated model will try to reuse (parts of) the model
  matrix and thus not recompute the entire model matrix as would
  otherwise happen by a direct call to \command{glppm}. This can reduce
  computations if the model uses many linear filter specifications as
  the computation of the model matrix is then often far more demanding
  than the subsequent optimization.

  Note that the current implementation of \code{update} is in general not able to
  identify the correct columns in the model matrix if the original model
  includes interaction terms and will most likely result in a
  recomputation of the model matrix. If interactions are included and
  the model matrix is not recomputed a warning is then issued.   
}


\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{

\code{\link{update}}, \code{\link{update.formula}}, \code{\link{glppm}}, \code{\linkS4class{PointProcessModel}}

}

\keyword{models}
