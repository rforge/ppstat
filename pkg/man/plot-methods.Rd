\name{plot-methods}
\alias{plot,numeric,PointProcessModel-method}

\title{plot}
\description{
  Plots the linear filter functions in for a \code{PointProcessModel}.
}
\usage{
\S4method{plot}{numeric,PointProcessModel}(x,y,trans=NULL,alpha=0.05,...)
}
\arguments{
  \item{x}{a \code{numeric} vector giving the evaluation points of the
    linear filter(s)}
  \item{y}{an object of class \code{PointProcessModel}.}
  \item{trans}{a function or character string. The function used to transform the linear filter functions
    before plotting}
  \item{alpha}{The 1-confidence level. Default \code{0.05}}
  \item{...}{Additional arguments to the call}

}
\value{

  A \code{\link{ggplot}} object.
  
}
\details{
  
  The (estimated) linear filter functions in the
  \code{PointProcessModel} are extracted, evaluated in the \code{x}
  values and plotted. Approximate pointwise level \code{1-alpha} confidence intervals
  are also plotted.

  Setting \code{trans} to e.g. \code{exp} the linear filters are
  exponentiated before plotting. This is in particular relevant
  for the \code{log}-link function, as this will then provide the
  multiplicative effect of the filter on the intensity as an alternative
  to the \code{log}-additive effect. 
}

\seealso{
  
  \code{\linkS4class{PointProcessModel}}, \code{\link{glppm}}, \code{\link{ggplot}}
  
}

\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}
\keyword{aplot}


