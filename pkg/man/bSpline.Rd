\name{bSpline}
\alias{bSpline}

\title{B-Splines}
\description{
  \code{bSpline} computes the B-spline basis matrix.
}
\usage{
bSpline(x, knots, ..., sym = FALSE, trunc = NULL)
}

\arguments{
  
  \item{x}{a \code{numeric} vector. The predictor variable.}
  
  \item{knots}{a \code{numeric} vector. The knots or breakpoints that
    define the spline basis.}

  \item{sym}{a \code{logical}. Should the basis functions be symmetric
  around 0. Default value is \code{FALSE}.}
  
  \item{trunc}{a \code{numeric} of length 1 or 2. The spline basis is
    truncated below \code{trunc} or to the interval from \code{trunc[1]}
    to \code{trunc[2]} if \code{trunc} has length 2. Default value is
    NULL, which means no truncation.}
  
  \item{\dots}{additional arguments passed to \code{splineDesign}, for
    instance, \code{ord} which gives the order of the spline functions,
    see \code{\link{splineDesign}}.}
  
}
\details{
  
  \code{bSpline} is a convenience wrapper function for
  \code{splineDesign} as an alternative to \code{bs}. It has defaults
  and some useful parameters suitable for formula specification of
  linear filter function bases when using
  \code{\link{pointProcessModel}}.
  
}
\value{
  
  A matrix of basis function evaluations. Zero-columns returned from
  \code{splineDesign} are discarded.

} \author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{\code{\link{splineDesign}}, \code{\link{bs}}, \code{\link{pointProcessModel}}}
\examples{
bSpline(seq(1,10), 2:8)
bSpline(seq(-10,10), 2:8, sym = TRUE)
bSpline(seq(-10,10), seq(-10,10,2), trunc = c(-4,4))
}
\keyword{math}
