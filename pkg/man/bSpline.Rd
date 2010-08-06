\name{bSpline}
\alias{bSpline}

\title{B-Splines}
\description{
  \command{bSpline} computes the B-spline basis matrix.
}
\usage{
bSpline(x, knots, ..., trunc=NULL, whichColumns = NULL)
}

\arguments{
  \item{x}{The predictor variable.}
  \item{knots}{The knots or breakpoints that define the spline basis.}
  \item{trunc}{The spline basis is truncated to 0 for x <=
  trunc. Default truncation below 0. Set equal to NULL to avoid truncation.}
  \item{whichColumns}{Optional selection of of subset of the basis functions}
  \item{\dots}{Additional arguments for \command{splineDesign}, for
    instance \code{ord} which gives the order of the spline function,
    see \command{\link{splineDesign}}}

}
\details{
  \command{bSpline} is a convenience wrapper function for
  \command{splineDesign} as an alternative to \command{bs}.
  
}
\value{
 A matrix of basis function evaluations.
}
\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{\command{\link{splineDesign}}, \command{\link{bs}}}
\examples{
bSpline(1:10,2:8)
bSpline(1:10,2:8,whichColumns=2:3)
}
\keyword{math}
