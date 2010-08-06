\name{ppstat-package}
\alias{ppstat-package}
\alias{ppstat}
\docType{package}
\title{ Point Process Statistics
}
\description{

This package implements statistical methods for one-dimensional marked
point process model.

}
\details{
\tabular{ll}{
Package: \tab ppstat\cr
Type: \tab Package\cr
Version: \tab 0.5 \cr
Date: \tab 2010-07-01 \cr
License: \tab GPL version 2 or newer \cr
LazyLoad: \tab yes\cr
}

\code{\link{glppm}} is the main function for fitting a generalized
linear point process model to a dataset, which is stored as an object
of class \code{\linkS4class{ProcessData}}. A ProcessData object contains two
components, one object of class \code{\linkS4class{MarkedPointProcess}}
containing the data of the modeled marked point process as well
as any additional covariate marked point processes, and one
object of class \code{\linkS4class{ContinuousProcess}} containing 
potential continuous time covariate processes. The latter object
will always hold the data as a series of observations at
discrete (time) points that are not necessarily equidistant.

The data structures for continuous and marked point processes are
implemented in the package \code{\link{processdata}}.

}
\author{
Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}, with contributions from 
	Lisbeth Carstensen \email{lisbeth@binf.ku.dk}.

Maintainer: \email{Niels.R.Hansen@math.ku.dk}
}
\references{

Hansen, N. R. \emph{Penalized maximum likelihood estimation for generalized linear point processes .} Working paper.

Hansen, N. R. and Carstensen, L. \emph{Modeling genomic organization with 
generalized linear point process models}. Working paper.

}

\examples{

### glppm(...)

}
\keyword{models}
