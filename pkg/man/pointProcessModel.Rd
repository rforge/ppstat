\name{pointProcessModel}
\alias{pointProcessModel}
\alias{pointProcessSmooth}
\alias{ppmFit}

\title{
  Construction and fitting of a generalized linear point process model.
}
\description{
  
  The functions \code{pointProcessModel} and \code{pointProcessSmooth}
  construct an object of class \code{PointProcessModel} and fits a
  generalized linear point process model where the intensity function is
  given in terms of a linear predictor process; a linear combination of
  linear filters. The linear filters are specified through a symbolic
  description of the linear predictor process. The function
  \code{ppmFit} does that actual fitting. 
  
}
\usage{

pointProcessModel(formula, data, family, support, Delta, basisPoints,
  Omega, coefficients, fixedCoefficients = list(), modelMatrix = TRUE,
  fit = modelMatrix, varMethod = 'Fisher', ...)

pointProcessSmooth(formula, data, family, support, Delta, basisPoints,
  coefficients, fixedCoefficients = list(), fit = TRUE,
  modelMatrix = TRUE, varMethod = 'Fisher', ...)

ppmFit(model, control=list(), ...)
}
\arguments{

  \item{formula}{an object of class \code{\link{formula}}: a symbolic
    description of the model to be fitted. The details of model
    specification are given under \sQuote{Details}.}

  \item{data}{an object of class \code{\linkS4class{ProcessData}} containing
    the point process data as well as any continuous process data.}

  \item{family}{an object of class
  \code{\linkS4class{Family}}. Specification of the general model
  family containing the specification of the \eqn{\phi}{phi} function,
  which links the linear predictor process to the predictable intensity process.}
  
  \item{support}{a \code{numeric} vector. Specifies the support of the
    filter functions as the interval from \code{support[1]} to
    \code{support[2]}. If \code{support} is of length 1 the support is the
    interval from 0 to \code{support[1]}. }

  \item{Delta}{a \code{numeric}. Basis functions are evaluated at
    \code{Delta}-grid values in the support.}
  
  \item{basisPoints}{A \code{numeric} vector. Specifies the points where
    the basis functions are evaluated. Either this argument must be
    specified or \code{support} together with \code{Delta}.}
  
  \item{Omega}{a square matrix resulting in a quadratic penalization on
    the parameter vector. Default is \code{NULL} meaning no penalization.}
  
  \item{coefficients}{an optional specification of the initial parameters
    used for the numerical optimization.}
  
  \item{fixedCoefficients}{a list describing parameters that should be
    fixed.  The specification is given in terms of a list containing two
    entries named \code{which} and \code{values} that give the position
    in the parameter vector and the value, respectively, of the
    parameters that should be fixed.}
  
  \item{modelMatrix}{a \code{logical}. Specifies if the model matrix is
    to be computed. Default is \code{TRUE}.}

  \item{fit}{a \code{logical}. Specifies if the model is to be
    fitted. Default value is that the model is only fitted if the model
    matrix is computed.}
  
  \item{varMethod}{a \code{character}. Specifies the method used for
    estimation of the variance matrix of the parameter
    estimators. Currently the default value, \code{'Fisher'}, and
    \code{'none'} are implemented, with \code{'Fisher'} using the
    inverse of the observed Fisher information and \code{'none'}
    preventing the estimation of the variance matrix.}

  \item{basisEnv}{an optional \code{environment} holding precomputed
  basis function evaluations.}

  \item{model}{a point process model as an object of class
    \code{\linkS4class{PointProcessModel}}.}

  \item{control}{a list of control parameters for \code{optim}.}

  \item{\dots}{additional parameters that are be passed on to \code{optim}.}
}
\details{
 \code{pointProcessModel} is the main function for estimation of a generalized
  linear point process model. The function sets up the model 
  specified in terms of a formula, and, if required, computes the  
  model matrix and fits the model to the data.

  When calling \code{pointProcessModel} the three arguments
  \code{formula}, \code{data} and \code{family} must always be
  specified. In addition, either the two arguments \code{support} and
  \code{Delta} or the argument \code{basisPoints} must be specified. 

  A column stored in the \code{value}-slot in the
  \code{continuousProcess} can enter in the formula specification by
  name just as if we did an ordinary regression on the variable stored
  in the column. So can the \code{id} and \code{position}
  slot. The variables from the \code{continuousProcess} enter directly,
  as specified by the formula, into the intensity through
  the inverse link function of the linear predictor process. 
  Interactions can also be specified, e.g. if there is a column
  named \code{factorX} the inclusion of \code{id*factorX} in the
  formula specification will result in an \code{id}-specific 
  parameter for the effect of \code{factorX} on the
  intensity. Specifications of the intensity that only relies on 
  variables from the \code{continuousProcess} result in an
  inhomogeneous Poisson process model.

  Other possible terms in the formula specification are terms involving
  one or more of the \code{markTypes} in the
  \code{markedPointProcess}. The response must be one of the
  \code{markTypes}.

  For the \code{Hawkes}-family the
  point process(es) enter in the model specification through linear
  filters specified by functions. Expansions of the linear filter
  functions in terms of basis functions is possible, using
  e.g. \code{bSpline} to generate an expansion in terms of
  B-splines. See \code{\link{Family}} for details on the family and
  \eqn{\phi}{phi} function specification.

  \code{ppmFit} does the actual estimation of the model. It is usually not
  called explicitly by the user. The current implementation does
  numerical minimization of the minus-log-likelihood function 
  using \code{optim} with the L-BFGS-B-algorithm.
  
  Control parameters for the numerical optimization may be passed to the
  \code{optim} function via the \code{control} argument for
  \code{ppmFit}. The maximal number of iterations for the
  L-BFGS-B-algorithm is set to 1000 when used in \code{ppmFit} in
  contrast to the default for \code{optim}, which is 100. See
  \code{\link{optim}} for details. A \code{control}
  argument in a call of \code{pointProcessModel} is automatically passed
  on to the call of \code{ppmFit}.

  The L-BFGS-B algorithm allows for box-constraints on the parameters. A
  typical constraint is positivity of all parameters. This is required
  if the \eqn{\phi}{phi} function is the identity (choose
  \code{link="identity"} in the family specification). The lower and
  upper bounds are specified in the call by setting the additional arguments
  \code{lower} and \code{upper} that are passed to the \code{optim}
  algorithm.  See \code{\link{optim}} for details.
}
\value{
  An object of class \code{\linkS4class{PointProcessModel}} containing the
  estimated parameters.

  Use \code{summary} to get a standard summary output from the
  fitted model.

  Plotting the model gives a plot of the estimated linear 
  filter function(s).
}
\references{
  Hansen, N. R. \emph{Penalized maximum likelihood estimation for
  generalized linear point processes}. 	arXiv:1003.0848v1. \url{http://arxiv.org/abs/1003.0848} 
}
\author{
Niels Richard Hansen \email{Niels.R.Hansen@math.ku.dk}.

Maintainer: \email{Niels.R.Hansen@math.ku.dk}
}

\seealso{
  \code{\linkS4class{PointProcessModel}},  \code{\link{update}},
  \code{\link{optim}}, \code{bSpline} 
}
\examples{
}
\keyword{models}
