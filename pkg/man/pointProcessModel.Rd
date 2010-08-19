\name{pointProcessModel}
\alias{pointProcessModel}
\alias{glppmFit}

\title{
  Construction and fitting of a generalized linear point process model.
}
\description{
  This function constructs an object of class \code{PointProcessModel}
  and fits a generalized linear point process model
  where the intensity function is given in terms of a linear predictor
  process; a linear combination of linear filters. The linear filters
  are specified through a symbolic description of the linear predictor
  process.  
}
\usage{
pointProcessModel(formula, data, family, support, Delta, basisPoints,
  Omega, coefficients, fixedCoefficients = list(), fit = TRUE,
  modelMatrix = TRUE, varMethod = "Fisher", ...)

glppmFit(model, control=list(), ...)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}}: a symbolic description of the model to be
    fitted. The details of model specification are given under
    \sQuote{Details}.}
  \item{data}{an object of class \code{\linkS4class{ProcessData}} containing
    the point process data as well as any continuous process data.}
  \item{family}{specification of the \eqn{\phi}{phi} function.}
  \item{support}{a \code{numeric} vector. Specifies the
  support of the filter functions as the interval from \code{support[1]}
  to \code{support[2]}. If \code{support} is of length 1 the support is
  the interval from 0 to \code{support[1]}. }
  \item{Delta}{a \code{numeric}. Basis functions are evaluated at
    \code{Delta}-grid values in the support.}
  \item{basisPoints}{A \code{numeric} vector. Specifies the points where
  the basis functions are evaluated. Either this argument must be
  specified or \code{support} together with \code{Delta}.}
  \item{Omega}{a square matrix resulting in a quadratic penalization on
    the parameter vector. Default is \code{NULL} meaning no penalization.}
  \item{coefficients}{an optional specification of the initial parameters
    used for the numerical optimization.}
  \item{fixedCoefficients}{a list describing parameters that should be fixed in
    the form of a list containing two entries named \code{which} and
    \code{values} that give the position in the parameter vector and the
    value, respectively, of the parameters that should be fixed.}
  \item{fit}{a \code{logical}. Specifies if the model is to be fitted.}
  \item{modelMatrix}{a \code{logical}. Specifies if the model matrix is
    to be computed.}
  \item{model}{a point process model as an object of class
    \code{\linkS4class{PointProcessModel}}}
  \item{control}{a list of control parameters for \command{optim}}
  \item{\dots}{additional parameters that will be passed to \command{optim}}
}
\details{
 \command{pointProcessModel} is the main function for estimation of a generalized
  linear point process model. This function sets up the model when
  specified in terms of a formula and computes the required point
  process model of class \code{\linkS4class{PointProcessModel}}.
  Numerical minimization of the minus-log-likelihood function relies on
  \command{optim} with the L-BFGS-B-algorithm.

  A column stored in the \code{value}-slot in the
  \code{continuousProcess} can enter in the formula specification by
  name just as if we did an ordinary regression on the variable stored
  in the column. So can the \code{Id} and \code{position}
  slot. The variables from the \code{continuousProcess} enter directly,
  as specified by the formula, into the intensity through
  the inverse link function of the linear predictor process. 
  Interactions can also be specified, e.g. if there is a column
  named \code{factorX} the inclusion of \code{Id*factorX} in the
  formula specification will result in an \code{Id}-specific 
  parameter for the effect of \code{factorX} on the
  intensity. Specifications of the intensity that only relies on 
  variables from the \code{continuousProcess} result in an
  inhomogeneous Poisson process model.

  Other possible terms in the formula specification are terms involving
  one or more of the \code{markTypes} in the
  \code{markedPointProcess}. The response must be one of the
  \code{markTypes}.

  For the \code{Hawkes}-family (currently only one implemented) the
  point process(es) enter in the model specification through linear
  filters specified by functions. Expansions of the linear filter
  functions in terms of basis functions is possible, using
  e.g. \code{bSpline} to generate an expansion in terms of
  B-splines. See \code{\link{Family}} for details on the family and
  \eqn{\phi}{phi} function specification.

  \command{glppmFit} does the actual estimation of the model. It is usually not
  called explicitly by the user.
  
  Control parameters for the numerical optimization may be passed to the
  \command{optim} function via the \code{control} argument for
  \command{glppmFit}. The maximal number of iterations for the
  L-BFGS-B-algorithm is set to 1000 when used in \command{glppmFit} in
  contrast to the default for \command{optim}, which is 100. See
  \code{\link{optim}} for details. A \code{control}
  argument in a call of \code{pointProcessModel} is automatically passed
  on to the call of \command{glppmFit}.

  The L-BFGS-B algorithm allows for box-constraints on the parameters. A
  typical constraint is positivity of all parameters. This is required
  if the \eqn{\phi}{phi} function is the identity (choose
  \code{link="identity"} in the Family specification). The lower and
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
  generalized linear point processes}. 	arXiv:1003.0848v1 
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
