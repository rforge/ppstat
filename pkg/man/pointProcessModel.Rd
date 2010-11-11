\name{pointProcessModel}
\alias{pointProcessModel}
\alias{pointProcessSmooth}
\alias{ppmFit}
\alias{ppmFit,PointProcessModel-method}

\title{
  Construction and fitting of a generalized linear point process model.
}
\description{
  
  The functions \code{pointProcessModel} and \code{pointProcessSmooth}
  construct objects of class \code{PointProcessModel} and
  \code{PointProcessSmooth}, respectively, and fit a
  generalized linear point process model where the intensity function is
  given in terms of a linear predictor process; a linear combination of
  linear filters. The linear filters are specified through a symbolic
  description of the linear predictor process. The function
  \code{ppmFit} does that actual fitting. 
  
}
\usage{

pointProcessModel(formula, data, family, support, N = 200, Delta,
  basisPoints, Omega, coefficients, fixedCoefficients = list(),
  modelMatrix = TRUE, fit = modelMatrix, varMethod = 'Fisher',
  basisEnv, ...)

pointProcessSmooth(formula, data, family, support, lambda = 1,
  allKnots = FALSE, N = 200, Delta,  basisPoints, coefficients,
  fixedCoefficients = list(), fit = TRUE, varMethod = 'Fisher',
  basisEnv, ...)

ppmFit(model, control=list(), ...)
}
\arguments{

  \item{formula}{an object of class \code{\link{formula}}. A symbolic
    description of the model to be fitted. See \sQuote{Details}.}

  \item{data}{an object of class \code{\linkS4class{MarkedPointProcess}} containing
    the point process data as well as any continuous process data.}

  \item{family}{an object of class
  \code{\linkS4class{Family}}. Specification of the general model
  family containing the specification of the \eqn{\phi}{phi} function,
  which links the linear predictor process to the predictable intensity process.}
  
  \item{support}{a \code{numeric} vector. Specifies the support of the
    filter functions as the interval from \code{support[1]} to
    \code{support[2]}. If \code{support} is of length 1 the support is the
    interval from 0 to \code{support[1]}. This argument must be
    specified if \code{basisPoints} is not.}

  \item{lambda}{a \code{numeric}. The tuning parameter for the
  penalization term used in \code{pointProcessSmooth}. Default value is
  1.}

\item{allKnots}{a \code{logical}. Should all knots for the spline basis expansion be used in
  \code{pointProcessSmooth}. If not, the same thinning algorithm is used as in \code{smoothing.spline}. Default value is
  FALSE.}

  \item{N}{a \code{numeric}. The number of basis function evaluations
    used in the support. Default value 200.}

  \item{Delta}{a \code{numeric}. Basis functions are evaluated at
    \code{Delta}-grid values in the support. If missing, \code{Delta} is
    set to the length of the support divided by \code{N}. If specified,
    overrides the use of \code{N}.}
  
  \item{basisPoints}{a \code{numeric} vector. Specifies the points where
    the basis functions are evaluated. If missing, the basis functions
    are evaluated equidistantly at \code{Delta}-grid points in the
    support.}
  
  \item{Omega}{an optional square matrix resulting in a quadratic penalization of
    the parameter vector. If missing, no penalization term is added to
    the minus-log-likelihood when the estimates are computed.}
  
  \item{coefficients}{an optional specification of the initial parameters
    used for the numerical optimization.}
  
  \item{fixedCoefficients}{a \code{list} describing parameters that
    should be fixed.  The specification is given in terms of two entries
    named \code{which} and \code{values}, which give the position in the
    parameter vector and the value, respectively, of the parameters that
    should be fixed.}
  
  \item{modelMatrix}{a \code{logical}. Specifies if the model matrix is
    to be computed. Default is \code{TRUE}.}

  \item{fit}{a \code{logical}. Specifies if the model is to be
    fitted. Default value is that the model is only fitted if the model
    matrix is computed. If the \code{formula} does not have a response
    variable the model is not fitted - disregarding the value of \code{fit}.}
  
  \item{varMethod}{a \code{character}. Specifies the method used for
    estimation of the variance matrix of the parameter
    estimators. Currently the default value, \code{'Fisher'}, and
    \code{'none'} are implemented, with \code{'Fisher'} meaning that the
    inverse of the observed Fisher information is used and \code{'none'}
    meaning that the variance matrix is not estimated.}

  \item{basisEnv}{an optional \code{environment} holding precomputed
    basis function evaluations.}

  \item{model}{a point process model as an object of class
    \code{\linkS4class{PointProcessModel}}.}

  \item{control}{a \code{list} of control parameters for \code{optim}.}

  \item{\dots}{additional parameters that are passed on to \code{optim}.}
}
\details{
 \code{pointProcessModel} is the main function for estimation of a generalized
  linear point process model. The function sets up the model 
  specified in terms of a formula, and, if required, computes the  
  model matrix and fits the model to the data.

  \code{pointProcessSmooth} is a wrapper to \code{pointProcessModel}
  that treats terms in the formula of the form \code{s(.)} in a special
  way. Such terms are replaced by an automatic basis expansion in terms
  of B-splines, and the corresponding parameters will be penalized in
  the standard way to ensure smoothness. 

  When calling \code{pointProcessModel} the three arguments
  \code{formula}, \code{data} and \code{family} must always be
  specified. In addition, either the argument \code{support} (and
  optionally \code{N} or \code{Delta}) or the argument
  \code{basisPoints} needs to be specified.

  The \code{formula} specifies together with the \code{family} the
  conditional intensity for a point process in terms of the history of
  the point process itself as well as additional processes and
  variables. The \code{formula} alone specifies a \dfn{linear predictor
  process}, which consists of a linear combination of stochastic
  processes determined by the terms in the formula. The parameters to be
  estimated in the model are the coefficients in this linear
  combination. The stochastic processes are computed and stored in the
  model matrix. The \code{family} specifies a \eqn{\phi}{phi} function
  that transforms the linear predictor process into the conditional
  intensity. 

  The possibilities for specifying the linear predictor process in terms of the formula are
  quite extensive. For instance, any variable corresponding to
  continuous process observations, unit specific variables, the unit
  identification variable or the position/time variable can enter in the
  formula specification by name just as if we did an ordinary regression
  on these variables.  The variables enter directly, as specified by the
  formula, into the intensity through the inverse link function of the
  linear predictor process.  Interactions can be specified, e.g. if
  there is a column named \code{V} and the unit identification variable
  is \code{id}, the inclusion of \code{id*V} in the formula
  specification will result in unit specific parameters for the
  effect of \code{V} on the intensity. Standard rules for formula
  specification applies, e.g. \code{id*V} expands to \code{id + V +
  id:V}, and \code{id:V} by itself results in a different
  parametrization. In addition, variables can be transformed using either
  a single function or basis expansions. A model formula specification
  of the intensity that only relies on the variables listed above
  results in a (conditional) inhomogeneous Poisson process model.

  Other possible terms in the formula specification are terms involving
  the mark variables, that is, the point process data. The response (the left hand side of the formula)
  must always be one of the mark variables. If the response is
  \emph{not} specified, the model is set up, and the model matrix is
  computed - if requested - but the model is not fitted. 

  How a point process, as given by a term in the formula involving
  a mark variable, enters in the model is partly determined by the family.
  For the \code{Hawkes}-family, point processes enter in the model
  through linear filters specified by filter functions. That is, if points for
  the mark variable \code{M} are observed at \eqn{s_1, \ldots s_n}{s_1, ...,
    s_n}, and \code{f(M)} is a term in the formula, the corresponding
  linear filter 

  \deqn{\sum_{i : t-s_i \in \mathrm{support}} f(t-s_i)}{sum_{i: t-s_i \in support} f(t-s_i)}

  then enters as a term in the linear predictor process.  Expansions of
  a filter function in terms of basis functions is possible, using
  e.g. \code{bSpline} to generate an expansion in terms of
  B-splines. See also \code{\linkS4class{Family}} for details on the
  family and \eqn{\phi}{phi} function specification. The term \code{M}
  results in a linear filter function on the support, whereas the term
  \code{const(M)} results in a constant function on the support. Thus
  an affine filter function can be obtained by including \code{I(const(M)
  + M)} in the formula.

  As described above, there is a fundamental difference between the ways 
  a point process and a continuous process enter in the model
  specification. A continuous process enters \emph{as is}
  whereas a point process enters through a linear filter. This is
  believed to be natural in most cases - the point history needs to be
  aggregated through a linear filter before it enters in the
  specification of the conditional intensity, whereas it is the
  instantaneous value of the continuous process that determines the
  conditional intensity. However, it may be desirable to apply a linear
  filter to a continuous process as well. This is possible by including
  terms involving a continuous process variable name concatenated with a
  ".d". That is, if \code{V} is a continuous process variable and
  \code{f(V.d)} is a term in the formula, the corresponding
  linear filter 

  \deqn{V_0 + \int_{\mathrm{support}} f(s) (V_{(t-s)-} - V_0) \mathrm{d}s
  }{V(0) + int_{support} f(s) (V((t-s)-) - V(0))ds}

  then enters as a term in the linear predictor process. As for point
  processes, the filter function can be expanded using any set of basis
  functions. 

  \code{ppmFit} does the actual estimation of the model. It is usually
  not called directly by the user. The current implementation does
  numerical minimization of the minus-log-likelihood function using
  \code{optim} - in general using the BFGS-algorithm. However, if the identity
   \eqn{\phi}{phi} function is used (choose \code{link = "identity"} in the family
  specification) the L-BFGS-B-algorithm is used
  instead, with the default lower bound \code{sqrt(Machine$double.eps)} on the
  intercept parameter and 0 on other parameters.  
  
  Control parameters for the numerical optimization may be passed to the
  \code{optim} function via the \code{control} argument for
  \code{ppmFit}. The maximal number of iterations for the
  algorithm is set to 1000 when used in \code{ppmFit} in
  contrast to the default for \code{optim}, which is 100. See
  \code{\link{optim}} for details. A \code{control}
  argument in a call of \code{pointProcessModel} is automatically passed
  on to the call of \code{ppmFit}.

  The \code{optim} method can always be specified via a named
  \code{method} argument, which is passed to \code{optim}. The L-BFGS-B
  algorithm allows for box-constraints on the parameters. This may be
  required, e.g. if \eqn{\phi}{phi} is the identity, as mentioned
  above. The lower and upper bounds are specified in the call by setting
  the additional arguments \code{lower} and \code{upper} that are passed
  to the \code{optim} function.  See \code{\link{optim}} for details.
  Note that for estimates computed using box-constraints with one or
  more parameter value(s) on the boundary of the parameter space, the
  estimated standard errors based on the inverse Fisher information are
  unreliable. This holds true for the default method with the identity
  \eqn{\phi}{phi} function as well.}

\value{ An object of class \code{\linkS4class{PointProcessModel}}
  containing the estimated parameters.

  The function \code{pointProcessSmooth} returns an object of class
  \code{\linkS4class{PointProcessSmooth}}, which is an extension of
  \code{PointProcessModel}, but currently with an identical content. 
  
  Use \code{summary} to get a standard summary output from the
  fitted model.

  Use \code{termPlot} on the resulting model to plot  the estimated
  filter function(s).
}
\note{

  The functionality of  \code{pointProcessSmooth} is limited. It has one
  standard scheme for the computation of the B-spline basis functions,
  and there is currently no automatic selection of the tuning
  parameter. Moreover, the various functions inherited from class
  \code{PointProcessModel} that can be applied to an
  object of class \code{PointProcessSmooth} may produce misleading
  information. There is e.g. no correction of the degrees of freedom,
  nor of the AIC statistics owing to the fact that this is a smoothed
  and penalized fit.

  The current implementation supports linear filters of continuous
  processes through the \code{.d} notation in the formula. However, the
  data structure allows for non-equidistant sampling, and the filter
  computations are done using numeric integration and not
  FFT. Consequently, the computation of the model matrix for large data
  sets can be prohibitively slow.

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
  \code{\link{optim}}, \code{\link{bSpline}}, \code{\link{termPlot}}.
}
\examples{
data(toyData)
toyPPM <- pointProcessModel(BETA ~ const(ALPHA),
                            data = toyData,
                            family = Hawkes(),
                            support = 2)
summary(toyPPM)
\dontrun{## Plot of the filter function
termPlot(toyPPM)}

toyPPM <- pointProcessModel(BETA ~ cut(ALPHA, c(0,1,2)),
                            data = toyData,
                            family = Hawkes(),
                            support = 2)
summary(toyPPM)
\dontrun{## Plot of the filter function
termPlot(toyPPM)}

## Removing the baseline intensity.
toyPPM <- update(toyPPM, .~. -1)
summary(toyPPM)

## Different baselines for each unit. This 'update' adds a term
## not originally in the model, thus the model matrix will be
## recomputed. 
toyPPM <- update(toyPPM, .~. + id)
summary(toyPPM)


}
\keyword{models}
