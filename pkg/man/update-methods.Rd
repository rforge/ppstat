\name{update-methods}
\alias{update,PointProcessModel-method}
\alias{update}
\title{update}
\description{
  
  The function \command{update} will update and refit a
  \code{PointProcessModel}.
  
}
\usage{

\S4method{update}{PointProcessModel}(model, formula, warmStart = TRUE,
fixedPar = list(), ...)

}
\arguments{
  \item{model}{the \code{PointProcessModel} object to be updated.}
  \item{formula}{Changes to the formula - see 'update.formula' for
    details.}
  \item{warmStart}{a \code{logica}l. Specifies if the optimization uses the previously
    estimated values of the parameters in the larger model as initial
    values if possible. Default is
    \code{TRUE}.}
  \item{fixedCoefficients}{a list describing parameters that should be fixed in
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
  This function can be used to update and refit a model of class
  \code{\linkS4class{PointProcessModel}} as an alternative to constructing a
  new model object by calling the function \command{pointProcessModel} directly.

  If possible, the updated model will reuse (parts of) the model
  matrix and thus not recompute the entire model matrix as would
  happen by a direct call to \command{pointProcessModel}. Reuse will happen if the updated
  model is a submodel of the original model and will reduce
  computation time as well as memory usage. For large data
  sets and/or many linear predictors the model matrix becomes very
  large, and the time to compute the matrix can be considerable. Most
  notably, by using \code{update} the model matrix is copied by
  reference from the original model to the submodel, and thus does not
  take up memory twice. When investigating several models that are all
  submodels of a single supermodel it is therefore recommended to fit
  the supermodel first and use \code{update} to fit the different
  submodels. 

  The current implementation of \code{update} works by using the assignment of
  columns in the model matrix to terms in the formula. This assignment
  is given by the attribute \code{assign} to the model matrix. A new model
  is identified as a submodel if all terms in the new model are terms in
  the original model. The updated model matrix consists of the columns
  that are assigned to a term in the new model. An important consequence
  is that if the original model includes interactions the updated model
  may not be parametrized as expected or as if the model was specified
  using \command{pointProcessModel}. The updated model may not even be
  as expected. In this case a warning is given.
}


\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{

\code{\link{update}}, \code{\link{update.formula}}, \code{\link{glppm}}, \code{\linkS4class{PointProcessModel}}

}

\keyword{models}
