\name{computeBasis}
\alias{computeBasis}
\alias{computeBasis-methods}
\alias{computeBasis,PointProcessModel,character,character-method}
\title{
Computes function evaluations
}
\description{

  This function computes function evaluations and stores them in a
  given environment.  
  
}
\usage{
\S4method{computeBasis}{PointProcessModel,character,character}
}
\arguments{
  \item{model}{an object of class \code{"PointProcessModel"}.}

  \item{term}{a \code{"character"}. Usually a term obtained from a
    formula, which gives the function to be evaluated.
    Coerced into the formula \code{~ term - 1}.}

  \item{varLabels}{an optional \code{"character"} vector. Holds the valid variable
    names for variables that enter in the formula. The basis evaluations
    are only computed and stored if the variables in the term are in
    this is vector. The check is bypassed if the argument is missing.}
}
\details{

  The function is called when building the model matrix for an object of
  class \code{"PointProcessModel"}.

  The function checks if the variables in the term are valid, that is,
  if they are among the variable names in \code{varLabels}. Then it
  checks if the environment in \code{model} holding the basis
  evaluations is locked, and if it is, if the basis evaluations for the
  term are present. If so, the function returns \code{TRUE}. It is an
  error if the term is not present (should not happen, an error message
  is given). Otherwise the basis evaluations are computed and stored
  in the appropriate environment in the model. 
  
}
\value{
A \code{"logical"}, which is \code{TRUE} if the variables are in the
\code{varLabels} and \code{FALSE} otherwise.

This function is called for its side effect of computing and storing the basis
evaluations in the environment. If the environment is locked the
evaluations are not recomputed.

}

\author{
Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
}

\seealso{\code{\linkS4class{"PointProcessModel"}}
}

\keyword{math}
