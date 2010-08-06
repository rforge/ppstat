\name{const}
\Rdversion{1.1}
\alias{const}

\title{
Indicator function multiplied by a constant.
}
\description{

  This function is a simple wrapper intended to be used in a
  formula specification for generalized linear point process models.

}
\usage{
const(x, y, c = 1)
}

\arguments{
  \item{x}{
  a numeric.  
  }
  \item{y}{
  a numeric.
}
  \item{c}{
  a numeric.  
  }
}

\value{
The function returns a numeric.   
The result is \code{c} if \code{x <= y} and 0 otherwise. If \code{x} and
\code{y} are vectors of equal length the result is a corresponding
vector where the entries in \code{x} and \code{y} are compared. 
}

\author{
  Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}
}


\examples{

const(10,5)
const(5,10)
const(1:10,rep(5,10),2)

}
\keyword{ math }
