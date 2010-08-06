\name{Family-class}
\docType{class}
\alias{Family-class}
\alias{initialize,Family-method}
\alias{makeLink}

\title{Class "Family"}
\description{	 

An object of class \code{Family} specifies the family and 
corresponding \eqn{\phi}{phi} function for a generalized linear point process
model. 

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Family", name,
  link, phi, Dphi, D2phi)}.


}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"}.}
    \item{\code{link}:}{Object of class \code{"character"}.}
    \item{\code{phi}:}{Object of class \code{"function"}.}
    \item{\code{Dphi}:}{Object of class \code{"function"}.}
    \item{\code{D2phi}:}{Object of class \code{"function"}.}
  }
}		
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Family")}.}
  }   
}
    
\section{Details}{
  One supported family is the \code{Hawkes} family with
  the \code{log}-link function:

  \code{new("Family",name="Hawkes",link="log")}

  Other implemented standard link functions are, "logit", "cloglog", "identity", "root", "logaffine". 
  The latter two requires the additional parameter c, which by default
  is 0. For "root" this gives
  \deqn{\phi(x) = x^{c+1}1(x > 0)}{phi(x) = x^{c+1}1(x > 0)}
  and for "logaffine"
  \deqn{\phi(x) =  \exp(x)1(x < c) + \exp(c)(x - c + 1)1(x > c).}{phi(x) =  exp(x)1(x < c) + exp(c)(x - c + 1)1(x > c).}
  Using "identity" requires a constrained parameter space. 
  
  }
    
\author{Niels Richard Hansen, \email{Niels.R.Hansen@math.ku.dk}}

\seealso{
	 \code{\linkS4class{PointProcessModel}} and  \code{\link{glppm}}. 
}
\examples{
showClass("Family")
}
\keyword{classes}
