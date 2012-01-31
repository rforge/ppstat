\name{NEWS}
\title{ppstat News}

\section{Changes in processdata version 0.8}{
  \subsection{New features}{
    \itemize{
      \item The class \code{MultivariatePointProcess} has been
      introduced, which contains a list of
      \code{PointProcessModel}s. Each element in the list is a
      model of one coordinate given one or more of the other
      coordinates, and in total the list comprises a multivariate model
      of point processes. 
      \item Methods such as \code{summary}, \code{termPlot},
      \code{stepInformation} and \code{ppmFit} have been implemented for
      the \code{MultivariatePointProcess} class.
      \item A \code{getLocalIndependenceGraph} for class
      \code{MultivariatePointProcess} has been introduced that returns
      an \code{igraph} object representing the local (conditional)
      independence graphs for the multivariate point process object.
      \item The \code{pointProcessModel} function now interprets a
      vector of variables on the left hand side of a formula as the
      specification of a multivariate point process model in which case
      the function return an object of class
      \code{MultivariatePointProcess}.
    }
  }
   \subsection{Other Changes}{
     \itemize{
       \item The implementation of the linear filters has been modified
       slightly concerning the treatment of the boundary, and it has
       been reimplemented using the \pkg{Rcpp} interface. 
       \item A battery of unit tests based on \pkg{RUnit} has been added. 
    }
    }
}


