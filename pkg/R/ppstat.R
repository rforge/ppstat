###======================================================================
### ppstat - point process statistics
###======================================================================

## bSpline replaces bs for convenience and provides easier computations
## of B-spline bases that we are using. Does not issue a warning when
## evaluated outside of knots. This will in fact be the rule rather
## than the exception.

bSpline <- function(x,knots,...,trunc=0,whichColumns=NULL){
  if(length(knots) <= 4) stop("Need at least 5 knots")
  design <- splineDesign(knots,x,...,outer.ok=TRUE)
  if(!is.null(whichColumns)) design <- design[,whichColumns]
  if(!is.null(trunc)) design[x <= trunc,] <- 0
  return(design)
}

## Truncated exponential

tExp <- function(x,tLevel=1e-4){
  y <- exp(x)
  y[y < tLevel] <- 0
  return(y)
}

## Constant function

const <- function(x,y,c=1){
  c*as.numeric(x<=y)
}

### TODO: This should be changed to an S4 method and specialized to the
### different point process classes.

print.summary.glppm <-  function (x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    if(length(x$fixedCoefficients$which) == 1) {
      cat("\nNote that", rownames(x$coefficients)[x$fixedCoefficients$which], "was fixed\n\n")
    } else if(length(x$fixedCoefficients$which) > 1) {
      cat("\nNote that", paste(rownames(x$coefficients)[x$fixedCoefficients$which],collapse=", "), "were fixed\n\n")
    }
    if(x$penalization) cat("\nWarning: Parameters estimated with a quadratic penalty term.\nCorrections of the degrees of freedom is currently not implemented\nand AIC is not corrected either. The estimated standard errors may be\ninaccurate and the univariate test statistics can be misleading.\n\n")
    cat("\nMinus-log-likelihood: ", format(x$mll, digits = max(4, digits + 1)), " on ",x$df, " degrees of freedom",
        "\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), 
        "\n\n", "Number of function evaluations: ", x$iter[1],
        "\n","Number of gradient evaluations: ", x$iter[2],
        "\n", sep = "")
    if(x$convergence!=0) cat("\nWarning: Algorithm did not converge.\n         'optim' convergence status:",x$convergence,"\n")  
    cat("\n")
    invisible(x)
}


