###======================================================================
### ppstat - point process statistics
###======================================================================

## bSpline replaces bs for convenience and provides easier computations
## of B-spline bases that we are using. Does not issue a warning when
## evaluated outside of knots. This will in fact be the rule rather
## than the exception for point process usages.

bSpline <- function(x, knots, ..., sym = FALSE, trunc = NULL) {
  if(length(knots) <= 4) stop("Need at least 5 knots")
  if(sym) {
    design <- list()
    
    if(any(x>=0))
      design[[1]] <- bSpline(x = x[x>=0], knots = knots, ..., trunc = NULL, sym = FALSE)

    if(any(x < 0))
      design[[2]] <- bSpline(x = -x[x<0], knots = knots, ..., trunc = NULL, sym = FALSE)
    design <- do.call("rbind",design)
  } else {
    design <- splineDesign(knots,x,...,outer.ok=TRUE)
  }
  if(!is.null(trunc)) {
    if(length(trunc) == 1) {
      design[x <= trunc, ] <- 0
    } else if(length(trunc) == 2) {
      design[x <= trunc[1] | x > trunc[2], ] <- 0
    } else {
      stop("The truncation argument 'trunc' must be either a single numeric or a vector of length 2")
    }
  }
  design <- design[,apply(design,2,function(s) any(s != 0))]
  return(design)
}

## Truncated exponential

tExp <- function(x, tLevel=1e-4){
  y <- exp(x)
  y[y < tLevel] <- 0
  return(y)
}

## Constant function

const <- function(x, y, c=1){
  c*as.numeric(x<=y)
}

### TODO: This should be changed to an S4 method and specialized to the
### different point process classes.

print.summary.ppm <-  function (x, digits = max(3, getOption("digits") - 3), 
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


