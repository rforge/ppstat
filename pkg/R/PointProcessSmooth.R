pointProcessSmooth <- function(
                               formula,
                               data,
                               family,
                               support,
                               lambda = 1,
                               allKnots = FALSE,
                               N = 200,
                               Delta,
                               basisPoints,
                               coefficients,
                               fixedCoefficients = list(),
                               fit = TRUE,
                               varMethod = 'Fisher',
                               basisEnv,
                               ...) {
  
  call <- match.call()
  argList <- as.list(call)[-1]
  argList$fit <- FALSE

  if(missing(basisPoints))
    {
      if(!(missing(support)))
        {
          if(length(support) == 1)
            support <- c(0,max(support[1],0))
        } else {
          stop("Must specify either 'support' or 'basisPoints'.")
        }
    } else {
      support = range(basisPoints)
    }
  ### Construction of the basis expansions
  terms <- terms(formula, "s")
  specials <- attr(terms, "specials")$s
  specialVar <- lapply(as.list(attr(terms, "variables"))[1+specials], all.vars)
  ## if(is.list(specialVar))
  ##   stop("Wrong specification of some smoother term.")
  ## if(length(specials) != length(specialVar))
  ##   stop("Some smoother term is applied to more than one variable.")

  if(attr(terms, "response") == 0)
    stop("No response variable specified.") 

  response <- all.vars(as.list(attr(terms, "variables"))[[1+attr(terms, "response")]])
  termLabels <- attr(terms, "term.labels")
  specialTerms <- which(apply(attr(terms, "factor")[specials, , drop = FALSE] > 0, 2, any))

  if(allKnots) {
    strategy <- "all"
  } else {
    strategy <- "log"
  }
  
  knots <- list()
  fList <- list()
  ## TODO: Clean this up! Ask on R-devel if there is
  ## a problem with environments for functions ....
  ## Is this simply lazy evaluation ... ?
  termFunction <- function(knots) {
    fknots <- knots
    function(x) bSpline(x, knots = fknots)
  }
  for(i in seq_along(specialVar)) {
    x <- getPointPosition(data)[getMarkType(data) %in% response]
    y <- getPointPosition(data)[getMarkType(data) %in% specialVar[[i]]]
    knots[[specialTerms[i]]] <- computeKnots(x, y, support, specialVar[[i]], strategy)
    term <- paste("fList[[", specialTerms[i], "]](", specialVar[[i]], ")", sep = "")
    fList[[specialTerms[i]]] <- termFunction(knots = knots[[specialTerms[i]]])
    termLabels[specialTerms[i]] <- term
  }

  formula <- reformulate(termLabels, response = response)

  argList$formula <- formula
  
  model <- do.call("pointProcessModel", argList)
  ## TODO: Modify colnames for the model matrix. 
  nrCoef <- dim(getModelMatrix(model))[2]
  Omega <- matrix(0, ncol = nrCoef, nrow = nrCoef)

  for(i in seq_along(specialTerms)) {
    penCoef <- which(attr(getModelMatrix(model), "assign") == specialTerms[i])
    d <- length(penCoef)
    s1 <- s2 <- s3 <- s4 <- numeric(d)
    s <-  .Fortran("sgram", as.double(s1), as.double(s2), as.double(s3), as.double(s4), as.double(knots[[specialTerms[i]]]), as.integer(d))
    pen <- matrix(0, d, d)
    diag(pen) <- s[[1]]
    pen[seq(2,d*d,d+1)] <- pen[seq(d+1,d*d,d+1)] <- s[[2]][1:(d-1)]
    pen[seq(3,d*(d-1),d+1)] <- pen[seq(2*d+1,d*d,d+1)] <- s[[3]][1:(d-2)]
    pen[seq(4,d*(d-2),d+1)] <- pen[seq(3*d+1,d*d,d+1)] <- s[[4]][1:(d-3)]
    Omega[penCoef, penCoef] <- pen
    }

  model@Omega <- lambda*Omega
  model@penalization <- TRUE
  
  if(fit) 
    model <- ppmFit(model, ...)
  
  model@call <- call
  model <- as(model, "PointProcessSmooth")
  return(model)
}

computeKnots <- function(x, y, support, variables, strategy = "log", method = "s", ...) {
  ## TODO: Implement this in C!?
  differences <- outer(x, y, '-')
  differences <- differences[differences > support[1] & differences < support[2]]
  differences <- unique(differences)
  ## TODO: Implement different strategies for "thinning". This one is taken from
  ## smooting.spline directly ...
  sknotl <- function(x, nk = "log")
    {
      ## if (!all.knots)
      ## return reasonable sized knot sequence for INcreasing x[]:
      n.kn <- function(n) {
        ## Number of inner knots
        if(n < 50L) n
        else trunc({
          a1 <- log( 50, 2)
          a2 <- log(100, 2)
          a3 <- log(140, 2)
          a4 <- log(200, 2)
          if	(n < 200L) 2^(a1+(a2-a1)*(n-50)/150)
          else if (n < 800L) 2^(a2+(a3-a2)*(n-200)/600)
          else if (n < 3200L)2^(a3+(a4-a3)*(n-800)/2400)
          else  200 + (n-3200)^0.2
        })
      }
      n <- length(x)
      if(isTRUE(nk == "log")){
        nk <- n.kn(n)
      } else if(isTRUE(nk == "all")) {
        nk <- n
      } 
      else if(!is.numeric(nk)) stop("'nknots' must be numeric <= n")
      else if(nk > n)
        stop("Cannot use more inner knots than unique 'x' values.")
      c(rep(x[1L], 3L), x[seq.int(1, n, length.out= nk)], rep(x[n], 3L))
    }
  
  knots <- sknotl(c(support[1], sort(differences), support[2]), nk = strategy) 
  return(knots)
}

## TODO: New summary function for an object of class 'PointProcessSmooth'.
## TODO: New update function. Model matrix needs to be recomputed if we change response

setMethod("summary","PointProcessSmooth",
          function(object,...) {
            callNextMethod()
          }
          )
