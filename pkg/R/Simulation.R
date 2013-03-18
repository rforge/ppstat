### Implementation of the multivariate Ogata thinning algorithm.

Ogata <- function(n = 1, lambda, h, hMax, A = Inf, seed = NULL, ...) {
  if (!is.null(seed))
    set.seed(seed)
  T <- matrix(0, n, 2)
  colnames(T) <- c("time", "mark")
  t <- 0
  i <- 0
  lamb <- lambda(0, h = h, ...)  ## Vector
  K <- sum(lamb)
  M <- length(lamb)
  Tlocal <- vector("list", M)
  while(i < n) {
    S <- rexp(1, K)
    U <- runif(1)
    t <- t + S
    tmA <- t - A
    for (m in seq_len(M)) {
      j <- 1
      while (j <= length(Tlocal[[m]]) && Tlocal[[m]][j] < tmA)
        j <- j + 1
      if (j > 1)
        Tlocal[[m]] <- Tlocal[[m]][-(1:(j-1))]
    }
    lamb <- lambda(t = t, T = Tlocal, h = h, ...)  ## Vector
    if(U < sum(lamb)/K) {
      i <- i + 1
      ## ALternatively, the next sampling step could be done
      ## using the U already simulated. 
      m <- sample.int(M, 1, prob = lamb)
      T[i, ] <- c(t, m)
      Tlocal[[m]] <- c(Tlocal[[m]], t)
    }
    K <- sum(lambda(t = t, T = Tlocal, h = hMax, predict = FALSE, ...))
  }
  T
} 

hawkesRate <- function(t, T = list(), h, predict = TRUE, Delta = 1, 
                       beta0 = rep(1, M), phi = function(x) pmax(x, 0), 
                       warn = TRUE, ...) {
  if (!is.list(T))
    stop("Argument 'T' must be a list.")
  t <- t[1]
  if (!is.list(h) && length(h) != length(T)) 
    stop("Argument 'h' is not a list of the same length af 'T'.")
  M <- length(h)
  iMax <- length(h[[1]][[1]])
  lamb <- beta0
  for (k in seq_len(M)) {
    if (length(T) == 0 || length(T[[k]]) == 0) {
      i <- FALSE
    } else {
      if (t < T[[k]][length(T[[k]])] || (t == T[[k]][length(T[[k]])] && predict)) 
        stop("Time argument 't' smaller than largest 'T'. Possible explosion.")
      i <- floor((t - T[[k]])/Delta + 1.5)
      if (i[length(i)] > iMax) {
        if (warn)
          warning("Point history truncated.", call. = FALSE)
        i <- i[i <= iMax]
      }
    }
    for (m in seq_len(M)) {
      if (!is.null(h[[m]][[k]]))
        lamb[m] <- lamb[m] + sum(h[[m]][[k]][i])
    }
  }
  phi(lamb)
}

setMethod("simulate", "MultivariatePointProcess",
          function(object, nsim = 1, seed = NULL, ...) {
            Delta <- object@models[[1]]@Delta
            A <- max(object@models[[1]]@support)
            linearFilters <- getLinearFilter(object)
            varNames <-  unique(unlist(lapply(linearFilters, 
                                              function(x) names(x)[-1]))
                                ) 
            M <- length(linearFilters) 
            h <- vector("list", M)
            hMax <- vector("list", M)
            beta0 <- numeric(M)
            for(m in seq_len(M)) {
              h[[m]] <- hMax[[m]] <- vector("list", length(varNames))
              names(h[[m]]) <- names(hMax[[m]]) <- varNames
              N <- nrow(linearFilters[[m]])
              for(k in seq_along(varNames)) {
                if(varNames[k] %in% names(linearFilters[[m]])) {
                  h[[m]][[varNames[k]]] <- linearFilters[[m]][, varNames[k]]
                } else {
                  h[[m]][[varNames[k]]] <- rep(0, N)
                }
                hMax[[m]][[varNames[k]]] <- rev(cummax(rev(h[[m]][[varNames[k]]])))
                
              }
              beta0[m] <- coefficients(object@models[[m]])["(Intercept)"]
            }
            Ogata(nsim, 
                  lambda = hawkesRate, 
                  h = h,
                  hMax = hMax,
                  A = A,
                  Delta = Delta,
                  beta0 = beta0,
                  seed = seed
            )
          }
)
