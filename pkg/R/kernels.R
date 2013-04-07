sobolevKernel <- function(s, r, t = 1, sub = NULL) {
  smin <- pmin(s, r)/t
  r <- pmax(s, r)/t
  s <- smin
  if (is.null(sub)) {
    1 + s * r / 48 +
      (1 - r) * s * (2 * r - r * r - s * s) / 6
  } else if (sub == 0) {
    1 + s * r / 48
  } else {
    (1 - r) * s * (2 * r - r * r - s * s) / 6
  }
}

gaussianKernel <- function(s, r, t = 1, c = 1) {
  s <- s/t
  r <- r/t
    exp( - c * (s - r)^2)
}
