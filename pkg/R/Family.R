### TODO: Implement a family "Event". This family should be responsible for
### models where we assume that events form different "individuals" all have 
### the same "zero-base" and where we allow for an at risk process.

### TODO: Implement a family "Gibbs" where we allow for anticipation in the 
### stochastic integrals.

setMethod("initialize", "Family",
          function(.Object, name = "Hawkes", link = "identity", c = 0, phi = NULL, Dphi = NULL, D2phi = NULL){
            .Object@name <- name
            if(is.null(phi)) {
              tmpLink <- makeLink(link,c=c)
              .Object@phi <- tmpLink$linkinv
              .Object@Dphi <- tmpLink$Dlinkinv
              .Object@D2phi <- tmpLink$D2linkinv
              .Object@link <- link
            } else {
              .Object@phi <- phi
              .Object@Dphi <- Dphi
              .Object@D2phi <- D2phi
              .Object@link <- NA
            }
            return(.Object)
          }
        )

Hawkes <- function(link = "root", ...) {
  return(new("Family", "Hawkes", link = link,...))
}

Gibbs <- function(link = "root", ...) {
  return(new("Family", "Gibbs", link = link, ...))
}

setMethod("family", "Family",
          function(object,...) {
            return(object@name)
          }
          )

makeLink <- function (link,c=0) 
{
    switch(link, logit = {
        linkfun <- function(mu) .Call("logit_link", mu, PACKAGE = "stats")
        linkinv <- function(eta) .Call("logit_linkinv", eta, 
            PACKAGE = "stats")
        Dlinkinv <- function(eta) .Call("logit_mu_eta", eta, PACKAGE = "stats")
        D2linkinv <- function(eta) .Call("logit_mu_eta", eta, PACKAGE = "stats")*(1-.Call("logit_linkinv", eta, PACKAGE = "stats"))
        valideta <- function(eta) TRUE
    }, cloglog = {
        linkfun <- function(mu) log(-log(1 - mu))
        linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)), 
            1 - .Machine$double.eps), .Machine$double.eps)
        Dlinkinv <- function(eta) {
            eta <- pmin(eta, 700)
            pmax(exp(eta-exp(eta)), .Machine$double.eps)
        }
        D2linkinv <- function(eta) {
            eta <- pmin(eta, 700)
            pmax((1-exp(eta))*exp(eta-exp(eta)), .Machine$double.eps)
        }
        valideta <- function(eta) TRUE
    }, identity = {
        linkfun <- function(mu) mu
        linkinv <- function(eta) eta
        Dlinkinv <- function(eta) rep.int(1, length(eta))
        D2linkinv <- function(eta) rep.int(0,length(eta))
        valideta <- function(eta) all(eta > 0)
    }, log = {
        linkfun <- function(mu) log(mu)
        linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
        Dlinkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
        D2linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, log2 = {
        linkfun <- function(mu) log(mu, 2)
        linkinv <- function(eta) pmax(2^eta, .Machine$double.eps)
        Dlinkinv <- function(eta) pmax(2^eta, .Machine$double.eps)
        D2linkinv <- function(eta) pmax(2^eta, .Machine$double.eps)
        valideta <- function(eta) TRUE
    }, root = {
      if(c==0) {
        linkfun <- function(mu) mu
        linkinv <- function(eta) pmax(eta*(eta > 0), .Machine$double.eps)
        Dlinkinv <- function(eta) rep(1,length(eta))*(eta > 0)
        D2linkinv <- function(eta) rep(0,length(eta))*(eta > 0) 
        valideta <- function(eta) TRUE
      } else {
        linkfun <- function(mu) mu^{1/(c+1)}
        linkinv <- function(eta) pmax(eta^(c+1)*(eta > 0), .Machine$double.eps^(c+1))
        Dlinkinv <- function(eta) (c+1) * pmax(eta^c*(eta > 0), .Machine$double.eps^c)
        D2linkinv <- function(eta) (c+1)*c*pmax(eta^{c-1}*(eta > 0), .Machine$double.eps^(c-1))
        valideta <- function(eta) TRUE
      }        
    }, logaffine = {
        linkfun <- function(mu) log(mu)(mu <= exp(c)) + (exp(-c)*mu + c - 1)*(mu > exp(c))
        linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)*(eta <=c) + exp(c)*(eta -c + 1)*(eta > c)
        Dlinkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)*(eta <=c) + exp(c)*(eta > c)
        D2linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)*(eta <=c)
        valideta <- function(eta) TRUE
    }, stop(sQuote(link), " link not recognised"))
    structure(list(linkfun = linkfun, linkinv = linkinv, Dlinkinv = Dlinkinv, D2linkinv = D2linkinv,
        valideta = valideta, name = link), class = "link-glppm")
}
