#' @name tsqweibulldistr
#' @rdname tsqweibulldistr
#' @title The q-Weibull distribution
#'
#' @description
#' Density, distribution function, quantile functions and random generation for the \emph{q}-Weibull distribution with shapes \code{shape1} and \code{shape2} and scale \code{scale}. 
#' 
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param shape1,shape2,scale shape and scale parameters, the latter defaulting to 1.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \emph{P[X ≤ x]}, otherwise, \emph{P[X > x]}.
#' 
#' @details 
#' If \code{shape1} is equal to 1, the \emph{q}-Weibull distribution becomes the Weibull distribution.
#' 
#' @returns 
#' \code{dtsqweibull} gives the density, \code{ptsweibull} gives the distribution function, \code{qtsexp} gives the quantile function, and \code{rtsweibull} generates random deviates.
#' 
#' The length of the result is determined by \code{n} for \code{rtsweibull}, and is the maximum of the lengths of the numerical arguments for the other functions.
#' 
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' 
#' @note 
#' The cumulative hazard \emph{H(t) = -log(1 - F(t))} is \code{-ptsweibull(t, q, r, lower = FALSE, log = TRUE)}.
#' 
#' @references 
#' Picoli, S. Jr.; Mendes, R. S.; Malacarne, L. C. (2003). "q-exponential, Weibull, and q-Weibull distributions: an empirical analysis". Physica A: Statistical Mechanics and Its Applications. 324 (3): 678–688. arXiv:cond-mat/0301552. Bibcode:2003PhyA..324..678P. doi:10.1016/S0378-4371(03)00071-2. S2CID 119361445.
#' 
#' @source 
#' \code{[dpq]tsweibull} are calculated directly from the definitions, \code{rtsweibull} uses inversion.
#' 
#' @seealso 
#' \code{\link[=tsqlog]{tsqlog()}} for the \emph{q}-logarithm and \emph{q}-exponential
#' 
#' @examples 
#' x <- c(0, rlnorm(50))
#' all.equal(dtsqweibull(x, shape1 = 1, shape2 = 2, scale = pi), dweibull(x, shape = 2, scale = pi))
#' all.equal(ptsqweibull(x, shape1 = 1, shape2 = 2, scale = pi), pweibull(x, shape = 2, scale = pi))
#' all.equal(qtsqweibull(x/11, shape1 = 1, shape2 = 2, scale = pi), qweibull(x/11, shape = 2, scale = pi))
NULL

#' @rdname tsqweibulldistr
#' @export
dtsqweibull <- function (x, shape1, shape2, scale = 1, log = FALSE)
{
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  N <- length(x)
  pdf <- numeric(N)
  pdf <- (( 2 - shape1) * shape2 * x^(shape2 - 1) / scale^shape2) *
    tsqexp(-(x / scale)^shape2, shape1)
  if (log.arg) {
    return(log(pdf))
  } else {
    return(pdf)
  }   
}

#' @rdname tsqweibulldistr
#' @export
ptsqweibull <- function (q, shape1, shape2, scale = 1, lower.tail=TRUE, log.p = FALSE) {
  rate1m <- scale/(2 - shape1)^(1/shape2)
  shape1m <- 1 / (2 - shape1)
  cdf <- tsqexp(-(q / rate1m)^shape2, shape1m)
  if(log.p) {
    qf = log(qf)
  }
  if(lower.tail) {
    return(1-cdf)
  }
  else {
    return(cdf)
  }
}

#' @rdname tsqweibulldistr
#' @export
qtsqweibull <- function (p, shape1, shape2, scale = 1, lower.tail=TRUE, log.p = FALSE) {
  rate1m <- scale/(2 - shape1)^(1/shape2)
  shape1m <- 1 / (2 - shape1)
  if(lower.tail) {
    p <- 1 - p
  }
  qf <- rate1m*(-tsqlog(p, shape1m))^(1/shape2)
  if(log.p) {
    qf = log(qf)
  }
  return(qf)
}

#' @rdname tsqweibulldistr
#' @export
rtsqweibull <- function (n, shape1, shape2, scale = 1) {
  u <- runif(n)
  shape1m <- 1 / (2 - shape1)
  rand <- (-tsqlog(u, shape1m) / (2 - shape1))^(1 / shape2) * scale
  return(rand)
}
