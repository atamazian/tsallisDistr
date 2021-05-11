#' @name tsqexpdistr
#' @rdname tsqexpdistr
#' @title The q-exponential distribution
#'
#' @description
#' Density, distribution function, quantile functions and random generation for the \emph{q}-exponential distribution with shape \code{shape} and rate \code{rate}. 
#' 
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param shape,rate shape and rate parameters, the latter defaulting to 1.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \emph{P[X ≤ x]}, otherwise, \emph{P[X > x]}.
#' 
#' @details 
#' If \code{shape} and \code{rate} are not specified, they assume the default value of 1.
#' 
#' If \code{shape} is equal to 1, the \emph{q}-exponential distribution becomes the exponential.
#' 
#' @returns 
#' \code{dtsqexp} gives the density, \code{ptsqexp} gives the distribution function, \code{qtsexp} gives the quantile function, and \code{rtsqexp} generates random deviates.
#' 
#' The length of the result is determined by \code{n} for \code{rtsqexp}, and is the maximum of the lengths of the numerical arguments for the other functions.
#' 
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' 
#' @note 
#' The cumulative hazard \emph{H(t) = -log(1 - F(t))} is \code{-ptsqexp(t, q, r, lower = FALSE, log = TRUE)}.
#' 
#' @references 
#' Tsallis, C. Nonadditive entropy and nonextensive statistical mechanics-an overview after 20 years. Braz. J. Phys. 2009, 39, 337–356.
#' Picoli, S. Jr.; Mendes, R. S.; Malacarne, L. C. (2003). "q-exponential, Weibull, and q-Weibull distributions: an empirical analysis". Physica A: Statistical Mechanics and Its Applications. 324 (3): 678–688. arXiv:cond-mat/0301552. Bibcode:2003PhyA..324..678P. doi:10.1016/S0378-4371(03)00071-2. S2CID 119361445.
#' 
#' @source 
#' \code{[dpq]tsqexp} are calculated directly from the definitions, \code{rtsqexp} uses inversion.
#' 
#' @seealso 
#' \code{\link[=tsqlog]{tsqlog()}} for the \emph{q}-logarithm and \emph{q}-exponential
#' 
#' @examples 
#' x <- c(0, rlnorm(50))
#' all.equal(dtsqexp(x, shape = 1), dexp(x))
#' all.equal(ptsqexp(x, shape = 1, rate = 1/pi), pexp(x, rate = 1/pi))
#' all.equal(qtsqexp(x/11, shape = 1, rate = 1/pi), qexp(x/11, rate = 1/pi))
NULL

#' @rdname tsqexpdistr
#' @export
dtsqexp <- function(x, shape, rate = 1, log = FALSE) {
  if (shape >= 2) {
    stop("bad parameter value: 'shape' must be < 2")
  }
  if (rate <= 0) {
    stop("bad parameter value: 'rate' must be > 0")
  }
  if (!is.logical(log.arg <- log) || length(log) != 1) {
    stop("bad input for argument 'log'")
  }   
  rm(log)
  pdf <- (2 - shape) * rate*tsqexp(-x*rate, shape)
  if (log.arg) {
    return(log(pdf))
  } else {
    return(pdf)
  } 
}

#' @rdname tsqexpdistr
#' @export
ptsqexp <- function(q, shape, rate = 1, lower.tail = TRUE,
                  log.p = FALSE) {
  if (shape >= 2) {
    stop("bad parameter value: 'shape' must be < 2")
  }
  if (rate <= 0) {
    stop("bad parameter value: 'rate' must be > 0")
  }
  if (!is.logical(log.arg <- log.p) || length(log.p) != 1) {
    stop("bad input for argument 'log'")
  }   
  
  if (!is.logical(lt.arg <- lower.tail) || length(lower.tail) != 1) {
    stop("bad input for argument 'lower.tail'")
  }   
  
  shape1 <- 1 / (2 - shape)
  cdf <- tsqexp(- rate * q / shape1, shape1)
  if (lt.arg) {
    cdf <- 1 - cdf
  }		
  if (log.arg) {
    return(log(cdf))
  } else {
    return(cdf)
  } 
}

#' @rdname tsqexpdistr
#' @export
qtsqexp <- function(p, shape, rate = 1) {
  if (shape >= 2) {
    stop("bad parameter value: 'shape' must be < 2")
  }
  if (rate <= 0) {
    stop("bad parameter value: 'rate' must be > 0")
  }
  shape1 <- 1 / (2 - shape)
  qf <- -shape1 * tsqlog(1 - p, shape1) / rate
  return(qf)				
}

#' @rdname tsqexpdistr
#' @export
rtsqexp <- function(n, shape, rate = 1) {
  if (shape >= 2) {
    stop("bad parameter value: 'shape' must be < 2")
  }
  if (rate <= 0) {
    stop("bad parameter value: 'rate' must be > 0")
  }
  alpha <- runif(n)
  shape1 <- 1 / (2 - shape)
  ans <- -shape1 * tsqlog (1 - alpha, shape1) / rate
  return(ans)	
}
