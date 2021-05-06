#' @name tsqgaussdistr
#' @rdname tsqgaussdistr
#' @title The q-Gaussian distribution
#'
#' @description
#' Density, distribution function, quantile functions and random generation for the \emph{q}-Gaussian distribution with shape \code{shape} and scale \code{scale}. 
#' 
#' @param x vector of quantiles
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param shape shape parameter
#' @param location location parameter
#' @param scale  scale parameter
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' 
#' @details 
#' If \code{shape1} is equal to 1, the \emph{q}-Gaussian distribution becomes the normal distribution.
#' 
#' @returns 
#' \code{dtsqgauss} gives the density, and \code{rtsqgauss} generates random deviates.
#' 
#' The length of the result is determined by \code{n} for \code{rtsqgauss}, and is the maximum of the lengths of the numerical arguments for the other functions.
#' 
#' The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used.
#' 
#' @references 
#' Umarov, Sabir; Tsallis, Constantino; Steinberg, Stanly (2008). "On a q-Central Limit Theorem Consistent with Nonextensive Statistical Mechanics" (PDF). Milan J. Math. Birkhauser Verlag. 76: 307–328. doi:10.1007/s00032-008-0087-y. S2CID 55967725.
#'  W. Thistleton, J.A. Marsh, K. Nelson and C. Tsallis, Generalized Box–Muller method for generating q-Gaussian random deviates, IEEE Transactions on Information Theory 53, 4805 (2007)
#' 
#' @source 
#' \code{[dtsqgauss} is calculated directly from the definitions, \code{rtsqgauss} uses inversion.
#' 
#' @seealso 
#' \code{\link[=tsqlog]{tsqlog()}} for the \emph{q}-logarithm and \emph{q}-exponential
#' 
#' @examples 
#' x <- c(0, rlnorm(50))
#' all.equal(dtsqgauss(x, shape = 1, scale = 2), dnorm(x, mean = 0, sd = 2))
NULL

#' @rdname dtsqgauss
#' @export
dtsqgauss <- function(x, shape, location = 0, scale = 1, log = FALSE) {
  if(shape == 1) {
    nc <- sqrt(pi)
  } else {
    nc <- (sqrt(pi)*gamma((3-shape)/(2*shape-2)))/(sqrt(shape-1)*
                                                         gamma(1/(shape-1)))
  }
  pdf <- sqrt(scale)*q.exp(-scale*(x-location)^2, shape)/nc
  
  if (log.arg) {
    return(log(pdf))
  } else {
    return(pdf)
  }   
}

rtsqgauss <- function(n, shape, location = 0, scale = 1) {
  u <- runif(n)
  v <- runif(n)
  shape1 <- (1 + shape) / (3 - shape)
  x <- sqrt(-2 * q.log(u, shape1)) * cos(2 * pi * v)
  y <- sqrt(-2 * q.log(u, shape1)) * sin(2 * pi * v)
  xm <- location + x * scale  #/ sqrt(scale * (3 - shape))
  ym <- location + y * scale  #/ sqrt(scale * (3 - shape))
  rand <- c(xm, ym)

  return(rand)
}
