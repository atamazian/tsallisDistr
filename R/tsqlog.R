#' @name tsqlog
#' @rdname tsqlog
#' @title q-logarithms and q-exponentials
#'
#' @description
#' \code{tsqlog} computes the q-logarithm function.
#' 
#' \code{tsqexp} computes the q-exponential function.
#' 
#' @param x: a numeric vector.
#' @param q: a \emph{q} parameter
#' 
#' @returns
#' A vector of the same length as \code{x} containing the transformed values. \code{tsqlog} and \code{tsqexp} gives the same results as \code{log} and \code{exp} if \code{q} is 1.
#' 
#' @references
#' Tsallis, Constantino (1994). "What are the numbers that experiments provide?". Química Nova. 17: 468.
#' @examples 
#' tsqlog(tsqexp(3));
#' tsqlog(tqexp(5, q = 2), q = 2);
NULL

#' @rdname tsqlog
#' @export
tsqlog <- function(x, q = 1) {
  if (q == 1) {
    y <- log(x)
  } else {
    y <- (x^(1 - q) - 1) / (1 - q)
  }
	return(y)
}

#' @rdname tsqlog
#' @export
tsqexp <- function(x, q = 1) {
  if (q == 1) {
    y <- sapply(x, exp)
  } else {
    xc <- (1 + (1 - q) * x > 0)
    y <- (1 + (1 - q) * x[xc])^(1 / (1 - q))
    y[!xc & q < 1] <- 0
    y[!xc & q > 1] <- Inf
  }
  return(y)
}
