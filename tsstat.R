# Copyright (C) Araik Tamazian, 2019


q.log <- function(x, q = 1) {
  # Computes the q-logarithm.
  #
  # Args:
  #   x: vector whose q-logarithm is to be calculated.
  #   qpar: q-parameter
  #
  # Returns:
  #   q-logarithm of x.
  if (q == 1) {
    y <- log(x)
  } else {
    y <- (x^(1 - q) - 1) / (1 - q)
  }
	return(y)
}

q.exp2 <- function(x, q = 1) {
  if (q == 1) {
    y <- exp(x)
  } else {
    if (1 + (1 - q) * x > 0) {
      y <- (1 + (1 - q) * x)^(1 / (1 - q))
    } else {
      ifelse(q < 1, 0, Inf)
    }
  }
  return(y)
}


q.exp <- function(x, q = 1) {
  # Computes the q-exponential.
  #
  # Args:
  #   x: vector whose q-exponential is to be calculated.
  #   q: parameter
  #
  # Returns:
  #   q-exponential of x.
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

# Q-exponential distribution --------------------------------------------------

dqexp <- function(x, q.shape = 1, rate = 1, log = FALSE) {
  # Computes the density function of the q-exponential distribution.
  #
  # Args:
  #   x: vector of quantiles.
  #   q.shape: shape parameter
  #   rate: rate parameter
  #   log: logical; if TRUE, probabilities p are given as log(p)
  #
  # Returns:
  #   Density function of q-exponential distribution.
  if (!is.logical(log.arg <- log) || length(log) != 1) {
    stop("bad input for argument 'log'")
  }   
  rm(log)
	if (q.shape >= 2 | rate <= 0) {
		pdf <- NaN
	} else {	
		pdf <- (2 - q.shape) * rate*q.exp(-x*rate, q.shape)
	}
  if (log.arg) {
    return(log(pdf))
  } else {
    return(pdf)
  } 
}

pqexp <- function(q, shape = 1, scale = 1, lower.tail = TRUE,
                  log.p = FALSE) {
  q.shape <- shape
  rate <- scale
  # Computes the distribution function of q-exponential distribution.
  #
  # Args:
  #   q: vector of quantiles.
  #   q.shape: shape parameter
  #   rate: rate parameter
  #   log.p: logical; if TRUE, probabilities p are given as log(p)
  #
  # Returns:
  #   Distribution function of q-exponential distribution.
  if (!is.logical(log.arg <- log.p) || length(log.p) != 1) {
    stop("bad input for argument 'log'")
  }   

  if (!is.logical(lt.arg <- lower.tail) || length(lower.tail) != 1) {
    stop("bad input for argument 'lower.tail'")
  }   

	if (q.shape >= 2 | rate <= 0) {
		cdf <- NaN
	} else {
		q.shape1 <- 1 / (2 - q.shape)
		cdf <- q.exp(- rate * q / q.shape1, q.shape1)
    if (lt.arg) {
      cdf <- 1 - cdf
    }		
	}
  if (log.arg) {
    return(log(cdf))
  } else {
    return(cdf)
  } 
}

qqexp <- function(p, q.shape = 1, rate = 1) {
  # Computes the quantile function of the q-exponential distribution.
  #
  # Args:
  #   p: vector of probabilities
  #   q.shape: shape parameter
  #   rate: rate parameter
  #   log.p: logical; if TRUE, probabilities p are given as log(p)
  #
  # Returns:
  #   Quantile function of q-exponential distribution.
	if (shape >= 2 | rate <= 0) {
		qfunc <- NaN
	}
	else {	
		q.shape1 <- 1 / (2 - q.shape)
		qf <- -q.shape1 * q.log(1 - p, q.shape1) / rate
	}
	return(qf)				
}

rqexp <- function(n, q.shape = 1, rate = 1) {
  # Generates random numbers with the q-exponential distribution.
  #
  # Args:
  #   n: number of observations.
  #   q.shape: shape parameter.
  #   rate: rate parameter.
  #
  # Returns:
  #   Random numbers with the q-exponential distribution.
	if (q.shape >= 2 | rate <= 0) {
		ans <- NaN
	}
	else {	
		alpha <- runif(n)
		q.shape1 <- 1 / (2 - q.shape)
		ans <- -q.shape1 * q.log (1 - alpha, q.shape1) / rate
	}
	return(ans)	
}

dqweibull <- function (x, qpar = 1, shape = 1, scale = 1, 
                       nc = 2 - qpar, log = FALSE)
{
  # Computes the density function of the q-weibull distribution.
  #
  # Args:
  #   x: vector of quantiles.
  #   qpar: q-parameter
  #   shape2: shape parameter
  #   scale: scale parameter
  #   log: logical; if TRUE, probabilities p are given as log(p)
  #
  # Returns:
  #   q-weibull density function of vector x.
  if (!is.logical(log.arg <- log) || length(log) != 1)
    stop("bad input for argument 'log'")
  rm(log)
  N <- length(x)
  pdf <- numeric(N)
  pdf <- (nc * shape * x^(shape - 1) / scale^shape) *
    q.exp(-(x / scale)^shape, qpar)
  if (log.arg) {
    return(log(pdf))
  } else {
    return(pdf)
  }   
}

pqweibull <- function (q, qpar = 1, shape = 1, scale = 1, lower.tail=FALSE) {
  qparm <- 1 / (2 - qpar)
  cdf <- q.exp(-(q / scale)^shape * (2 - qpar), qparm)
  if(lower.tail) {
    return(1-cdf)
  }
  else {
    return(cdf)
  }
}

rqweibull <- function (n, qpar = 1, shape = 1, scale = 1, 
                       nc = 2 - qpar) {
  # Generates q-weibull distributed random numbers.
  #
  # Args:
  #   n: number of observations.
  #   qpar: q-parameter.
  #   shape2: shape parameter.
  #   scale: scale parameter.
  #   nc: normalization constant. 
  #
  # Returns:
  #   n q-weibull distributed random numbers.
  u <- runif(n)
  qparm <- 1 / (2 - qpar)
  ncm <- nc / (2 - qpar)
  rand <- (-q.log(u / ncm, qparm) / (2 - qpar))^(1 / shape) * scale
  return(rand)
}

rqgauss <- function(n, location = 0, q.shape = 1, sigma = 1) {
  # Generates q-gaussian distributed random numbers.
  #
  # Args:
  #   n: number of observations.
  #   location: location parameter
  #   q.shape: shape parameter (q).
  #   scale: scale parameter. 
  #
  # Returns:
  #   n q-gaussian distributed random numbers.
  u <- runif(n)
  v <- runif(n)
  q.shape1 <- (1 + q.shape) / (3 - q.shape)
  x <- sqrt(-2 * q.log(u, q.shape1)) * cos(2 * pi * v)
  y <- sqrt(-2 * q.log(u, q.shape1)) * sin(2 * pi * v)
  xm <- location + x * sigma  #/ sqrt(scale * (3 - q.shape))
  ym <- location + y * sigma  #/ sqrt(scale * (3 - q.shape))
  rand <- c(xm, ym)

  return(rand)
}

rlqgauss <- function(n, location = 0, q.shape = 1, scale = 1) {
  # Generates q-gaussian distributed random numbers.
  #
  # Args:
  #   n: number of observations.
  #   location: location parameter
  #   q.shape: shape parameter (q).
  #   scale: scale parameter. 
  #
  # Returns:
  #   n q-gaussian distributed random numbers.
  u <- runif(n)
  v <- runif(n)
  q.shape1 <- (1 + q.shape) / (3 - q.shape)
  x <- sqrt(-2 * q.log(u, q.shape1)) * cos(2 * pi * v)
  y <- sqrt(-2 * q.log(u, q.shape1)) * sin(2 * pi * v)
  xm <- exp(location + x  / sqrt(scale * (3 - q.shape)))
  ym <- exp(location + y  / sqrt(scale * (3 - q.shape)))
  rand <- c(xm, ym)
  
  return(rand)
}

dqgauss <- function(x, location, shape, scale) {
#  nc[shape < 1] <- (2*sqrt(pi)*gamma(1/(1-shape)))/((3-shape)*sqrt(1-shape)*
#                                           gamma((3-shape)/(2-2*shape)))
#  nc[shape == 1] <- sqrt(pi)
  if(shape == 1) {
    nc <- sqrt(pi)
  } else {
    nc <- (sqrt(pi)*gamma((3-shape)/(2*shape-2)))/(sqrt(shape-1)*
                                                         gamma(1/(shape-1)))
  }
  
  sqrt(scale)*q.exp(-scale*(x-location)^2, shape)/nc
}

pqgauss <- function(q, location, shape, scale, lower.tail=TRUE) {
  if (lower.tail) {
    cumsum(dqgauss(q, location, shape, scale))
  } else {
    1-cumsum(dqgauss(q, location, shape, scale))
  }
}
