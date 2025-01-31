#' The Generalized Poisson Distribution
#' 
#' Density function for the generalized Poisson distribution. 
#' 
#' @param x a vector of (non-negative integer) quantiles
#' @param lambda1 a single numeric value for parameter \code{lambda1} with \eqn{lambda1 > 0}
#' @param lambda2 a single numeric value for parameter \code{lambda2} with \eqn{0 \le lamdba2 < 1}.  
#'                When \code{lambda2=0}, the generalized Poisson distribution 
#'                reduces to the Poisson distribution
#' 
#' @details
#' The generalized Poisson distribution has the density
#' \deqn{ p(x) = \lambda_1 (\lambda_1 + \lambda_2 \cdot x)^{x-1} 
#'   \frac{ \exp(-\lambda_1-\lambda_2 \cdot x) )}{x!}}{%
#'         p(x) = lambda1 (lambda1 + lambda2 x)^(x-1)  exp(-lambda1-lambda2 x) )/x!}
#'   for \eqn{x = 0,1,2,\ldots},b
#'   with \eqn{\mbox{E}(X)=
#'   \frac{\lambda_1}{1-\lambda_2}}{E(x)=lambda1/(1-lambda2)} and variance 
#'   \eqn{\mbox{var}(X)=\frac{\lambda_1}{(1-\lambda_2)^3}}{var(x)=lambda1/(1-lambda2)^3}.
#'   
#' @references Joe, H., Zhu, R. (2005). Generalized poisson distribution: the property of 
#' mixture of poisson and comparison with negative binomial distribution. 
#' Biometrical Journal \bold{47}(2):219--229. 
#' 
#' @author Based on Joe and Zhu (2005). Implementation by  Vitali Witowski (2013).
#' 
#' @seealso \code{\link{pgenpois}}, \code{\link{rgenpois}}; 
#'  \link{Distributions} for other standard distributions, 
#'  including \code{\link{dpois}} for the Poisson distribution.
#' 
#' @keywords distribution
#'
#' @return
#'  \code{\link{dgenpois}} gives the density of the generalized Poisson distribution.
#' @export
#'
#' @examples
#' dgenpois(x = seq(0,20), lambda1 = 10, lambda2 = 0.5) 
#' pgenpois(q = 5, lambda1 = 10, lambda2 = 0.5) 
#' hist(rgenpois(n = 1000, lambda1 = 10, lambda2 = 0.5) )

dgenpois <- function(x, lambda1, lambda2)
{ 
	
  if (length(x) < max(length(lambda1), length(lambda2)))
  {
  	x <- c(rep(x, times = max(length(lambda1), length(lambda2))))
  }
  if (length(lambda1) < max(length(x), length(lambda2)))
  {
  	lambda1 <- c(rep(lambda1, times = max(length(x), length(lambda2))))
  }
  if (length(lambda2) < max(length(x), length(lambda1)))
  {
  	lambda2 <- c(rep(lambda2, times = max(length(x), length(lambda1))))
  }
  a <- NULL
  for (j in 1:max(length(x), length(lambda1), length(lambda2)))
  { 
  	if(x[j] < 2)
  { 
  	b = (lambda1[j] * (lambda1[j] + x[j] * lambda2[j])^(x[j] - 1) * (exp(-(
    lambda1[j] + x[j] * lambda2[j])))) / (factorial(x[j]))
  } else  { 
  	  f1 <- (lambda1[j] + x[j] * lambda2[j])
      g1 <- exp(-lambda2[j])
      e <- 1
      for (i in 2:x[j])
      { 
      	e <- e * ((f1 * g1) / i)
      }
      d1 <- lambda1[j] * exp(-lambda1[j]) * g1
      b <- e * d1
      if (is.na(b))
      {
      	b <- 4.940656e-324
      }
      if(b == Inf)
      {
      	b <- 4.940656e-324
      }
      if(b == -Inf)
      {
      	b <- 4.940656e-324
      }
    }
    a <- c(a,b)
  }
return(a)
}
