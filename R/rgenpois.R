#' The Generalized Poisson Distribution
#' 
#' Density, distribution function and random generation function for the generalized 
#' Poisson distribution. 
#'
#' @param n number of observations
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
#' @seealso \code{\link{pgenpois}}, \code{\link{dgenpois}}; 
#'  \link{Distributions} for other standard distributions, 
#'  including \code{\link{dpois}} for the Poisson distribution.
#' 
#' @keywords distribution
#'
#' @return
#'  \code{\link{rgenpois}} generates random deviates of the generalized Poisson distribution.
#' @export
#'
#' @examples
#' dgenpois(x = seq(0,20), lambda1 = 10, lambda2 = 0.5) 
#' pgenpois(q = 5, lambda1 = 10, lambda2 = 0.5) 
#' hist(rgenpois(n = 1000, lambda1 = 10, lambda2 = 0.5) )

rgenpois <-function(n, lambda1, lambda2)
{
  random_genpois <- numeric(n)
  for (i in 1:n) 
  {
    temp_random_genpois <- 0
    random_number <- runif(1)
    kum <- dgenpois(0, lambda1 = lambda1, lambda2 = lambda2)
    while(random_number > kum) 
    {	
      temp_random_genpois <- temp_random_genpois + 1
      kum <- kum + dgenpois(temp_random_genpois, lambda1 = lambda1, lambda2 = lambda2)
    }
    random_genpois[i] <- temp_random_genpois
  }
  return(random_genpois)
}