#' Estimation by Directly Maximizing the log-Likelihood 
#' 
#' Estimates the parameters of a (stationary) discrete-time hidden Markov model by 
#' directly maximizing the log-likelihood of the model using the \link[stats]{nlm}-function. 
#' See MacDonald & Zucchini (2009, Paragraph 3) for further details.
#'
#' @param x  a vector object containing the time-series of observations that are assumed 
#'           to be realizations of the (hidden Markov state dependent) observation process 
#'           of the model.
#' @param m  interger ;a (finite) number of states in the hidden Markov chain
#' @param delta a vector object containing starting values for the marginal probability 
#'              distribution of the \code{m} states of the Markov chain at the time point 
#'              \code{t=1}. This implementation of the algorithm uses the stationary 
#'              distribution as delta. 
#' @param gamma a matrix (\code{nrow=ncol=m}) containing starting values for the transition 
#'              matrix of the hidden Markov chain
#' @param distribution_class a single character string object with the abbreviated name of 
#'           the \code{m} observation distributions of the Markov dependent observation 
#'           process. The following distributions are supported by this algorithm: 
#'           Poisson (\code{pois}); generalized Poisson (\code{genpois}); 
#'           normal (\code{norm}, discrete log-Likelihood not applicable by this algorithm)
#' @param distribution_theta a list object containing starting values for the parameters 
#'           of the \code{m} observation distributions that are dependent on the hidden 
#'           Markov state
#' @param DNM_limit_accuracy a single numerical value representing the convergence 
#'          criterion of the direct numerical maximization algorithm using the 
#'          \link[stats]{nlm}-function. Default value is \code{0.001}.
#' @param DNM_max_iter a single numerical value representing the maximum number of 
#'          iterations of the direct numerical maximization using the 
#'          \link[stats]{nlm}-function. Default value is \code{50}.
#' @param DNM_print a single numerical value to determine the level of printing of the 
#'          \code{nlm}-function.  See \code{nlm}-function for further informations. 
#'          The value \code{0} suppresses, that no printing will be outputted. 
#'          Default value is \code{2} for full printing.
#'
#' @return \code{direct_numerical_maximization } returns a list containing the estimated 
#'   parameters of the hidden Markov model and other components.
#'   \describe{
#'   \item{x}{input time-series of observations.}
#'   \item{m}{input number of hidden states in the Markov chain.}
#'   \item{logL}{a numerical value representing the logarithmized likelihood calculated 
#'         by the \code{\link{forward_backward_algorithm}}.}
#'   \item{AIC}{a numerical value representing Akaike's information criterion for the 
#'         hidden Markov model with estimated parameters.}
#'   \item{BIC}{a numerical value representing the Bayesian information criterion for the 
#'         hidden Markov model with estimated parameters.}
#'   \item{delta}{a vector object containing the estimates for the marginal probability 
#'        distribution of the \code{m} states of the Markov chain at time-point 
#'        point \code{t=1}.}
#'   \item{gamma}{a matrix containing the estimates for the transition matrix of the 
#'        hidden Markov chain.}
#'   \item{distribution_theta}{a list object containing estimates for the parameters of 
#'        the \code{m} observation distributions that are dependent on the hidden Markov state.}
#'   \item{distribution_class}{input distribution class.}
#'   }
#'   
#' @references MacDonald, I. L., Zucchini, W. (2009) \emph{Hidden Markov Models for Time 
#'             Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
#' 
#' @author The basic algorithm of a Poisson-HMM is provided by MacDonald & Zucchini 
#'   (2009, Paragraph A.1). Extension and implementation by Vitali Witowski (2013).
#'   
#' @seealso \code{\link{HMM_based_method}}, \code{\link{HMM_training}}, 
#'   \code{\link{Baum_Welch_algorithm}}, 
#'   \code{\link{forward_backward_algorithm}}, 
#'   \code{\link{initial_parameter_training}} 
#' @keywords ts 
#' @export
#'
#' @examples
#' x <- c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,
#'   14,16,13,11,11,10,12,19,23,25,24,23,20,21,22,22,18,7,
#'   5,3,4,3,2,3,4,5,4,2,1,3,4,5,4,5,3,5,6,4,3,6,4,8,9,12,
#'   9,14,17,15,25,23,25,35,29,36,34,36,29,41,42,39,40,43,
#'   37,36,20,20,21,22,23,26,27,28,25,28,24,21,25,21,20,21,
#'   11,18,19,20,21,13,19,18,20,7,18,8,15,17,16,13,10,4,9,
#'   7,8,10,9,11,9,11,10,12,12,5,13,4,6,6,13,8,9,10,13,13,
#'   11,10,5,3,3,4,9,6,8,3,5,3,2,2,1,3,5,11,2,3,5,6,9,8,5,
#'   2,5,3,4,6,4,8,15,12,16,20,18,23,18,19,24,23,24,21,26,
#'   36,38,37,39,45,42,41,37,38,38,35,37,35,31,32,30,20,39,
#'   40,33,32,35,34,36,34,32,33,27,28,25,22,17,18,16,10,9,
#'   5,12,7,8,8,9,19,21,24,20,23,19,17,18,17,22,11,12,3,9,
#'   10,4,5,13,3,5,6,3,5,4,2,5,1,2,4,4,3,2,1) 
#'   
#'   # Assumptions (number of states, probability vector, 
#'   # transition matrix, and distribution parameters)
#'   
#'   m <- 4
#'   delta <- c(0.25,0.25,0.25,0.25)
#'   gamma <- 0.7 * diag(m) + rep(0.3 / m)
#'   distribution_class <- "pois"
#'   distribution_theta <- list(lambda = c(4,9,17,25))
#'   
#'   # Estimation of a HMM using the method of 
#'   # direct numerical maximization
#'   \donttest{
#'     trained_HMM_with_m_hidden_states <- 
#'       direct_numerical_maximization(x = x, 
#'                                     m = m, 
#'                                 delta = delta, 
#'                                 gamma = gamma, 
#'                    distribution_class = distribution_class,
#'                          DNM_max_iter = 100,
#'                    distribution_theta = distribution_theta)
#'                    
#'     print(trained_HMM_with_m_hidden_states)
#'  }   

direct_numerical_maximization <-
function(x, m,  delta, gamma, distribution_class, distribution_theta, DNM_limit_accuracy = 0.001, DNM_max_iter = 50, DNM_print = 2)
{
  # This algorithm hasn't been extended for the use with the discrete likelihood, yet 
  discr_logL <- FALSE     
  discr_logL_eps <- 0.5   
  
  
  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k <- 1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois" | distribution_class == "bivariate_pois")
  {
    k <- 2
  }	
  
  
  DNM_log_n2w <- function(np)
  {
    wp <- log(np)
    return(wp)
  }
  
  
  DNM_exp_w2n <- function(wp)
  {
    np <- exp(wp)
    return(np)
  }
  
  
  DNM_logit_n2w <- function(np)
  {
    wp=log(np / (1 - np))
    return(wp)
  }
  
  
  DNM_invlogit_w2n <- function(wp)
  {
    np=exp(wp) / (1 + exp(wp))
    return(np)
  } 
  
  
  DNM_n2w <- function(m, gamma, distribution_class, distribution_theta)
  {	
    
    if (distribution_class == "pois")
    {	
      
      tlambda <- distribution_theta$lambda
      tlambda <- DNM_log_n2w(tlambda)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tlambda,tgamma)
    }
    
    if (distribution_class == "norm")
    {	
      tmean <- distribution_theta$mean
      tsd <- distribution_theta$sd
      tsd <- DNM_log_n2w(tsd)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tmean, tsd, tgamma)
    }
    
    if (distribution_class == "genpois")
    {	
      tlambda1 <- distribution_theta$lambda1
      tlambda1 <- DNM_log_n2w(tlambda1)
      tlambda2 <- distribution_theta$lambda2
      tlambda2 <- DNM_logit_n2w(tlambda2)
      tgamma <- NULL
      if (m > 1)
      {
        foo <- log(gamma / diag(gamma))
        tgamma <- as.vector(foo[!diag(m)])
      }
      parameter_vector <- c(tlambda1, tlambda2, tgamma)
    }
    return(parameter_vector)	
  }
   
   
  DNM_w2n <- function(parameter_vector, gamma, distribution_class, distribution_theta, m)
  {
    if (distribution_class == "pois")
    {
      distribution_theta$lambda <- DNM_exp_w2n(parameter_vector[1:m])
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((1 * m) + 1):((1 * m)+(m * m - m))]
        gamma <- gamma/apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m)-gamma+1),rep(1,m))
    }
    
    if (distribution_class == "norm")
    {
      distribution_theta$sd <- DNM_exp_w2n(parameter_vector[(m + 1):(m + m)])
      distribution_theta$mean <- parameter_vector[1:m]
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((2 * m) + 1):(2 * m+(m * m - m))]
        gamma <- gamma/apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m)-gamma+1),rep(1,m))
    }
    
    if (distribution_class == "genpois")
    {
      distribution_theta$lambda2 <- DNM_invlogit_w2n(parameter_vector[(m + 1):(m + m)])
      distribution_theta$lambda1 <- DNM_exp_w2n(parameter_vector[1:m])
      gamma <- diag(m)
      epar <- exp(parameter_vector)
      if (m > 1)
      {
        gamma[!gamma] <- epar[((2 * m)+1):(2 * m+(m * m - m))]
        gamma <- gamma / apply(gamma,1,sum)
      }
      delta <- solve(t(diag(m) - gamma + 1), rep(1,m))
    }
    return(list(delta=delta,
                gamma=gamma,
                distribution_theta = distribution_theta))
  }

  
  function_neglogL <- function(x, p, distribution_class, distribution_theta, m)
  {
    underflow_error_1=FALSE
    underflow_error_2=FALSE
    pos_logL_error=FALSE
    worse_logL_error=FALSE
    np <- DNM_w2n(parameter_vector = p, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    fb <-  try( forward_backward_algorithm(x = x, gamma = np$gamma, delta = np$delta, distribution_class = distribution_class, distribution_theta = np$distribution_theta,  discr_logL = discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)
    if (inherits(fb, "try-error"))
    {  
      underflow_error_1 <- TRUE
      fb$logL <- -Inf
    }
    neglogL <- (-1) * fb$logL
    
  return(neglogL)	
  }
  
  if (distribution_class == "pois")
  {	
    if (any(distribution_theta$lambda == 0))
    { 
      distribution_theta$lambda <- ifelse(!distribution_theta$lambda == 0, distribution_theta$lambda, 1.1)
    }
    
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class,  distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    delta <- par$delta
    
    gamma <- par$gamma
    
    distribution_theta <- par$distribution_theta
    
    estimated_mean_values <- distribution_theta$lambda
  }
  
  if (distribution_class == "norm")
  {	
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m = m, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    delta=par$delta
    
    gamma=par$gamma
    
    distribution_theta=par$distribution_theta
    
    estimated_mean_values <- distribution_theta$mean
  }
  
  if (distribution_class == "genpois")
  {	
    if (any(distribution_theta$lambda1 == 0))
    { 
      distribution_theta$lambda1 <- ifelse(!distribution_theta$lambda1 == 0, distribution_theta$lambda1, 0.1)
    }
    
    parameter_vector0 <- DNM_n2w(m = m, gamma = gamma, distribution_class = distribution_class, distribution_theta = distribution_theta)
    
    z <- nlm(f = function_neglogL, p = parameter_vector0, x = x, m = m, distribution_class = distribution_class,  distribution_theta = distribution_theta, print.level = DNM_print, gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
    
    logL <- - z$minimum
    
    par <- DNM_w2n(parameter_vector = z$estimate, m=m, distribution_class = distribution_class,distribution_theta = distribution_theta)
    
    delta <- par$delta
    gamma <- par$gamma
    distribution_theta <- par$distribution_theta
    estimated_mean_values <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2)
  }
  
  AIC <- AIC_HMM(logL = logL, m = m, k = k) 
  BIC <- BIC_HMM(size = length(x), logL = logL, m = m, k = k) 
  
  return(list(x = x, 
              m = m, 
              logL = logL, 
              AIC = AIC, 
              BIC = BIC, 
              delta = delta, 
              gamma = gamma, 
              distribution_class = distribution_class, 
              distribution_theta = distribution_theta, 
              estimated_mean_values=estimated_mean_values))
}
