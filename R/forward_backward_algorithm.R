#' Calculating Forward and Backward Probabilities and Likelihood
#' 
#' The function calculates the logarithmized forward and backward probabilities and 
#' the logarithmized likelihood for a discrete time hidden Markov model, as defined 
#' in MacDonald & Zucchini (2009, Paragraph 3.1- Paragraph 3.3 and Paragraph 4.1).
#'
#' @param x a vector object containing the time-series of observations that are assumed 
#'          to be realizations of the (hidden Markov state dependent) observation process 
#'          of the model.
#' @param delta a vector object containing values for the marginal probability distribution 
#'    of the \code{m} states of the Markov chain at the time point \code{t=1}.
#' @param gamma a matrix (\code{nrow=ncol=m}) containing values for the transition matrix 
#'    of the hidden Markov chain.
#' @param distribution_class a single character string object with the abbreviated name of 
#'    the \code{m} observation distributions of the Markov dependent observation process. 
#'    The following distributions are supported: Poisson (\code{pois}); 
#'    generalized Poisson (\code{genpois}); normal (\code{norm}); geometric (\code{geom}).
#' @param distribution_theta a list object containing the parameter values for the \code{m} 
#'    observation distributions of the observation process that are dependent on the 
#'    hidden Markov state.
#' @param discr_logL a logical object. It is \code{TRUE} if the discrete log-likelihood 
#'    shall be calculated (for \code{distribution_class="norm"} ) instead of the general 
#'    log-likelihood.  See MacDonald & Zucchini (2009, Paragraph 1.2.3) 
#'    for further details.  Default is \code{FALSE}.
#' @param discr_logL_eps a single numerical value to approximately determine the discrete 
#'    log-likelihood for a hidden Markov model based on normal distributions 
#'    (for \code{"norm"}). The default value is \code{0.5}.  
#'    See MacDonald & Zucchini (2009, Paragraph 1.2.3) for further details. 
#'
#' @return \code{forward_backward_algorithm} returns a list containing the logarithmized 
#'    forward and backward probabilities and the logarithmized likelihood.
#'  \describe{
#'  \item{log_alpha}{a (T,m)-matrix (when T indicates the length/size of the observation 
#'        time-series and m the number of states of the HMM) containing the logarithmized 
#'        forward probabilities.}
#'  \item{log_beta}{a (T,m)-matrix (when T indicates the length/size of the observation 
#'       time-series and m the number of states of the HMM) containing the logarithmized 
#'       backward probabilities.}
#'  \item{logL}{a single numerical value representing the logarithmized likelihood.}
#'  \item{logL_calculation}{a single character string object which indicates how 
#'        \code{logL} has been calculated (see Zucchini (2009) 
#'        Paragraph 3.1-3.4, 4.1, A.2.2, A.2.3 for further details).}
#'  }
#'  
#' @references MacDonald, I. L., Zucchini, W. (2009) \emph{Hidden Markov Models for Time 
#'    Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
#'    
#' @author The basic algorithm for a Poisson-HMM is provided by MacDonald & Zucchini 
#'   (2009, Paragraph A.2.2). Extension and implementation by Vitali Witowski (2013).
#'   
#' @seealso \code{\link{HMM_based_method}}, \code{\link{HMM_training}},
#'    \code{\link{Baum_Welch_algorithm}}, \code{\link{direct_numerical_maximization}},
#'    \code{\link{initial_parameter_training}}   
#'    
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
#' # Calculating logarithmized forward/backward probabilities 
#' # and logarithmized likelihood
#' 
#' \donttest{
#' forward_and_backward_probabilities_and_logL <- 
#' forward_backward_algorithm (x = x, 
#'                         delta = delta, 
#'                         gamma = gamma, 
#'            distribution_class = distribution_class, 
#'            distribution_theta = distribution_theta)
#'            
#'  print(forward_and_backward_probabilities_and_logL)
#'  }
forward_backward_algorithm <- function(x, delta, gamma, distribution_class, 
                                       distribution_theta, discr_logL = FALSE, 
                                       discr_logL_eps = 0.5)
{
  size <- length(x)
  m <- length(delta)
  fb_svtpoue <- fb_small_value_to_prevent_overflow_and_underflow_errors <- 4.940656e-142   

  function_discr_log_L_p_norm <- function(x, mean, sd, discr_logL_eps)
  {
    foo <- pnorm((x + discr_logL_eps), mean = mean, sd = sd) - 
           pnorm((x - discr_logL_eps), mean = mean, sd = sd)	
  return(foo)
  }    
    
# Calcualtion of probabilities ---------------------------------------------------

  if (distribution_class == "pois")
  {	
    probabilities <- outer(X = x, Y = distribution_theta$lambda, FUN = dpois)
  }
  
  
  if (distribution_class == "geom")
  {  
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dgeom, 
                              prob = distribution_theta$prob))
  }
  

  if (distribution_class == "genpois")
  {
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dgenpois, 
                              lambda1 = distribution_theta$lambda1, 
                              lambda2 = distribution_theta$lambda2))
  }
  
  
  if (distribution_class == "norm" & discr_logL == FALSE)
  {
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN = dnorm, mean = distribution_theta$mean, sd = distribution_theta$sd))
  }
  
  
  if (distribution_class == "norm" & discr_logL == TRUE)
  {  
    probabilities <-  matrix(x, ncol = m, nrow = size)
    probabilities <-  t(apply(X = probabilities, MARGIN = 1, FUN =  function_discr_log_L_p_norm, mean = distribution_theta$mean, sd = distribution_theta$sd, discr_logL_eps = discr_logL_eps))
  }
  
      
  if (distribution_class == "bivariate_pois")
  {   
    size <- length(x[,1])
    
    probabilities <- matrix(0, ncol = m, nrow = size)
    for (i in 1:size)
    {
      for (j in 1:m)
      {
        probabilities[i,j] <- dpois(x[i,1], lambda=distribution_theta$lambda_1[j]) * dpois(x[i,2],lambda=distribution_theta$lambda_2[j])
      } 
    }	
  }

    probabilities <- ifelse(!is.na(probabilities), probabilities, fb_svtpoue)
    probabilities <- ifelse(!probabilities <= 0, probabilities, fb_svtpoue) 	
    probabilities <- ifelse(!probabilities > 1, probabilities, fb_svtpoue) 	
    probabilities <- ifelse(!probabilities == Inf, probabilities, fb_svtpoue)	
    probabilities <- ifelse(!probabilities == -Inf, probabilities, fb_svtpoue) 
    

  log_alpha <- matrix(fb_svtpoue, nrow = size, ncol = m)
  foo <- delta * probabilities[1,]
  sum_of_foo <- sum(foo) + fb_svtpoue
  scaled_logL <- log(sum_of_foo)
  foo <- foo / sum_of_foo
  log_alpha[1,] <- scaled_logL + log(foo)
  for (i in 2:size)
  {
    foo <- foo %*% gamma * probabilities[i,]
    sum_of_foo <- sum(foo) + fb_svtpoue
    scaled_logL <- scaled_logL + log(sum_of_foo)
    foo <- foo / sum_of_foo
    log_alpha[i,] <- scaled_logL + log(foo)
  }
  # Calculating log_L via alpha_T -------------------------------------
  logL_calculated_with_alpha_T <- scaled_logL
  
  

  log_beta <- matrix(fb_svtpoue, nrow = size, ncol = m)
  log_beta[size,] <- rep(0,m)
  foo <- rep(1 / m, m)
  scaled_logL <- log(m)
  for (i in (size-1):1)
  {
    foo <- gamma %*% (probabilities[i+1,] * foo)
    log_beta[i,] <- log(foo) + scaled_logL
    sum_of_foo <- sum(foo) + fb_svtpoue
    foo <- foo / sum_of_foo
    scaled_logL <- scaled_logL + log(sum_of_foo)
  }
  
  # Calculating log_L via beta_1 -----------------------------------------------------
  logL_calculated_with_beta_1 <- scaled_logL
  
  

  logL_calculated_with_alpha_t_and_beta_t <- 0
  middle_t <- round(size / 2)
  c <- max(log_alpha[middle_t,])
  for (i in 1:m)
  { 
  	logL_calculated_with_alpha_t_and_beta_t <- logL_calculated_with_alpha_t_and_beta_t + exp( log_alpha[middle_t,i] + log_beta[middle_t,i] - c)
  }
  logL_calculated_with_alpha_t_and_beta_t <- log(logL_calculated_with_alpha_t_and_beta_t) + c
  
  logL <- logL_calculated_with_alpha_t_and_beta_t
  logL_calculation <- "logL calculated with alpha_middle_t and beta_middle_t"
  if (logL == -Inf | logL == Inf | is.na(logL))
  {   
    logL_calculation <- "logL calculated with alpha_T"
    logL <- logL_calculated_with_alpha_T	
  }
  if (logL == -Inf | logL == Inf | is.na(logL))
  {
    logL <- logL_calculated_with_beta_1
    logL_calculation <- "logL calculated with beta_1"
  }

return(list(log_alpha = log_alpha, 
            log_beta = log_beta, 
            logL = logL, 
            logL_calculation = logL_calculation))			
}
