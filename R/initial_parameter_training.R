#' Algorithm to Find Plausible Starting Values for Parameter Estimation
#' 
#' The function computes plausible starting values for both the Baum-Welch algorithm and 
#' the algorithm for directly maximizing the log-Likelihood.  Plausible starting values 
#' can potentially diminish problems of (i) numerical instability and (ii) not finding 
#' the global optimum.
#'
#' @param x a vector object containing the time-series of observations that are assumed to be realizations of the (hidden Markov state dependent) observation process of the model.
#' @param m a (finite) number of states in the hidden Markov chain.
#' @param distribution_class a single character string object with the abbreviated name of 
#'    the $m$ observation distributions of the Markov dependent observation process.  
#'    The following distributions are supported:  Poisson (\code{pois}); 
#'    generalized Poisson (\code{genpois}, only available for 
#'    \code{training_method="numerical"}); normal (\code{norm})).
#' @param n a single numerical value specifying the number of samples to find the best 
#'    starting value for the training algorithm.  Default value is \code{100}.
#' @param discr_logL a logical object.  Default is \code{FALSE} for the general 
#'    log-likelihood, \code{TRUE} for the discrete log-likelihood 
#'    (for \code{distribution_class = "norm"}).
#' @param discr_logL_eps a single numerical value, used to approximate the discrete 
#'    log-likelihood for a hidden Markov model based on nomal distributions 
#'    (for \code{"norm"}).  The default value is \code{0.5}.
#'
#' @details
#' From our experience, parameter estimation for long time-series of observations 
#' (\code{T>1000}) or observation values \code{>1500} tend to be numerical unstable and 
#' does not necessarily find a global maximum.  Both problems can eventually be 
#' diminished with plausible starting values.  Basically, the idea behind 
#' \code{initial_parameter_training} is to sample randomly \code{n} sets of \code{m} 
#' observations from the time-series \code{x}, as means (\code{E}) of the state-dependent 
#' distributions. This \code{n} samplings of \code{E}, therefore induce \code{n} sets of 
#' parameters (\code{distribution_theta}) for the HMM without running a (slow) parameter 
#' estimation algorithm. Furthermore, \code{initial_parameter_training} calculates the 
#' log-Likelihood for all those \code{n} sets of parameters.  The set of parameters with 
#' the best Likelihood are outputted as plausible starting values.
#' 
#' (Additionally to the \code{n} sets of randomly chosen observations as means, the 
#' \code{m} quantiles of the observations are also checked as plausible means within 
#' this algorithm.)
#' 
#' @author Vitali Witowski (2013).
#' @seealso \code{\link{Baum_Welch_algorithm}}, \code{\link{direct_numerical_maximization}},
#'    \code{\link{HMM_training}}
#'    
#' @return
#' The function \code{ initial_parameter_training } returns a list containing the 
#' following components:
#' \describe{
#' \item{m}{input number of states in the hidden Markov chain.}
#' \item{k}{a single numerical value representing the number of parameters of the defined 
#'      distribution class of the observation process.}
#' \item{logL}{logarithmized likelihood of the model evaluated at the HMM with given 
#'       starting values (\code{delta, gamma, distribution theta}) induced by \code{E}.}
#' \item{E}{randomly choosen means of the observation time-series \code{x}, used for the 
#'      observation distributions, for which the induced parameters 
#' (\code{delta, gamma, distribution theta}) produce the largest Likelihood.}
#' \item{distribution_theta}{a list object containing the plausible starting values for 
#'      the parameters of the \code{m} observation distributions that are dependent on 
#'      the hidden Markov state.}
#' \item{delta}{a vector object containing plausible starting values for the marginal 
#'      probability distribution of the \code{m} states of the Markov chain at the time 
#'      point \code{t=1}.  At the moment:\cr \code{delta = rep(1/m, times=m)}.}
#' \item{gamma}{a matrix (\code{nrow=ncol=m}) containing the plausible starting values 
#'      for the transition matrix of the hidden Markov chain.  At the moment:\cr 
#'      \code{gamma = 0.8 * diag(m) + rep(0.2/m, times=m)}.} 
#' }
#' @keywords iteration

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
#' # Finding plausibel starting values for the parameter estimation 
#' # for a generealized-Pois-HMM with m=4 states
#' 
#' m <- 4 
#' \donttest{
#'  plausible_starting_values <- 
#'  initial_parameter_training(x = x, 
#'                             m = m, 
#'            distribution_class = "genpois", 
#'                             n = 100)
#'                             
#'  print(plausible_starting_values)
#'  }  
#' 
initial_parameter_training <-
function(x, m, distribution_class, n = 100, discr_logL =  FALSE, discr_logL_eps = 0.5)
{	
  par_var = 0.5
  gamma <- 0.8 * diag(m) + rep(0.2 / m, m)
  delta <- rep(1 / m, m)
  
  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k <- 1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois")
  {
    k <- 2
  }	

  if (distribution_class == "pois") 
  {	
  	
    # Building array-object to save 2+n different sets of parameters
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E","lambda")
    
    
    
    # Bilding the first set of parameters via 'quantile'-based picking for 
    # the mean values E, and saving them into the first place of the array
    E <- quantile(x = x, probs = seq(0, 1, 1 / (m)), na.rm = FALSE, names = FALSE, type = 7)
    E=E[1:m]
    E <- sort(E)
    distribution_theta <- list(lambda = E)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$lambda[j]
    }
    
    
    
    # Bilding the second set of parameters via 'some kind of 
    # uniformly-distribution'-picked values as E, 
    # and saving them into the second place of the array
    E <- seq((min(x)+1), max(x), length = m)
    E <- sort(E)
    distribution_theta <- list(lambda = E)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$lambda[j]
    }
    
    
    # Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, 
    # and saving them into the 3:n + 2 places of the array
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1, replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        distribution_theta <- list(lambda = E)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta,discr_logL =  discr_logL, discr_logL_eps = 0.5), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }       
        logL <- fb$logL
        erg[1,1,i] <- logL		
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$lambda[j]
        }
      }
    }	
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E = erg[,2,L_min], distribution_theta = list(lambda = erg[,3,L_min]), delta = delta, gamma = gamma))	
  }
  
    
  
  
  if (distribution_class == "geom") 
  {	

    # Building array-object to save 2+n different sets of parameters
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E", "prob")



    # Bilding the first set of parameters via 'quantile'-based picking for the 
    # mean values E, and saving them into the first place of the array
    E <- quantile(x = x, probs = seq(0, 1, 1 / (m)), na.rm = FALSE,names = F, type = 7)
    E=E[1:m]
    E <- sort(E)
    prob <- 1 / E
    distribution_theta <- list(prob = prob)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$prob[j]
    }



    # Bilding the second set of parameters via 'some kind of 
    # uniformly-distribution'-picked values as E, and saving them into the 
    # second place of the array 
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    prob <- 1 / E
    distribution_theta <- list(prob = prob)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$prob[j]
    }



    # Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, 
    # and saving them into the 3:n + 2 places of the array
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m-1))
        {	
          xt <- xt[!xt==E[h]]
          E <- append(E, sample(x = xt, 1,replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        prob <- 1 / E
        distribution_theta <- list(prob = prob)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta,discr_logL =  discr_logL,discr_logL_eps = 0.5), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }
        logL <- fb$logL
        erg[1,1,i] <- logL		
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        } 
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$prob[j]
        }
      }
    }
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta = list(prob=erg[,3,L_min]), delta = delta, gamma = gamma))	
  }
    
  
  if (distribution_class == "norm") 
  {	
    
    
    # Building array-object to save 2+n different sets of parameters
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E", "mean", "sd")
   
   
   
    # Bilding the first set of parameters via 'quantile'-based picking for the 
    # mean values E, and saving them into the first place of the array
    E <- quantile(x = x, probs = seq(0, 1, 1/(m)), na.rm = FALSE,names = F, type = 7)
    E <- E[1:m]
    E <- sort(E)
    sd_obs_pi <- rep(sd(x) / m, times = m)
    distribution_theta <- list(mean = E, sd = sd_obs_pi)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$mean[j]
    }
    for (j in 1:m)
    {
      erg[j,4,1] <- distribution_theta$sd[j]
    }
 
 
 
    # Bilding the second set of parameters via 'some kind of 
    # uniformly-distribution'-picked values as E, and saving them into the 
    # second place of the array
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    sd_obs_pi <- rep(sd(x) / m, times = m)
    distribution_theta <- list(mean = E, sd = sd_obs_pi)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$mean[j]
    }
    for (j in 1:m)
    {
      erg[j,4,2] <- distribution_theta$sd[j]
    }




    # Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, 
    # and saving them into the 3:n + 2 places of the array
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {   
        xt <- x
        E <- sample(x = x, 1, replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1,replace = FALSE, prob = NULL))
        }
        E <- sort(E)		
        sd_obs_pi <- rep(sd(x) / m, times = m)
        distribution_theta <- list(mean = E, sd = sd_obs_pi)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }
        logL <- fb$logL
        erg[1,1,i] <- logL
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$mean[j]
        }	
        for (j in 1:m)
        {
          erg[j,4,i] <- distribution_theta$sd[j]
        }	
      }
    }		
    L_min <- which.max(erg[1,1,])
    return(list(m = m,k = k,logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta = list(mean = erg[,3,L_min], sd = erg[,4,L_min]), delta = delta, gamma = gamma))
  }
  

  if (distribution_class == "genpois") 
  {	

    # Building array-object to save 2+n different sets of parameters
    erg <- array(NA, dim = c(m, 1 + 1 + k, (n + 2)))    
    colnames(erg)<-c("logL", "E","lambda1","lambda2")
 
 
    # Bilding the first set of parameters via 'quantile'-based picking for the 
    # mean values E, and saving them into the first place of the array
    E <- quantile(x = x, probs = seq(0, 1, 1/(m)), na.rm = FALSE,names = F, type = 7)
    E <- E[1:m]
    E <- sort(E)
    E[E == 0] <- 1
    lambda2 <- rep(par_var, times = m)
    lambda1 <- E * (1 - lambda2)
    distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,1] <- logL
    for (j in 1:m)
    {
      erg[j,2,1] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,1] <- distribution_theta$lambda1[j]
    }
    for (j in 1:m)
    {
      erg[j,4,1] <- distribution_theta$lambda2[j]
    }



    # Bilding the second set of parameters via 'some kind of 
    # uniformly-distribution'-picked values as E, and saving them into the 
    # second place of the array
    E <- seq((min(x) + 1), max(x), length = m)
    E <- sort(E)
    E[E == 0] <- 1
    lambda2 <- rep(par_var, times = m)
    lambda1 <- E * (1 - lambda2)
    distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
    fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
    if (inherits(fb,"try-error"))
    {  
      fb$logL <- -Inf
    }
    logL <- fb$logL
    erg[1,1,2] <- logL
    for (j in 1:m)
    {
      erg[j,2,2] <- E[j]
    }
    for (j in 1:m)
    {
      erg[j,3,2] <- distribution_theta$lambda1[j]
    }
    for (j in 1:m)
    {
      erg[j,4,2] <- distribution_theta$lambda2[j]
    }



    # Bilding the 3:n + 2 (loop) sets of parameters via 'randomly'-picked values as E, 
    # and saving them into the 3:n + 2 places of the array
    if ( n > 0)
    {
      for (i in 3:(n + 2))
      {  
        xt <- x
        E=sample(x = x, 1,replace = FALSE, prob = NULL)
        for (h in 1:(m - 1))
        {	
          xt <- xt[!xt == E[h]]
          E <- append(E, sample(x = xt, 1, replace = FALSE, prob = NULL))
        }
        E <- sort(E)
        E[E == 0] <- 1
        lambda2 <- rep(par_var, times = m)
        lambda1 <- E * (1 - lambda2)
        distribution_theta <- list(lambda1 = lambda1, lambda2 = lambda2)
        fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, distribution_class = distribution_class, distribution_theta = distribution_theta, discr_logL =  discr_logL, discr_logL_eps = discr_logL_eps), silent = FALSE)	
        if (inherits(fb,"try-error"))
        {  
          fb$logL <- -Inf
        }       
        logL <- fb$logL
        erg[1,1,i] <- logL
        for (j in 1:m)
        {
          erg[j,2,i] <- E[j]
        }
        for (j in 1:m)
        {
          erg[j,3,i] <- distribution_theta$lambda1[j]
        }	
        for (j in 1:m)
        {
          erg[j,4,i] <- distribution_theta$lambda2[j]
        }		
      }
    }
    L_min <- which.max(erg[1,1,])
    return(list(m = m, k = k, logL = erg[1,1,L_min], E=erg[,2,L_min], distribution_theta=list(lambda1=erg[,3,L_min], lambda2=erg[,4,L_min]), delta = delta, gamma = gamma))
  }
  
  
}
