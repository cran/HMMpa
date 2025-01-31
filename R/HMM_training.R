#' Training of Hidden Markov Models
#' 
#' Function to estimate the model specific parameters 
#' (\code{delta, gamma, distribution_theta}) for a hidden Markov model, given a 
#' time-series and a user-defined distribution class. Can also be used for model 
#' selection (selecting the optimal number of states \code{m}). 
#' See Details for more information. 
#'
#' @param x a vector object of length \code{T} containing observations of a time-series 
#'    \code{x}, which are assumed to be realizations of the (hidden Markov state dependent) 
#'    observation process of the HMM.
#' @param distribution_class a single character string object with the abbreviated name of 
#'    the $m$ observation distributions of the Markov dependent observation process.  
#'    The following distributions are supported:  Poisson (\code{pois}); 
#'    generalized Poisson (\code{genpois}, only available for 
#'    \code{training_method="numerical"}); normal (\code{norm})).
#' @param min_m minimum number of hidden states in the hidden Markov chain. 
#'    Default value is \code{2}.
#' @param max_m maximum number of hidden states in the hidden Markov chain. 
#'    Default value is \code{6}.
#' @param n a single numerical value specifying the number of samples to find the best 
#'    starting values for the training algorithm.  Default value is \code{n=100}.
#' @param training_method a logical value indicating whether the Baum-Welch algorithm 
#'    (\code{"EM"}) or the method of direct numerical maximization (\code{"numerical"}) 
#'    should be applied for estimating the model specific parameters. 
#'    See \code{\link{Baum_Welch_algorithm}} and 
#'    \code{\link{direct_numerical_maximization}} for further details. 
#' @param discr_logL a logical object.  Default is \code{FALSE} for the general 
#'    log-likelihood, \code{TRUE} for the discrete log-likelihood 
#'    (for \code{distribution_class = "norm"}).
#' @param discr_logL_eps a single numerical value, used to approximate the discrete 
#'    log-likelihood for a hidden Markov model based on nomal distributions 
#'    (for \code{"norm"}).  The default value is \code{0.5}.
#' @param Mstep_numerical a logical object indicating whether the Maximization Step of 
#'    the Baum-Welch algorithm should be performed by numerical maximization. 
#'    Default is \code{FALSE}.
#' @param dynamical_selection a logical value indicating whether the method of dynamical 
#'    initial parameter selection should be applied (see Details).  Default is \code{TRUE}.
#' @param BW_max_iter  a single numerical value representing the maximum number of 
#'    iterations in the Baum-Welch algorithm.  Default value is \code{50}.
#' @param BW_limit_accuracy  a single numerical value representing the convergence 
#'    criterion of the Baum-Welch algorithm. Default value is is \code{0.001}.
#' @param BW_print  a logical object indicating whether the log-likelihood at each 
#'    iteration-step shall be printed. Default is \code{TRUE}.
#' @param DNM_max_iter  a single numerical value representing the maximum number of 
#'    iterations of the numerical maximization using the nlm-function 
#'    (used to perform the Maximization Step of the Baum-Welch-algorithm).  Default value is \code{50}.
#' @param DNM_limit_accuracy a single numerical value representing the convergence 
#'    criterion of the numerical maximization algorithm using the \link[stats]{nlm} 
#'    function (used to perform the Maximization Step of the Baum-Welch- algorithm).  
#'    Default value is \code{0.001}.
#' @param DNM_print a single numerical value to determine the level of printing of the 
#'    \code{nlm}-function.  See \code{nlm}-function for further informations. 
#'    The value \code{0} suppresses, that no printing will be outputted. Default 
#'    value is \code{2} for full printing.
#'
#' @details
#' More precisely, the function works as follows:\cr
#' \bold{Step 1:}
#' In a first step, the algorithm estimates the model specific parameters for different 
#' values of \code{m} (indeed for \code{min_m,...,max_m}) using either the function 
#' \code{\link{Baum_Welch_algorithm}} or \cr
#' \code{\link{direct_numerical_maximization}}. Therefore, the function first searches for 
#' plausible starting values by using the function \code{\link{initial_parameter_training}}. 
#' \bold{Step 2:}
#' In a second step, this function evaluates the AIC and BIC values for each HMM 
#' (built in Step 1) using the functions \code{\link{AIC_HMM}} and \code{\link{BIC_HMM}}. 
#' Then, based on that values, this function decides for the most plausible number of 
#' states \code{m} (respectively for the most appropriate HMM for the given time-series 
#' of observations). In case when AIC and BIC claim for a different \code{m}, the 
#' algorithm decides for the smaller value for \code{m} (with the background to have a 
#' more simplistic model). 
#' If the user is intereseted in having a HMM with a fixed number for \code{m}, 
#' \code{min_m} and \code{max_m} have to be chosen equally (for instance 
#' \code{min_m=4 = max_m} for a HMM with \code{m=4} hidden states).
#' To speed up the parameter estimation for each \eqn{m > m_min}, the user can choose the 
#' method of dynamical initial parameter selection.
#' If the method of dynamical intial parameter selection \bold{is not applied}, the 
#' function 
#' \code{\link{initial_parameter_training}} will be called to find plausible starting 
#' values for each state 
#' \eqn{ m \in \{min_m, \ldots, max_m\}}{ min_m \le m \le max_m}. \cr
#' If the method of dynamical intial parameter selection \bold{is applied}, 
#' then starting parameter values using the function \code{\link{initial_parameter_training}} 
#' will be found only for the first HMM (respectively the HMM with \code{m_min} states). 
#' The further starting parameter values for the next HMM (with \code{m+1} states and so on) 
#' are retained from the trained parameter values of the last HMM (with \code{m} 
#' states and so on). 
#' 
#' @return
#' \code{HMM_training} returns a list containing the following components:
#' \describe{
#' \item{trained_HMM_with_selected_m}{a list object containing the key data of the optimal 
#'      trained HMM (HMM with selected \code{m}) -- summarized output of the 
#'      \code{\link{Baum_Welch_algorithm}} or \cr \code{\link{direct_numerical_maximization}} 
#'      algorithm, respectively.}
#' \item{list_of_all_initial_parameters}{a list object containing the plausible starting 
#'      values for all HMMs (one for each state \code{m}).}
#' \item{list_of_all_trained_HMMs}{a list object containing all trained m-state-HMMs. 
#'      See \code{\link{Baum_Welch_algorithm}} or \code{\link{direct_numerical_maximization}} 
#'      for \code{training_method="EM"} or 
#' \code{training_method="numerical"}, respectively.}
#' \item{list_of_all_logLs_for_each_HMM_with_m_states}{a list object containing all 
#'      logarithmized Likelihoods of each trained HMM.}
#' \item{list_of_all_AICs_for_each_HMM_with_m_states}{a list object containing the AIC 
#'      values of all trained HMMs.}
#' \item{list_of_all_BICs_for_each_HMM_with_m_states}{a list object containing the BIC 
#'      values of all trained HMMs.}
#' \item{model_selection_over_AIC}{is logical.  \code{TRUE}, if model selection was based 
#'      on AIC and \code{FALSE}, if model selection was based on BIC.}
#' }
#' 
#' @references MacDonald, I. L.,  Zucchini, W. (2009) 
#'    \emph{Hidden Markov Models for Time Series: An Introduction Using R}, 
#'    Boca Raton: Chapman & Hall.
#'    
#' @author Vitali Witowski (2013)
#' @seealso \code{\link{initial_parameter_training}}, \code{\link{Baum_Welch_algorithm}},
#'    \code{\link{direct_numerical_maximization}}, \code{\link{AIC_HMM}}, 
#'    \code{\link{BIC_HMM}}
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
#' # Train a poisson hidden Markov model using the Baum-Welch 
#' # algorithm for different number of states m=2,...,6
#' \donttest{
#'  trained_HMMs <- 
#'   HMM_training(x = x, 
#'            min_m = 2, 
#'            max_m = 6, 
#'      distribution_class = "pois", 
#'      training_method = "EM")
#'      
#' # Various output values for the HMM 
#' names(trained_HMMs)
#' 
#' # Print details of the most plausible HMM for the given 
#' # time-series of observations
#' print(trained_HMMs$trained_HMM_with_selected_m)
#' 
#' # Print details of all trained HMMs (by this function) 
#' # for the given time-series of observations
#' print(trained_HMMs$list_of_all_trained_HMMs)
#' 
#' # Print the BIC-values of all trained HMMs for the given 
#' # time-series of observations  
#' print(trained_HMMs$list_of_all_BICs_for_each_HMM_with_m_states)
#' 
#' # Print the logL-values of all trained HMMs for the 
#' # given time-series of observations  
#' print(trained_HMMs$list_of_all_logLs_for_each_HMM_with_m_states)
#' }   
#' 
HMM_training <-
function(x, distribution_class, min_m = 2, max_m = 6, 
         n = 100, training_method = "EM", discr_logL = FALSE, discr_logL_eps = 0.5, 
         Mstep_numerical = FALSE, dynamical_selection = TRUE, 
         BW_max_iter = 50, BW_limit_accuracy = 0.001, BW_print = TRUE,
         DNM_max_iter = 50, DNM_limit_accuracy = 0.001, DNM_print = 2)
{
  par_var <- 0.5 
  how_many_HMMs <- length(min_m:max_m) + (min_m - 1)
  
  # Define the arrays to save the estimates for each m-state-HMM -----
  list_of_all_initial_parameters <- pairlist(NA)
  list_of_all_trained_HMMs <- pairlist(NA)
  list_of_all_logLs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  list_of_all_AICs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  list_of_all_BICs_for_each_HMM_with_m_states <- rep(NA, times = min_m - 1)
  
  for (i in min_m:max_m) 
  {	
  	# print each step to see for which m the HMM is trained 
  	print(paste('HMM with m =', toString(i)))
    
    # intial parameter training/setting -----
    if (dynamical_selection == FALSE) 
    {
      temp_initial_parameters <- initial_parameter_training(n = n, x = x, m = i, distribution_class = distribution_class,  discr_logL= discr_logL, discr_logL_eps = discr_logL_eps)
    }			
    if (dynamical_selection == TRUE) 
    {
      if (i == min_m) 
      {
        temp_initial_parameters <- initial_parameter_training(n = n, x = x, m = i, distribution_class = distribution_class, discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)
      }		
      if (i > min_m) 
      {	
        temp_initial_parameters$m <- i				
        temp_initial_parameters$delta <- c(temp_trained_HMM_regarding_distribution_theta$delta , 1 / i)
        temp_initial_parameters$delta <- temp_initial_parameters$delta / sum(temp_initial_parameters$delta) 				
        temp_initial_parameters$gamma <- matrix(0.2 / i, nrow = i, ncol = i)
        temp_initial_parameters$gamma[i,i] <- temp_initial_parameters$gamma[i,i] + 0.8
        for (j in 1:(i - 1)) 
        {
          temp_initial_parameters$gamma[j,1:(i - 1)] <- temp_trained_HMM_regarding_distribution_theta$gamma[j,]
        }
        for (j in 1:i) 
        {
          temp_initial_parameters$gamma[j,] <- temp_initial_parameters$gamma[j,] / sum(temp_initial_parameters$gamma[j,])
        }				
        temp_distance <- x
        for (j in 1:length(x)) 
        {
          temp_distance[j] <- min(abs(x[j] - temp_trained_HMM_regarding_distribution_theta$estimated_mean_values))
        }
        temp_next_estimated_mean_value <- x[which.max(temp_distance)]
        temp_initial_parameters$E <- c(temp_trained_HMM_regarding_distribution_theta$estimated_mean_values, temp_next_estimated_mean_value)				
        if (distribution_class == "pois") 
        {
          temp_initial_parameters$distribution_theta <- list(lambda=temp_initial_parameters$E)
        }
        if (distribution_class == "genpois") 
        {
          temp_initial_parameters$distribution_theta <- list(lambda1 = temp_initial_parameters$E * (1 - par_var) , lambda2 = rep(par_var, times = i))
        }
        if (distribution_class == "norm") 
        {
          temp_initial_parameters$distribution_theta <- list(mean = temp_initial_parameters$E, sd = rep(sd(x) / i, times = i))
        }
        if (distribution_class == "geom") 
        {
          temp_initial_parameters$distribution_theta <- list(prob = 1 / temp_initial_parameters$E)
        }				
        temp_initial_parameters$logL <- forward_backward_algorithm(x = x, delta = temp_initial_parameters$delta, gamma = temp_initial_parameters$gamma, distribution_class = distribution_class, distribution_theta = temp_initial_parameters$distribution_theta,  discr_logL = discr_logL, discr_logL_eps = discr_logL_eps)$logL
      }
    }		
    list_of_all_initial_parameters[[i]] <- temp_initial_parameters
    
    # choosing either the Baum-Welch algorithm or the method of 
    # directly maximizing the likelihood 
    if (training_method == "EM") 
    {
      temp_trained_HMM_regarding_distribution_theta <- Baum_Welch_algorithm(x = x, m = temp_initial_parameters$m, gamma = temp_initial_parameters$gamma, delta = temp_initial_parameters$delta, distribution_class = distribution_class, distribution_theta = temp_initial_parameters$distribution_theta, discr_logL = discr_logL, Mstep_numerical = Mstep_numerical, discr_logL_eps = discr_logL_eps,  BW_max_iter = BW_max_iter, BW_limit_accuracy = BW_limit_accuracy, DNM_max_iter = DNM_max_iter, DNM_limit_accuracy = DNM_limit_accuracy,  DNM_print = DNM_print, BW_print = BW_print)
    }
    if (training_method == "numerical") 
    {
      temp_trained_HMM_regarding_distribution_theta = direct_numerical_maximization(x = x, m = temp_initial_parameters$m, distribution_class = distribution_class, gamma = temp_initial_parameters$gamma, distribution_theta = temp_initial_parameters$distribution_theta, DNM_max_iter = DNM_max_iter, DNM_limit_accuracy = DNM_limit_accuracy, DNM_print = DNM_print)
    }		
    list_of_all_trained_HMMs[[i]] <- temp_trained_HMM_regarding_distribution_theta
    
    list_of_all_logLs_for_each_HMM_with_m_states <- c(list_of_all_logLs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$logL)
    list_of_all_AICs_for_each_HMM_with_m_states <- c(list_of_all_AICs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$AIC)
    list_of_all_BICs_for_each_HMM_with_m_states <- c(list_of_all_BICs_for_each_HMM_with_m_states, temp_trained_HMM_regarding_distribution_theta$BIC)
  }	
  
  # Selection of the HMM with the most plausible m

  AIC_selection_of_m <- which.min(list_of_all_AICs_for_each_HMM_with_m_states)
  BIC_selection_of_m <- which.min(list_of_all_BICs_for_each_HMM_with_m_states)
  
  model_selection_over_AIC <- TRUE
  selected_m <- AIC_selection_of_m
  if(AIC_selection_of_m > BIC_selection_of_m)
  {
    selected_m <- BIC_selection_of_m
    model_selection_over_AIC <- FALSE
  }
  
  trained_HMM_with_selected_m <- list_of_all_trained_HMMs[[selected_m]]
  
  
  trained_HMM_with_selected_m <- list(x = trained_HMM_with_selected_m$x, 
                                      m = trained_HMM_with_selected_m$m, 
                                      logL = trained_HMM_with_selected_m$logL, 
                                      AIC = trained_HMM_with_selected_m$AIC, 
                                      BIC = trained_HMM_with_selected_m$BIC, 
                                      delta = trained_HMM_with_selected_m$delta, 
                                      gamma = trained_HMM_with_selected_m$gamma, 
                                      distribution_class = trained_HMM_with_selected_m$distribution_class, 
                                      distribution_theta= trained_HMM_with_selected_m$distribution_theta, 
                                      estimated_mean_values= trained_HMM_with_selected_m$estimated_mean_values)
                                      
 
print(trained_HMM_with_selected_m)
  
return(list(trained_HMM_with_selected_m = trained_HMM_with_selected_m, 
            list_of_all_initial_parameters = list_of_all_initial_parameters, 
            list_of_all_trained_HMMs = list_of_all_trained_HMMs, 
            list_of_all_logLs_for_each_HMM_with_m_states = list_of_all_logLs_for_each_HMM_with_m_states, 
            list_of_all_AICs_for_each_HMM_with_m_states = list_of_all_AICs_for_each_HMM_with_m_states, 
            list_of_all_BICs_for_each_HMM_with_m_states = list_of_all_BICs_for_each_HMM_with_m_states, 
            model_selection_over_AIC=model_selection_over_AIC))
}
