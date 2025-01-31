#' Estimation Using the Baum-Welch Algorithm
#' 
#' Estimates the parameters of a (non-stationary) discrete-time hidden Markov model. 
#' The Baum-Welch algorithm is a version of the EM (Estimation/Maximization) algorithm.  
#' See MacDonald & Zucchini (2009, Paragraph 4.2) for further details. 
#'
#' @param x a vector object containing the time-series of observations that are assumed to 
#'          be realizations of the (hidden Markov state dependent) observation process of the model.
#' @param m integer; a (finite) number of states in the hidden Markov chain.
#' @param delta  vector object containing starting values for the marginal probability 
#'         distribution of the \code{m} states of the Markov chain at the time 
#'         point \code{t=1} for the Baum-Welch algorithm.
#' @param gamma a matrix (\code{ncol=nrow=m}) containing starting values for the 
#'         transition matrix of the hidden Markov chain
#' @param distribution_class a single character string object with the abbreviated name of 
#'         the \code{m} observation distributions of the Markov dependent observation process. 
#'         The following distributions are supported: Poisson (\code{pois}); 
#'         generalized Poisson (\code{genpois}, parameter estimation via the 
#'         Baum-Welch algorithm is only supported if the M-step is performed numerically, 
#'         i.e. if \code{Mstep_numerical = TRUE}); normal (\code{norm})
#' @param distribution_theta a list object containing starting values for the parameters 
#'         of the \code{m} observation distributions of the observation process that are 
#'         dependent on the hidden Markov state.
#' @param discr_logL a logical object indicating whether the discrete log-likelihood 
#'        should be used (for \code{distribution_class="norm"}) for estimating the model 
#'        specific parameters instead of the general log-likelihood. 
#'        See MacDonald & Zucchini (2009, Paragraph 1.2.3) for further details.  
#'        Default value is \code{FALSE}.
#' @param discr_logL_eps a single numerical value to approximately determine the discrete 
#'        likelihood for a hidden Markov model based on nomal distributions 
#'        (for \code{"norm"}).  Default value is \code{0.5}.  
#'        See MacDonald & Zucchini (2009, Paragraph 1.2.3) for further details.
#' @param BW_max_iter a single numerical value representing the maximum number of iterations 
#'        in the Baum-Welch algorithm. Default value is \code{50}.
#' @param BW_limit_accuracy a single numerical value representing the convergence criterion 
#'        of the Baum-Welch algorithm. Default value is \code{0.001}.
#' @param BW_print a logical object indicating whether the log-likelihood at each 
#'        iteration-step shall be printed. Default value is \code{TRUE}.
#' @param Mstep_numerical a logical object indicating whether the Maximization Step of the 
#'        Baum-Welch algorithm shall be performed by numerical maximization using the 
#'        \link[stats]{nlm}-function.  Default value is \code{FALSE}.
#' @param DNM_limit_accuracy a single numerical value representing the convergence 
#'        criterion of the numerical maximization algorithm using the 
#'        \link[stats]{nlm}-function (used to perform the M-step of the 
#'        Baum-Welch-algorithm). Default value is \code{0.001}.
#' @param DNM_max_iter a single numerical value representing the maximum number of iterations
#'        of the numerical maximization using the \link[stats]{nlm}-function 
#'        (used to perform the M-step of the Baum-Welch-algorithm). Default value is \code{50}.
#' @param DNM_print a single numerical value to determine the level of printing of 
#'        the \code{nlm}-function.  See \code{nlm}-function for further informations. 
#'        The value \code{0} suppresses, that no printing will be outputted. 
#'        Default value is \code{2} for full printing.
#'
#' @return \code{Baum_Welch_algorithm} returns a list containing the estimated parameters 
#' of the hidden Markov model and other components. See MacDonald & Zucchini (2009, Paragraph 4.2) 
#' for further details on the calculated objects within this algorithm. 
#' \describe{
#' \item{x}{input time-series of observations.}
#' \item{m}{input number of hidden states in the Markov chain.}
#' \item{zeta}{a (T,m)-matrix (when T indicates the length/size of the observation 
#'      time-series and m the number of states of the HMM) containing probabilities 
#'      (estimates of the conditional expectations of the missing data given the 
#'      observations and the estimated model specific parameters) calculated by the algorithm. 
#'      See MacDonald & Zucchini (2009, Paragraph 4.2.2) for further details.}
#' \item{eta}{a (T,m,m)-dimensional-array (when T indicates the length of the observation 
#'      time-series and m the number of states of the HMM) containing probabilities 
#'      (estimates of the conditional expectations of the missing data given the 
#'      observations and the estimated model specific parameters) calculated by the algorithm. 
#'      See MacDonald & Zucchini (2009, Paragraph 4.2.2) for further details.}
#' \item{logL}{a numerical value representing the logarithmized likelihood calculated by 
#'       the \code{\link{forward_backward_algorithm}}.}
#' \item{iter}{number of performed iterations.}
#' \item{BIC}{a numerical value representing the Bayesian information criterion for the 
#'            hidden Markov model with estimated parameters.}
#' \item{delta}{a vector object containing the estimates for the marginal probability 
#'             distribution of the \code{m} states of the Markov chain at 
#'             time-point point \code{t=1}.}
#' \item{gamma}{a matrix containing the estimates for the transition matrix of the 
#'              hidden Markov chain.}
#' \item{...}{other input values (as arguments above). In the case that the algorithm 
#'            stops before the targeted accuracy or the maximum number of iterations has 
#'            been reached, further values are displayed and the estimates from the last 
#'            successful iteration step are saved.}}
#' @references Baum, L., Petrie, T., Soules, G., Weiss, N. (1970). A maximization technique 
#' occurring in the statistical analysis of probabilistic functions of markov chains. 
#' The annals of mathematical statistics, vol. \bold{41}(1), 164--171.
#' 
#' Dempster, A., Laird, N., Rubin, D. (1977). Maximum likelihood from incomplete data 
#' via the EM algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 
#' vol. \bold{39}(1), 1--38.
#' 
#' MacDonald, I. L.,  Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: 
#' An Introduction Using R}, Boca Raton: Chapman & Hall.     
#'
#' @author The basic algorithm for a Poisson-HMM is provided by MacDonald & Zucchini 
#'         (2009, Paragraph 4.2, Paragraph A.2.3).  Extension and implementation by 
#'         Vitali Witowski (2013).    
#'         
#' @seealso \code{\link{HMM_based_method}}, \code{\link{HMM_training}}, 
#' \code{\link{direct_numerical_maximization}},  \code{\link{forward_backward_algorithm}}, 
#' \code{\link{initial_parameter_training}}
#' 
#' @export
#' @keywords ts iteration 
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
#'   # Estimation of a HMM using the Baum-Welch algorithm
#'   \donttest{
#'     trained_HMM_with_m_hidden_states <- 
#'      Baum_Welch_algorithm(x = x, 
#'                           m = m, 
#'                       delta = delta, 
#'                       gamma = gamma,
#'          distribution_class = distribution_class, 
#'          distribution_theta = distribution_theta)
#'          
#'  print(trained_HMM_with_m_hidden_states)
#'  }
#'  
Baum_Welch_algorithm <- 
function(x, m, delta, gamma, distribution_class, distribution_theta, discr_logL = FALSE, 
    discr_logL_eps = 0.5, BW_max_iter = 50, BW_limit_accuracy = 0.001, BW_print=TRUE,
    Mstep_numerical = FALSE, DNM_limit_accuracy = 0.001, DNM_max_iter = 50, DNM_print = 2) 
{

  if (distribution_class == "pois" | distribution_class == "geom")
  {
    k=1
  }	
  if (distribution_class == "norm" | distribution_class == "genpois" |
      distribution_class == "bivariate_pois") 
  {
    k=2
  }	
  
  size <- length(x)	
  oldlogL <- -Inf
  reached_limit_of_accuracy <- FALSE
  underflow_error_1 <- FALSE
  underflow_error_2 <- FALSE
  pos_logL_error <- FALSE
  worse_logL_error <- FALSE

  
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
    wp = log(np / (1 - np) )
    return(wp)
  }
 
  DNM_invlogit_w2n <- function(wp)
  {
    np = exp(wp) /(1 + exp(wp))
    return(np)
  }
   
   
  negterm3 <- function(x, p, distribution_class, m, list_eta_zeta)
  { 
    size <- length(x)
    
    if (distribution_class == "pois") 
    {

      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dpois(x=x[tt], lambda=DNM_exp_w2n(p[i]))) 
        }
      }    
      
      if (is.na(term3) | term3 == Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }
      negterm3 <- -term3
    }
    
    
    if (distribution_class == "norm")
    {
      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dnorm(x[tt], mean=p[i], sd=DNM_exp_w2n(p[(m + i)])))
        }
      }
      if(is.na(term3) | term3 == Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }
      negterm3 <- -term3
    }
    
    if (distribution_class == "genpois") 
    {
      term3 <- 0
      for (i in 1:m)
      { 
        for (tt in 1:size)
        { 
          term3 <- term3 + list_eta_zeta$zeta[tt,i] * log( dgenpois(x[tt], lambda1=DNM_exp_w2n(p[i]), lambda2=DNM_invlogit_w2n(p[(m + i)])))
        }
      } 
      if(is.na(term3) | term3==Inf | term3 == -Inf)
      {
        term3 <- -Inf
      }  
      negterm3 <- -term3
    }	
    return(negterm3)
  }
  

# Training --------------------------------------------------------------------------
  for (l in 1:BW_max_iter)
  {
    
## Estimation Step (E-Step) ---------------------------------------------------------
   fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, 
            distribution_class = distribution_class, 
            distribution_theta = distribution_theta, 
            discr_logL = discr_logL, discr_logL_eps = discr_logL_eps), 
            silent = FALSE)    

    if (inherits(fb, "try-error")) 
    {  
       underflow_error_1 <- TRUE
       l <- l - 1
       logL <- oldlogL
       break
    }
    
    if (any(is.na(fb$log_alpha)) | any(is.na(fb$log_beta)))
    { 
      underflow_error_2 <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }
    
    
    logL <- fb$logL
    
    
    if (BW_print == TRUE) 
    { 
      print(list("##############################", 
            m = paste("train (EM) HMM with m =", toString(m)), 
            distribution = distribution_class, iteration = l, logL = logL, 
            "##############################"))
    }
    
    if (logL >= 0 ) 
    { 
      pos_logL_error <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }
    
    if (oldlogL > logL )     
    { 
      worse_logL_error <- TRUE
      l <- l - 1
      logL <- oldlogL
      break
    }   


    if (distribution_class == "pois")
    { 	
        zeta <- matrix(c(0), ncol=m, nrow = size)  
        zeta <- exp(fb$log_alpha + fb$log_beta - logL) 

        eta <- array(NA, dim = c((size-1), m, m))
        for (j in 1:m)
        { 
        	for (i in 1:m)
            { 
            	for (t in 2:size)
            	{ 
          			eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
          			    log( dpois(x[t], distribution_theta$lambda[j]) ) + fb$log_beta[t,j] - logL)
      	  		}
      		}
      	}    
      	      	 
    }
    
    
    
    if (distribution_class == "geom") 
    { 
        zeta <- matrix(c(0), ncol = m, nrow = size)  
        zeta <- exp( fb$log_alpha + fb$log_beta - logL) 

      eta <- array(NA, dim = c( (size - 1), m, m))
      for (j in 1:m) 
      	{ 
      		for (i in 1:m) 
      		{ 
      			for(t in 2:size)
      			{ 
      				eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
      				    log( dgeom(x[t], distribution_theta$prob[j]) ) + fb$log_beta[t,j] - logL)
      			}
      		}
      	}
      	
    }	
    
    
        
    if (distribution_class == "genpois")
    { 
      zeta <- matrix(c(0), ncol = m, nrow = size)  
      zeta <- exp( fb$log_alpha + fb$log_beta - logL) 
     
     	eta <- array(NA, dim = c((size - 1), m, m))
      	for (j in 1:m) 
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t-1, i, j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) 
      				    + log(dgenpois(x[t], distribution_theta$lambda1[j], distribution_theta$lambda2[j]) ) 
      				    + fb$log_beta[t,j] - logL)
      			}
      		}
      	}
      	
    }
    
    
    if (distribution_class == "norm")
    { 
      zeta <- matrix(c(0), ncol = m, nrow = size )  
      zeta <- exp( fb$log_alpha + fb$log_beta - logL) 
      
    	eta <- array(NA, dim = c((size-1), m, m))
      	for (j in 1:m)
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t-1,i , j] <- exp( fb$log_alpha[(t-1),i] + log(gamma[i,j]) + 
      				    log(dnorm(x[t], distribution_theta$mean[j], distribution_theta$sd[j]) ) 
      				    + fb$log_beta[t,j]  - logL)
      			}
      		}
      	}
      	
    }
    
    
    
    if (distribution_class == "bivariate_pois")
    { 	
      size=length(x[,1])
       
      zeta <- matrix(c(0), ncol = m, nrow=size)  
      zeta <- exp(fb$log_alpha + fb$log_beta - logL) 
      
      eta <- array(NA, dim = c((size - 1), m, m))
      	for (j in 1:m)
      	{ 
      		for (i in 1:m)
      		{ 
      			for (t in 2:size)
      			{ 
      				eta[t - 1, i, j] <- exp( fb$log_alpha[(t - 1),i] + log(gamma[i, j]) 
      				  + log(dpois(x[t, 1], distribution_theta$lambda_1[j]) ) 
      				  + log(dpois(x[t, 2], distribution_theta$lambda_2[j])) 
      				  + fb$log_beta[t, j] - logL)          
      			}
      		}
      	}    	
    }


  list_eta_zeta = list(zeta = zeta, eta = eta)
          

## Maximization Step (M-Step) ------------------------------------------------
    
 delta <- list_eta_zeta$zeta[1,]
    
 gamma <- matrix(0, ncol=m, nrow=m)
    for (i in 1:m)
    { 
    	for (j in 1:m)
    	{ 
    		sum_numerator <- 0
    		sum_numerator <- sum(list_eta_zeta$eta[, i, j]) 
    		sum_denominator <- 0
      		for (k in 1:m) 
      		{ 
        		sum_denominator <- sum_denominator + sum(list_eta_zeta$eta[, i, k]) 
      		}
      
      	gamma[i, j] <-  sum_numerator / sum_denominator      
    	}
    }
    
    
    
    if (distribution_class == "pois" & Mstep_numerical == FALSE) 
    { 
      
      lambda <- rep(0, times = m)
      
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[, i] * x)
        sum_denominator <- 0
        sum_denominator <- sum(list_eta_zeta$zeta[, i])
                
        lambda[i] <- sum_numerator / sum_denominator
      }
      
      estimated_mean_values <- lambda
      distribution_theta <- list(lambda = lambda)
    }
    
    
    if (distribution_class == "norm" & Mstep_numerical == FALSE)
    { 	
      mean <- rep(0, times = m)
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[, i] * x)
          
        sum_denominator <- 0
        sum_denominator <- sum(list_eta_zeta$zeta[, i])
     
      mean[i] <- sum_numerator / sum_denominator
      }
      
      sd <- rep(0, times = m)
      for (i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt,i] * ((x[tt] - mean[i])^2)
        }
                  
       sum_denominator <- 0
       sum_denominator <- sum(list_eta_zeta$zeta[, i])
       sd[i] <- sum_numerator / sum_denominator
       sd[i] <- sqrt(sd[i])
      }
      
      estimated_mean_values <- mean
      distribution_theta <- list(mean = mean, sd = sd)
    }
    
    
    if(distribution_class == "geom") 
    { 
      prob <- rep(0, times = m)
      for (i in 1:m) 
      { 
        sum_numerator <- 0
        sum_numerator <- sum(list_eta_zeta$zeta[,i] * x)
        

        sum_denominator <- 0
        for (tt in 1:size)
        {
          sum_denominator <- sum_denominator + list_eta_zeta$zeta[tt, i] * (x[tt] - 1)
        }
        
        prob[i] <- sum_numerator / sum_denominator
      }
      
      estimated_mean_values <- 1 / prob
      distribution_theta <- list(prob = prob)
    }
    
    
    if(distribution_class == "bivariate_pois")
    { 
      size <- length(x[, 1])	
      lambda_1 <- rep(0, times = m)
      lambda_2 <- rep(0, times = m)
      
      for(i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt, i] * x[tt, 1]
        }
    
        sum_denominator <- 0
        for (tt in 1:size)
        { 
          sum_denominator + list_eta_zeta$zeta[tt, i]
        }
        
        lambda_1[i] <- sum_numerator / sum_denominator
      }
      
      for (i in 1:m)
      { 
      	sum_numerator <- 0
        for (tt in 1:size)
        { 
        	sum_numerator <- sum_numerator + list_eta_zeta$zeta[tt, i] * x[tt, 2]
        }
        
        sum_denominator <- 0
        for (tt in 1:size)
        { 
          sum_denominator <- sum_denominator + list_eta_zeta$zeta[tt, i]
        }
        
        lambda_2[i] <- sum_numerator / sum_denominator
      }
      
      
      estimated_mean_values <- lambda_1
      
      distribution_theta <- list(lambda_1 = lambda_1, lambda_2 = lambda_2)
    }
    
 
    
    if(distribution_class == "pois" & Mstep_numerical == TRUE) 
    {
      trans_lambda <- DNM_log_n2w(distribution_theta$lambda)
      
      vector_of_parameters <- c(trans_lambda)   
      
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, print.level = DNM_print, 
                      gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
            
      est_lambda <- DNM_exp_w2n(minterm3$estimate)
      
      estimated_mean_values <- est_lambda 
      
      estimated_var <- est_lambda
      estimated_sd <- sqrt(estimated_var)
      
      distribution_theta <- list(lambda=est_lambda)
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))  
    }
    
    
    if (distribution_class == "norm" & Mstep_numerical == TRUE)
    {	
      trans_mean <- distribution_theta$mean
      
      trans_sd <- DNM_log_n2w(distribution_theta$sd)   
      
      vector_of_parameters <- c(trans_mean, trans_sd) 
        
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, print.level = DNM_print, 
                      gradtol = DNM_limit_accuracy, iterlim = DNM_max_iter)
      
      estimated_mean <- minterm3$estimate[1:m]
      
      sum_denominator <- minterm3$estimate[1:m]
      
      estimated_sd <- DNM_exp_w2n(minterm3$estimate[(m + 1):(2 * m)])
      
      estimated_mean_values <- estimated_mean
      
      estimated_var <- estimated_sd^2
      
      estimated_sd <- estimated_sd      
      
      distribution_theta <- list(mean = estimated_mean, sd = estimated_sd) 
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))
    }
    
    if (distribution_class == "genpois" & Mstep_numerical == TRUE) 
    {	
      trans_lambda1 <- DNM_log_n2w(distribution_theta$lambda1)
      
      trans_lambda2 <- DNM_logit_n2w(distribution_theta$lambda2)
      
      vector_of_parameters <- c(trans_lambda1, trans_lambda2)   
      
      minterm3 <- nlm(negterm3, distribution_class = distribution_class, m = m, 
                      p = vector_of_parameters, x = x, 
                      list_eta_zeta = list_eta_zeta, 
                      print.level = DNM_print, gradtol = DNM_limit_accuracy, 
                      iterlim = DNM_max_iter)
      
      estimated_lambda1 <- DNM_exp_w2n(minterm3$estimate[1:m])
      
      estimated_lambda2 <- DNM_invlogit_w2n(minterm3$estimate[(m + 1):(m + m)])
      
      estimated_mean_values <- estimated_lambda1 / (1 - estimated_lambda2)
      
      estimated_var <- estimated_lambda1 / ((1-estimated_lambda2)^3)
      
      estimated_sd <- sqrt(estimated_var)
      
      distribution_theta <- list(lambda1 = estimated_lambda1, lambda2 = estimated_lambda2)
      
      print(list(distribution_theta = distribution_theta, 
                 E = estimated_mean_values, 
                 var = estimated_var, 
                 sd = estimated_sd))
    }
    
    difference_old_logL_and_new_logL = abs(oldlogL - logL)
    
    if (difference_old_logL_and_new_logL < BW_limit_accuracy) 
    { 
      reached_limit_of_accuracy = TRUE
      break
    }
    
    oldlogL <- logL
  }


# Accessing AIC and BIC for the trained HMM  ----------------------------------------
  AIC <- AIC_HMM(logL = logL, m = m, k = k) 
  BIC <- BIC_HMM(size = size, logL = logL, m = m, k = k) 
  

# Results ---------------------------------------------------------------------------
  return(list(x = x,
              m = m,
              zeta = zeta,
              eta = eta,
              iter = l,
              logL = logL,
              AIC = AIC,
              BIC = BIC,
              delta = delta,
              gamma = gamma,
              distribution_class = distribution_class,
              distribution_theta = distribution_theta,
              estimated_mean_values = estimated_mean_values,
              reached_limit_of_accuracy = reached_limit_of_accuracy,
              underflow_error_1 = underflow_error_1,
              underflow_error_2 = underflow_error_2,
              pos_logL_error = pos_logL_error,
              worse_logL_error = worse_logL_error))
}
