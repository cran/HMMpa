local_decoding_algorithm <-
function(x, m, delta, gamma, distribution_class, distribution_theta, discr_logL = FALSE, discr_logL_eps = 0.5)
{

################################################################################
######################### Needed variables and functions #######################
################################################################################
  svtpoue <- small_value_to_prevent_overflow_and_underflow_errors <- 4.940656e-142
  
  size <- length(x)
      
################################################################################
#### Preparation for calcualtion of local most probable states #################
################################################################################  
fb <-  try( forward_backward_algorithm(x = x, gamma = gamma, delta = delta, 
         distribution_class = distribution_class, 
         distribution_theta = distribution_theta, discr_logL = discr_logL, 
         discr_logL_eps = discr_logL_eps), silent = FALSE)
  if (inherits(fb,"try-error"))
  {  
    stop("underflow-error inside function forward_backward_algorithm ")
  }
  
  if (distribution_class == "pois")
  {
    distribution_means <- distribution_theta$lambda		
  }
  
  if (distribution_class == "geom")
  {
    distribution_means <- 1 / distribution_theta$prob$prob					
  }
  
  if (distribution_class == "genpois")
  {	
    distribution_means <- distribution_theta$lambda1 / (1 - distribution_theta$lambda2) 
  }
  
  if (distribution_class == "norm" & discr_logL == FALSE)
  {	
    distribution_means <- distribution_theta$mean
  }
  
  if (distribution_class == "norm" & discr_logL == TRUE)
  {	
    distribution_means <- distribution_theta$mean
  }

################################################################################
################ Calculation of probabilities P(Z_t=i | X_1,...,X_t)  ##########
################################################################################  
  
  state_probabilities <- matrix(NA, ncol = m, nrow = size)
  for (i in 1:size)
  {
  	state_probabilities[i,] <- exp(fb$log_alpha[i,] + fb$log_beta[i,] - fb$logL)
  }

################################################################################
######################### The algorithm     ####################################
################################################################################
  
 decoding <- rep(NA,size)
 for (i in 1:size)
 {
 	decoding[i] <- which.max(state_probabilities[i,])
 }
 
################################################################################
############## Return results ##################################################
################################################################################
  
return(list(state_probabilities = state_probabilities, 
	          decoding = decoding))
}
