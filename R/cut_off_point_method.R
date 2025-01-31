#' Cut-Off Point Method for Assigning Physical Activity Patterns
#'
#' This function assigns an activity range to each observation of a time-series, 
#' such as for a sequence of impulse counts recorded by an accelerometer.  
#' The activity ranges are defined by thresholds called \dQuote{cut-off points}.  
#' Furthermore, bout periods are analysed (see Details for further informations). 
#' 
#' @param x a vector object of length \code{T} containing non-negative observations of a 
#'          time-series, such as a sequence of accelerometer impulse counts.
#' @param cut_points a vector object containing cut-off points to separate activity ranges. 
#'          For instance, the vector c(7,15,23) separates the four activity ranges 
#'          [0,7);[7,15);[15,23);[23,Inf).
#' @param names_activity_ranges an optional character string vector to name the activity 
#'          ranges induced by the cut-points.  This vector must contain one element more 
#'          than the vector \code{cut_points}.
#' @param hidden_PA_levels an optional vector object of length \code{T} containing a 
#'    sequence of the estimated hidden physical activity levels (i.e. means) underlying 
#'    the time-series of accelerometer counts.  Such a sequence can be extracted by 
#'    decoding a trained hidden Markov model.  The cut-point method classifies then each 
#'    count by its level in the hidden Markov chain that generates the physical activity 
#'    counts, and does not use the observed count value (see \code{\link{HMM_based_method}} 
#'    for further details).  Default is NA (for the traditional cut-point method).
#' @param bout_lengths a vector object (with even number of elemets) to define the 
#'    range of the bout intervals (see Details for the definition of bouts). For instance, 
#'    \code{bout_lengths = c(1,1,2,2,3,10,11,20,1,20)} defines the five bout intervals 
#'    [1,1] (1 count); [2,2] (2 counts); [3,10] (3-10 counts); [11,20] (11-20 counts); [1,20] 
#'    (1-20 counts - overlapping with other bout intervalls is possible). 
#'    Default value is \code{bout_lengths=NULL}.
#' @param plotting a numeric value between \code{0} and \code{5} (generates different outputs). 
#'    NA suppresses graphical output. Default value is \code{0}.\cr
#'    \code{0}: output 1-5 \cr
#'    \code{1}: summary of all results \cr
#'    \code{2}: time series of activity counts, classified into activity ranges \cr
#'    \code{3}: time series of bouts (and, if available, the sequence of the estimated 
#'    hidden physical activity levels, extracted by decoding a trained HMM, in green color)\cr
#'    \code{4}: barplots of absolute and relative frequencies of time spent in different 
#'    activity ranges\cr
#'    \code{5}: barplots of absolute frequencies of different bout intervals 
#'    (overall and by activity ranges )
#'    
#' @details
#' A bout is defined as a period of time spending a defined intensity of physical 
#' activities in a specified physical activity range, without switching to activity 
#' intensities in a different activity range.
#' 
#'
#' @return \code{ cut_off_point_method } returns a list containing the extracted sequence 
#' of activity ranges  and plots key figures. 
#' \describe{
#' \item{activity_ranges}{an array object containing the cut-off intervals that indicate the 
#'       activity ranges.}
#' \item{classification}{an integer vector containing the sequence of activity ranges 
#'      that were assigned to the observed time-series of accelerometer counts. 
#'      If \code{hidden_PA_levels=NA}, then \code{classification} is the output of the 
#'      traditional cut-point method, meaning that an activity range has been assigned to 
#'      each accelerometer count over its observed value actual position.  In case when 
#'      \code{hidden_PA_levels} is available, \code{classification} is the output of the 
#'      extendend cut-point method using hidden Markov models (see \code{\link{HMM_based_method}} 
#'      for further details). }
#' \item{classification_per_activity_range}{a pairlist object containing the classification 
#'      of the observed counts by the assigned activity range.}
#' \item{freq_acitvity_range}{table object containing the absolute frequencies of 
#'      classifications into activity ranges.}
#' \item{rel_freq_acitvity_range}{table object containing the relative frequencies of 
#'      classifications into activity ranges.}
#' \item{quantity_of_bouts}{overall number of bouts.}
#' \item{bout_periods}{an array including the bout length assigned to acitiy ranges.}
#' \item{abs_freq_bouts_el}{a pairlist object containing the absolute frequency of bout 
#'       length per epoch length (aggregated).}}
#' @author Vitali Witowski (2013)
#' @seealso \link{HMM_based_method}
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
#' # 1) Traditional cut point method -----------------------
#' # Assigning activity ranges to activity counts using 
#' # fictitious cut-off points that produce the four activity 
#' # ranges "sedentary"", "light"", "moderate"", and "vigorous". 
#' \donttest{
#'   solution_of_traditional_cut_off_point_method <- 
#'       cut_off_point_method(x = x, 
#'                            cut_points = c(5,15,23), 
#'                            names_activity_ranges = c("SED","LIG","MOD","VIG"), 
#'                            bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12,13,40,41,265,1,265), 
#'                            plotting = 0)
#'  print(solution_of_traditional_cut_off_point_method)
#'  }
#'  
#' # 2) Extension of the traditional cut_point method 
#' #     using HMMs                                     
#' # The following three steps define an extension of the 
#' # traditional cut-off method by first extracting the hidden 
#' # physical activity pattern behind the accelerometer counts 
#' # using a HMM (those three steps are basically combined in 
#' # the function HMM_based_method, see HMM_based_method for 
#' # further details and references): 
#' 
#' # Step 1 ---
#' # Train hidden Markov model for different number of 
#' # states m=2,...,6 and select the model with the most 
#' # plausible m
#'  \donttest{
#' m_trained_HMM <- 
#'   HMM_training(x = x, 
#'                min_m = 2, 
#'                max_m = 6, BW_print=FALSE,
#'                distribution_class = "pois")$trained_HMM_with_selected_m
#'  }
#'  
#'  # Step 2 ---
#'  # Decode the trained HMM (by using the 
#'  # Viterbi algorithm (global decoding)) to get the estimated 
#'  # sequence of hidden physical activity levels 
#'  # underlying the the accelerometer counts 
#'  
#'  # You have to compute 'm_trained_HMM' first (see Step 1)
#'  \donttest{
#'    global_decoding <- 
#'      HMM_decoding(x = x, 
#'                   m = m_trained_HMM$m, 
#'                   delta = m_trained_HMM$delta, 
#'                   gamma = m_trained_HMM$gamma, 
#'                   distribution_class = m_trained_HMM$distribution_class, 
#'                   distribution_theta = m_trained_HMM$distribution_theta,
#'                   decoding_method = "global")
#'   hidden_PA_levels <- global_decoding$decoding_distr_means
#'   }
#'   
#'   # Step 3 ---
#'   # Assigning activity ranges to activity counts using the 
#'   # information extracted by decoding the HMM for the counts 
#'   # (PA-levels) and fictitious cut-off points that produce 
#'   # four so-called activity ranges:"sedentary", "light", 
#'   # "moderate" and "vigorous":
#'   
#'   # You have to compute 'm_trained_HMM' and 'hidden_PA_levels' first (see above)
#'   \donttest{
#'     solution_of_HMM_based_cut_off_point_method <- 
#'         cut_off_point_method(x = x, 
#'                              hidden_PA_levels = hidden_PA_levels, 
#'                              cut_points = c(5,15,23), 
#'                              names_activity_ranges = c("SED","LIG","MOD","VIG"), 
#'                              bout_lengths = c(1,1,2,2,3,3,4,4,5,5,6,12,13,40,41,265,1,265), 
#'                              plotting=1)
#'  }

cut_off_point_method <-
function(x, cut_points, names_activity_ranges = NA, hidden_PA_levels = NA, bout_lengths = NULL, plotting= 0) 
{   

 if (is.null(bout_lengths))
 {
 	stop("Set variable 'bout_lengths' to use this function. See help-manual for further information. For example: bout_lengths=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,12,13,20,21,40,41,60,61,80,81,120,121,240,241,480,481,1440,1,1440)")
 }
	
  cut_points <- c(0, cut_points)
  number_of_activity_ranges <- length(cut_points)
  
  activity_ranges <- array(data = NA, dim = c(1,2,number_of_activity_ranges))
  colnames(activity_ranges) <- c("[lower boundary,","upper boundary)")
  
  # Function for bout-extraction (from the result of the classification) ----------------

  bouts <- function(counts, sequence_activity_range, cut_points, bout_lengthss)
  {
    
    for_plotting_a <- sequence_activity_range
    
    for (i in 1:length(cut_points))
    {
      for_plotting_a[for_plotting_a == i] <- cut_points[i]
    }
    for_plotting_a[sequence_activity_range == (length(cut_points) + 1)] <- max(counts)
    
    
    for_plotting_b <- sequence_activity_range
    for_plotting_b <- for_plotting_b - 1
    for (i in 1:length(cut_points))
    {
      for_plotting_b[for_plotting_b == i] <- cut_points[i]
    }
    for_plotting_b[for_plotting_b == 0] <- 0
    
    
    
    dimension_columns_of_next_arry <- max(length(bout_lengthss), length(x), max(bout_lengthss))
    if (dimension_columns_of_next_arry == Inf)
    {
      dimension_columns_of_next_arry <- max(bout_lengthss[!bout_lengthss == Inf], length(x))
    }
    
    bout_periods <- array(0, dim = c((length(cut_points) + 1), dimension_columns_of_next_arry, 1))
    
    if (all(!is.na(names_activity_ranges))) 
    {
      dimnames(bout_periods) <- list(activity_range = names_activity_ranges, 
                                     Countbouts = 1: dimension_columns_of_next_arry)
    } else {
    dimnames(bout_periods) <- list(activity_range = 1:(length(cut_points) + 1), 
                                   Countbouts = 1: dimension_columns_of_next_arry)
    }
    temp_sequence_activity_range <- c(sequence_activity_range, -Inf)
    
    for (j in 1:(length(cut_points) + 1))
    {
      temp_vec <- 0
      for (i in 1:length(temp_sequence_activity_range)) 
      {
        if (temp_sequence_activity_range[i] == j) 
        {
          temp_vec <- temp_vec + 1
        } else {
          temp_which_bout <- temp_vec
          bout_periods[j,temp_which_bout,1] <- bout_periods[j,temp_which_bout,1] + 1
          
          temp_vec <- 0  
        }
        
        if (temp_vec == length(temp_sequence_activity_range))
        {
          temp_which_bout <- temp_vec
          bout_periods[j,temp_which_bout,1] = bout_periods[j,temp_which_bout,1] + 1
          
          temp_vec <- 0
        }
      }
    }
     
    quantity_of_bouts <- sum(bout_periods)
    
    return(list(quantity_of_bouts = quantity_of_bouts,
                bout_periods = bout_periods,
                for_plotting_a = for_plotting_a, 
                for_plotting_b = for_plotting_b)) 
  }


  for (i in 1:number_of_activity_ranges)
  {
    activity_ranges[1,1,i] <- cut_points[i]
    activity_ranges[1,2,i] <- cut_points[i+1]
  }
  activity_ranges[1,2,number_of_activity_ranges] <- Inf
  
  classification <- rep(NA, times = length(x))
  
  if (any(is.na(hidden_PA_levels)))
  {
    for (i in 1:length(x))
    {	
      for (j in 1:number_of_activity_ranges)
      {
        if (x[i] >= activity_ranges[1,1,j] & x[i] < activity_ranges[1,2,j])
        {	
          classification[i] <- c(j)
        }
      }
    }
    classification_per_activity_range <- pairlist(NA)
    for (i in 1:number_of_activity_ranges)
    {
      classification_per_activity_range[[i]] <- x
      for (j in 1:length(classification_per_activity_range[[i]]))
      {
        if (!classification[j] == i) 
        {
          classification_per_activity_range[[i]][j] <- NA
        }
      }
    }
  }
  # If (hidden) Pa-levels are given, then the classification is conducted via the 
  # total magnitude of the hidden PA-level of a count 
  if (any(!is.na(hidden_PA_levels)))
  {   
    for (i in 1:length(x))
    {	
      for (j in 1:number_of_activity_ranges)
      {
        if (hidden_PA_levels[i] >= activity_ranges[1,1,j] & hidden_PA_levels[i] < activity_ranges[1,2,j])
        {	
          classification[i] <- c(j)
        }
      }
    }
    classification_per_activity_range <- pairlist(NA)
    for (i in 1:number_of_activity_ranges)
    {
      classification_per_activity_range[[i]] <- x
      for (j in 1:length(classification_per_activity_range[[i]]))
      {
        if (!classification[j] == i)
        {
          classification_per_activity_range[[i]][j] <- NA
        }
      }
    }
  }


  freq_acitvity_range <- table(c(seq(1, number_of_activity_ranges), classification))
  freq_acitvity_range <- freq_acitvity_range - 1
  rel_freq_acitvity_range <- freq_acitvity_range / sum(table(classification))
  
  if (all(!is.na(names_activity_ranges)))
  {
    rownames(freq_acitvity_range) <- names_activity_ranges
    rownames(rel_freq_acitvity_range) <- names_activity_ranges
    colnames(activity_ranges[1,,]) <- names_activity_ranges
    names(classification_per_activity_range ) <- names_activity_ranges
  }
  
  # Bout-extraction (from the result of the classification) ----------------------------

  bout_extraction <- bouts(counts = x, sequence_activity_range = classification, 
                           cut_points = cut_points[c(-1)], bout_lengthss = bout_lengths)
  
  quantity_of_bouts <- bout_extraction$quantity_of_bouts
  
  bout_periods <- bout_extraction$bout_periods
  
  
  bout_classes <- matrix(bout_lengths, nrow = 2, byrow = FALSE)
  names = NULL
  for (i in seq(1, dim(bout_classes)[2]))
  {
    if (bout_classes[1,i] == bout_classes[2,i])
    {
      names <- c(names, toString(bout_classes[1,i]))
    } else {
      names <- c(names, paste(toString(bout_classes[1,i]), '-', toString(bout_classes[2,i])) )
    }
  }

  colnames(bout_classes) <- names
  rownames(bout_classes) <- c("from", "to")
  
  
  abs_freq_bouts_el <- rep(0, times = length(names))
  for (i in 1:(length(names)))
  {
    abs_freq_bouts_el[i] <- sum(bout_periods[,seq(bout_classes[1,i], bout_classes[2,i]),1])
  }
  

  
  names(abs_freq_bouts_el) <- names
  
  abs_freq_bouts_el <- pairlist(all_activity_ranges = abs_freq_bouts_el)
  
  for (j in 1: number_of_activity_ranges)	    
  {	 
    abs <- rep(0, times = length(names))
    for (i in 1:(length(names)))
    {
      abs[i] <- sum(bout_periods[j,seq(bout_classes[1,i], bout_classes[2,i]),1])
    }
      
    names(abs) <- names
    
    if (all(!is.na(names_activity_ranges))) 
    {
      abs_freq_bouts_el[[names_activity_ranges[j]]] <- abs	
    } else {
      abs_freq_bouts_el[[j+1]] <- abs	
    }
  }
  

  if (!is.na(plotting))
  {	
  	# output plot 1-5 
    if (plotting == 0)
    {   
      par(mfrow = c(2,2))
      
      if (any(!is.na(hidden_PA_levels))) 
      {
      	plot(x, main = "classification (HMM)", xlab = "time", ylab= "counts", col = "white")
      } else {
         plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	
      
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i], '(', round(100 * rel_freq_acitvity_range[i], digits = 2), "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i, '(', round(100 * rel_freq_acitvity_range[i], digits=2), "%",')')))
        }
      }	
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges (HMM)"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      } else {
      		barplot(rel_freq_acitvity_range, main = c("fragmentation into acitvity ranges"), xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }	
  
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
   	  		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }		      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      if (any(!is.na(hidden_PA_levels)))
      {
     		 barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), 
     		         ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[1]], 
     		         names.arg = names, main = "bouts in relation to their length (HMM)")      
      } else {
     		 barplot(ylab = ("abs. frequency"), xlab = c("length of bout [epoch period]"), 
     		         ylim = c(0, max(abs_freq_bouts_el[[1]])), abs_freq_bouts_el[[1]], 
     		         names.arg = names, main = "bouts in relation to their length")      
      }		  
      

   
      par(mfrow = c(1,1))
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }		 
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, ylim = c(0, max(x) + max(x) * 0.1), main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }
      
          
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      
      par(mfrow = c(2,1))
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i], 
                                '(', round(100 * rel_freq_acitvity_range[i], digits=2), 
                                "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '(', round(100 * rel_freq_acitvity_range[i], digits=2), 
                                "%", ')')))
        }
      }		
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    	barplot(rel_freq_acitvity_range, 
      	    	        main = c("fragmentation into acitvity ranges (HMM)"), 
      	    	        xlab = "activity range", 
      	    	        ylab = "relative frequency", 
      	    	        names.arg = foo)
      } else {
      			barplot(rel_freq_acitvity_range,
      			        main = c("fragmentation into acitvity ranges"), 
      			        xlab = "activity range", 
      			        ylab = "relative frequency", 
      			        names.arg = foo)
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i], 
                                '( # ', round(freq_acitvity_range[i], digits=2), ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , '( # ', round( freq_acitvity_range[i], digits=2), ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
   	  		barplot(freq_acitvity_range, 
   	  		        main = c("fragmentation into acitvity ranges (HMM)"), 
   	  		        xlab = "activity range", 
   	  		        ylab = "absolute frequency", 
   	  		        names.arg = foo)	
      } else {
      		barplot(freq_acitvity_range, 
      		        main = c("fragmentation into acitvity ranges"), 
      		        xlab = "activity range", 
      		        ylab = "absolute frequency", 
      		        names.arg = foo)		
      }		
      
      
      par(mfrow = c(number_of_activity_ranges / 2 + 1, 2))
      
      for (i in 1:(number_of_activity_ranges + 1))
      {   
      	
      	if (any(!is.na(hidden_PA_levels)))
      	{
      		    barplot(ylab = ("abs. frequency"), 
      		            xlab = c("length of bout [epoch period]"), 
      		            ylim = c(0, max(abs_freq_bouts_el[[1]])), 
      		            abs_freq_bouts_el[[i]], names.arg = names, 
      		            main = "bouts in relation to their length (HMM)")
      	} else {
      			barplot(ylab = ("abs. frequency"), 
      			        xlab = c("length of bout [epoch period]"), 
      			        ylim = c(0, max(abs_freq_bouts_el[[1]])), 
      			        abs_freq_bouts_el[[i]], names.arg = names, 
      			        main = "bouts in relation to their length")
      	}
      

        if (all(!is.na(names_activity_ranges)))
        {
          legend("topright", bg = "white", 
                 c("all", names_activity_ranges)[i], fill = c("gray"))
        } else {
          legend("topright", bg = "white", 
                 c("all",seq(1, number_of_activity_ranges))[i], fill = c("gray"))	
        }
      }
      par(mfrow = c(1,1))  
    }
    
    
    # Summary of all results -------------------------------------------------------------
    if (plotting == 1)
    {   
      par(mfrow = c(2,2))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
         plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i], 
                                '(', round( 100* rel_freq_acitvity_range[i], digits = 2), 
                                "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, 
                   c(paste(i  , '(', round( 100* rel_freq_acitvity_range[i], digits = 2), 
                           "%", ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		barplot(rel_freq_acitvity_range, 
      		        main = c("fragmentation into acitvity ranges (HMM)"), 
      		        xlab = "activity range", 
      		        ylab = "relative frequency", 
      		        names.arg = foo)
      } else {
      		barplot(rel_freq_acitvity_range, 
      		        main = c("fragmentation into acitvity ranges"), 
      		        xlab = "activity range", 
      		        ylab = "relative frequency", 
      		        names.arg = foo)
      }		
      
      
      if (any(!is.na(hidden_PA_levels)))
      {
      		plot(x, main = c("# bouts (HMM)", quantity_of_bouts), xlab = "time", 
      		     ylab = "counts", col = "white")
      } else {
      		plot(x, main = c("# bouts", quantity_of_bouts), xlab = "time", ylab = "counts", col = "white")
      }		
      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      
      
            
      if (any(!is.na(hidden_PA_levels)))
      {
 			barplot(ylab = ("abs. frequency"), 
 			        xlab = c("length of bout [epoch period]"), 
 			        ylim = c(0, max(abs_freq_bouts_el[[1]])),
              abs_freq_bouts_el[[1]],
 			        names.arg = names, 
 			        main = "bouts in relation to their length (HMM)")
      } else {
 			barplot(ylab = ("abs. frequency"),
 			        xlab = c("length of bout [epoch period]"), 
 			        ylim = c(0, max(abs_freq_bouts_el[[1]])),
              abs_freq_bouts_el[[1]], 
 			        names.arg = names, 
 			        main = "bouts in relation to their length")
      }		    
      
      par(mfrow = c(1,1))    
    }	
    
    
    # time series of activity counts, classified into activity ranges ----------
    if (plotting == 2)
    {   
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	      plot(x, main = "classification (HMM)", xlab = "time", ylab = "counts", col = "white")
      } else {
      	plot(x, main = "classification", xlab = "time", ylab = "counts", col = "white")
      }	      
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
      for (i in 1: number_of_activity_ranges)
      {
        points(as.numeric(classification_per_activity_range[[i]]), col = i)	
      }
    }	
    
    
    # time series of bouts (and, if available, the sequence of the estimated 
    # hidden physical activity levels, extracted by decoding a trained HMM) 
    if (plotting == 3)
    {   
      par(mfrow = c(1,1))
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    plot(x, ylim = c(0, max(x) + max(x) * 0.1), 
      	         main = c("# bouts (HMM)", quantity_of_bouts), 
      	         xlab = "time", ylab = "counts", col = "white")
      } else {
      		plot(x, ylim = c(0, max(x) + max(x) * 0.1), 
      		     main = c("# bouts", quantity_of_bouts), 
      		     xlab = "time", ylab = "counts", col = "white")
      }
      
      
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_a, col = "black", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      points(bout_extraction$for_plotting_b, col = "white", type="h")
      
      if (any(!is.na(hidden_PA_levels)))
      {
        points(hidden_PA_levels, col = "green", type = "l")
      }
      
      abline(h = cut_points[c(-1)], col = "grey50", lty = "dashed")
    }
    
    
    # barplots of absolute and relative frequencies of time spent in different activity ranges 
    if (plotting == 4)
    {   
      par(mfrow = c(2,1))
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i] , 
                                '(', round( 100* rel_freq_acitvity_range[i], digits = 2), 
                                "%", ')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  , 
                                '(', round( 100* rel_freq_acitvity_range[i], digits = 2), 
                                "%", ')')))
        }
      }		
      
      if (any(!is.na(hidden_PA_levels)))
      {
      	    	barplot(rel_freq_acitvity_range, 
      	    	        main = c("fragmentation into acitvity ranges (HMM)"), 
      	    	        xlab = "activity range", 
      	    	        ylab = "relative frequency", 
      	    	        names.arg = foo)
      } else {
      			barplot(rel_freq_acitvity_range, 
      			        main = c("fragmentation into acitvity ranges"), 
      			        xlab = "activity range", ylab = "relative frequency", names.arg = foo)
      }
      
      
      if (all(!is.na(names_activity_ranges)))
      {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(names_activity_ranges[i],
                                '( # ',round(freq_acitvity_range[i], digits = 2),')')))
        }
      } else {
        foo <- NULL
        for (i in 1: number_of_activity_ranges)
        {
          foo <- c(foo, c(paste(i  ,'( # ',round( freq_acitvity_range[i], digits = 2),')')))
        }
      }		
      

      
      if (any(!is.na(hidden_PA_levels)))
      {
 			barplot(freq_acitvity_range, 
 			        main = c("fragmentation into acitvity ranges (HMM)"), 
 			        xlab = "activity range", 
 			        ylab = "absolute frequency", names.arg = foo)		
      } else {
 			barplot(freq_acitvity_range, 
 			        main = c("fragmentation into acitvity ranges"), 
 			        xlab = "activity range", 
 			        ylab = "absolute frequency", names.arg = foo)		
      }
      
      par(mfrow = c(1,1))  	
    }	
    
    
    # barplots of absolute frequencies of different bout intervals 
    # (overall and by activity ranges ) 
    if (plotting == 5)
    {   
      par(mfrow = c(number_of_activity_ranges/2+1,2))
      
      for (i in 1:(number_of_activity_ranges+1))
      {
      	
      if (any(!is.na(hidden_PA_levels)))
      {
      	    barplot(ylab = ("abs. frequency"), 
      	            xlab = c("length of bout [epoch period]"), 
      	            ylim = c(0, max(abs_freq_bouts_el[[1]])), 
      	            abs_freq_bouts_el[[i]], names.arg = names, 
      	            main = "bouts in relation to their length (HMM)")
      } else {
      		barplot(ylab = ("abs. frequency"), 
      		        xlab = c("length of bout [epoch period]"), 
      		        ylim = c(0, max(abs_freq_bouts_el[[1]])), 
      		        abs_freq_bouts_el[[i]], 
      		        names.arg = names, 
      		        main = "bouts in relation to their length")
      }
      
      
        if (all(!is.na(names_activity_ranges)))
        {
          legend("topright", bg = "white", c("all",names_activity_ranges)[i], fill =c("gray"))
        } else {
          legend("topright", bg = "white", c("all",seq(1, number_of_activity_ranges))[i], fill =c("gray"))	
        }
      }
      par(mfrow = c(1,1))  
    }
  }
  

print(list(activity_ranges = activity_ranges,
           classification = classification, 
           rel_freq_acitvity_range = rel_freq_acitvity_range,
           quantity_of_bouts = quantity_of_bouts, 
           abs_freq_bouts_el = abs_freq_bouts_el$all_activity_ranges))
  
  
return(list(activity_ranges = activity_ranges,
            classification = classification, 
            classification_per_activity_range = classification_per_activity_range,
            freq_acitvity_range = freq_acitvity_range, 
            rel_freq_acitvity_range = rel_freq_acitvity_range,
            quantity_of_bouts = quantity_of_bouts, 
            bout_periods = bout_periods,
            abs_freq_bouts_el = abs_freq_bouts_el))
}
