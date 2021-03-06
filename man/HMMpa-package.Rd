\name{HMMpa-package}
\alias{HMMpa-package}
\alias{HMMpa}
\docType{package}
\title{
Analysing Accelerometer Data Using Hidden Markov Models
}
\description{
This package provides functions for analyzing accelerometer outpout data (known as a time-series of (impulse)-counts) to quantify length and intensity of physical activity.

Usually, so called \emph{activity ranges} are used to classify an activity as \dQuote{sedentary}, \dQuote{moderate} and so on. \emph{Activity ranges} are separated by certain thresholds (\emph{cut-off points}). The choice of these cut-off points depends on different components like the subjects' age or the type of accelerometer device. 

Cut-off point values and defined activity ranges are important input values of the following analyzing tools provided by this package:

1. \strong{Cut-off point method} (\emph{assigns an activity range to a count given its total magintude}).
This traditional approach assigns an activity range to each count of the time-series independently of each other given its total magnitude.

2. \strong{HMM-based method} (\emph{assigns an activity range to a count given its underlying PA-level}).
This approach uses a stochastic model (the hidden Markov model or HMM) to identify the (Markov dependent) time-series of physical activity states underlying the given time-series of accelerometer counts. In contrast to the cut-off point method, this approach assigns activity ranges to the estimated PA-levels corresponding to the hidden activity states, and not directly to the accelerometer counts.


}
\details{
\tabular{ll}{
Package: \tab HMMpa\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2018-04-20\cr
License: \tab GPL-3 \cr
}

The new procedure for analyzing accelerometer data can be roughly described as follows:\cr
First, a hidden Markov model (HMM) is trained to estimate the number \code{m} of hidden physical activity states and the model specific parameters (\code{delta, gamma, distribution_theta}).  Then, a user-sepcified decoding algorithm decodes the trainded HMM to classify each accelerometer count into the \code{m} hidden physical activity states.  Finally, the estimated distribution mean values (PA-levels) corresponding to the hidden physical activity states are extracted and the accelerometer counts are assigned by the total magnitudes of their corresponding PA-levels to given physical activity ranges (e.g. "sedentary", "light", "moderate" and "vigorous") by the traditional cut-off point method.


}
\author{
Vitali Witowski, \cr
Ronja Foraita,
Leibniz Institute for Prevention Research and Epidemiology (BIPS) \cr

Maintainer: Ronja Foraita <foraita@leibniz-bips.de>
}
\references{
Baum, L., Petrie, T., Soules, G., Weiss, N. (1970). A maximization technique occurring in the statistical analysis of probabilistic functions of markov chains. The annals of mathematical statistics, vol. \bold{41}(1), 164--171.

Brachmann, B. (2011). Hidden-Markov-Modelle fuer Akzelerometerdaten. Diploma Thesis, University Bremen - Bremen Institute for Prevention Research and Social Medicine (BIPS).

Dempster, A., Laird, N., Rubin, D. (1977). Maximum likelihood from incomplete data via the EM algorithm. Journal of the Royal Statistical Society. Series B (Methodological), vol. \bold{39}(1), 1--38.

Forney, G.D. (1973). The Viterbi algorithm. Proceeding of the IEE, vol. \bold{61}(3), 268--278.

Joe, H., Zhu, R. (2005). Generalized poisson distribution: the property of mixture of poisson and comparison with negative binomial distribution. Biometrical Journal, vol. \bold{47}(2), 219--229. 

MacDonald, I. L., Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.

Viterbi, A.J. (1967). Error Bounds for concolutional codes and an asymptotically optimal decoding algorithm. Information Theory, IEEE Transactions on, vol. \bold{13}(2), 260--269.

Witowski, V., Foraita, R., Pitsiladis, Y., Pigeot, I., Wirsik, N. (2014) Using hidden Markov models to improve quantifying physical activity in accelerometer data - A simulation study. PLOS ONE. \bold{9}(12), e114089. 
http://dx.doi.org/10.1371/journal.pone.0114089


}
\note{
  We thank Moritz Hanke for his help in realizing this package.
}

\keyword{ ts }
\keyword{ iteration }

\examples{

################################################################
###     Example 1  (traditional approach)    ################### 
###     Solution of the cut-off point method ################### 
################################################################

################################################################
### Fictitious activity counts #################################
################################################################

x <- c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,
  14,16,13,11,11,10,12,19,23,25,24,23,20,21,22,22,18,7,
  5,3,4,3,2,3,4,5,4,2,1,3,4,5,4,5,3,5,6,4,3,6,4,8,9,12,
  9,14,17,15,25,23,25,35,29,36,34,36,29,41,42,39,40,43,
  37,36,20,20,21,22,23,26,27,28,25,28,24,21,25,21,20,21,
  11,18,19,20,21,13,19,18,20,7,18,8,15,17,16,13,10,4,9,
  7,8,10,9,11,9,11,10,12,12,5,13,4,6,6,13,8,9,10,13,13,
  11,10,5,3,3,4,9,6,8,3,5,3,2,2,1,3,5,11,2,3,5,6,9,8,5,
  2,5,3,4,6,4,8,15,12,16,20,18,23,18,19,24,23,24,21,26,
  36,38,37,39,45,42,41,37,38,38,35,37,35,31,32,30,20,39,
  40,33,32,35,34,36,34,32,33,27,28,25,22,17,18,16,10,9,
  5,12,7,8,8,9,19,21,24,20,23,19,17,18,17,22,11,12,3,9,
  10,4,5,13,3,5,6,3,5,4,2,5,1,2,4,4,3,2,1)  
\donttest{    
traditional_cut_off_point_method <- cut_off_point_method(x=x, 
 	cut_points = c(5,15,23), 
 	names_activity_ranges = c("SED","LIG","MOD","VIG"), 
 	bout_lengths = c(1,1,2,4,5,10,11,20,21,60,61,260),
 	plotting = 1)
}


################################################################ 
################################################################ 
###      Examples 2,3 and 4  (new approach)             ########
###      Solution of the HMM based cut-off point method ######## 
################################################################
###      Demonstrated both in three steps (Example 2)     ###### 
###      and condensed in one function (Examples 3 and 4) ######
################################################################

################################################################
### Example 2) Manually in three steps    ######################
################################################################

################################################################
### Fictitious activity counts #################################
################################################################
                                                                
x <- c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,           
  14,16,13,11,11,10,12,19,23,25,24,23,20,21,22,22,18,7,         
  5,3,4,3,2,3,4,5,4,2,1,3,4,5,4,5,3,5,6,4,3,6,4,8,9,12,         
  9,14,17,15,25,23,25,35,29,36,34,36,29,41,42,39,40,43,         
  37,36,20,20,21,22,23,26,27,28,25,28,24,21,25,21,20,21,        
  11,18,19,20,21,13,19,18,20,7,18,8,15,17,16,13,10,4,9,         
  7,8,10,9,11,9,11,10,12,12,5,13,4,6,6,13,8,9,10,13,13,         
  11,10,5,3,3,4,9,6,8,3,5,3,2,2,1,3,5,11,2,3,5,6,9,8,5,         
  2,5,3,4,6,4,8,15,12,16,20,18,23,18,19,24,23,24,21,26,         
  36,38,37,39,45,42,41,37,38,38,35,37,35,31,32,30,20,39,        
  40,33,32,35,34,36,34,32,33,27,28,25,22,17,18,16,10,9,         
  5,12,7,8,8,9,19,21,24,20,23,19,17,18,17,22,11,12,3,9,         
  10,4,5,13,3,5,6,3,5,4,2,5,1,2,4,4,3,2,1)                      
     
################################################################
## Step 1: Training of a HMM ##################################
##         for the given time-series of counts ################ 
##         ####################################################
##         More precisely: training of a poisson ##############
##         distribution based hidden Markov model for #########
##         number of states m=2,...,6 #########################
##         and selection of the model with the most ########### 
##         plausible m ########################################
################################################################
\donttest{
m_trained_HMM <- HMM_training(x = x, 
 min_m = 2, 
 max_m = 6, 
 distribution_class = "pois")$trained_HMM_with_selected_m  
 	 
###############################################################
## Step 2: Decoding of the trained HMM for the given ##########
##         time-series of accelerometer counts to extract #####
##         hidden PA-levels ###################################
###############################################################
hidden_PA_levels <- HMM_decoding(x = x, 
 	 m = m_trained_HMM$m, 
 	 delta = m_trained_HMM$delta, 
 	 gamma = m_trained_HMM$gamma, 
 	 distribution_class = m_trained_HMM
 	 $distribution_class, 
 	 distribution_theta = m_trained_HMM$
 	 distribution_theta)$decoding_distr_means

############################################################### 
## Step 3: Assigning of user-sepcified activity ranges ########
##         to the accelerometer counts via the total ##########
##         magnitudes of their corresponding ##################
##         hidden PA-levels ###################################
##         ####################################################
##         In this example four activity ranges ###############
##         (named as "sedentary", "light", "moderate" #########
##         and "vigorous" physical activity) are ##############
##         separated by the three cut-points 5, 15 and 23) ####
################################################################
HMM_based_cut_off_point_method <- cut_off_point_method(x = x, 
 	 hidden_PA_levels = hidden_PA_levels, 
 	 cut_points = c(5,15,23), 
 	 names_activity_ranges = c("SED","LIG","MOD","VIG"), 
 	 bout_lengths = c(1,1,2,4,5,10,11,20,21,60,61,260),
 	 plotting=1)      				 
}



###############################################################
## Example 3) In a single function (Step 1-3 of Example 2  ####
##            combined in one function)                    ####
###############################################################
	
###############################################################
## Fictitious activity counts #################################
###############################################################
                                                                
x <- c(1,16,19,34,22,6,3,5,6,3,4,1,4,3,5,7,9,8,11,11,           
  14,16,13,11,11,10,12,19,23,25,24,23,20,21,22,22,18,7,         
  5,3,4,3,2,3,4,5,4,2,1,3,4,5,4,5,3,5,6,4,3,6,4,8,9,12,         
  9,14,17,15,25,23,25,35,29,36,34,36,29,41,42,39,40,43,         
  37,36,20,20,21,22,23,26,27,28,25,28,24,21,25,21,20,21,        
  11,18,19,20,21,13,19,18,20,7,18,8,15,17,16,13,10,4,9,         
  7,8,10,9,11,9,11,10,12,12,5,13,4,6,6,13,8,9,10,13,13,         
  11,10,5,3,3,4,9,6,8,3,5,3,2,2,1,3,5,11,2,3,5,6,9,8,5,         
  2,5,3,4,6,4,8,15,12,16,20,18,23,18,19,24,23,24,21,26,         
  36,38,37,39,45,42,41,37,38,38,35,37,35,31,32,30,20,39,        
  40,33,32,35,34,36,34,32,33,27,28,25,22,17,18,16,10,9,         
  5,12,7,8,8,9,19,21,24,20,23,19,17,18,17,22,11,12,3,9,         
  10,4,5,13,3,5,6,3,5,4,2,5,1,2,4,4,3,2,1)                      

###############################################################
## Use a (m=4 state) hidden Markov model based on the #########
## generalized poisson distribution to assign an      #########
## activity range to the counts                       #########
###############################################################
## In this example three activity ranges                   ####
## (named as "light", "moderate" and "vigorous" physical   ####
## activity) are separated by the two cut-points 15 and 23 ####
###############################################################
\donttest{
HMM_based_cut_off_point_method <- HMM_based_method(x = x, 
 	 cut_points = c(15,23), 
 	 min_m = 4, 
 	 max_m = 4, 
 	 names_activity_ranges = c("LIG","MOD","VIG"), 
 	 distribution_class = "genpois", 
 	 training_method = "numerical",
 	 DNM_limit_accuracy = 0.05,
 	 DNM_max_iter = 10,
 	 bout_lengths = c(1,1,2,4,5,10,11,20,21,60,61,260),
 	 plotting = 1)
} 	 

###############################################################
## Example 4) In a single function (Step 1-3 of Example 2  ####
##            combined in one function)                    ####
##            (large and highly scatterd time-series)      ####
###############################################################
	 
################################################################
### Generate a large time-series of highly scattered counts ####
################################################################

x <- HMM_simulation(
 size = 1500, 
 m = 10,
 gamma = 0.93 * diag(10) + rep(0.07 / 10, times = 10),
 distribution_class = "norm", 
 distribution_theta = list(mean = c(10, 100, 200, 300, 450, 
 600, 700, 900, 1100, 1300, 1500), 
 sd=c(rep(100,times=10))), 
 obs_round=TRUE, 
 obs_non_neg=TRUE,
 plotting=5)$observations

################################################################
### Compare results of the tradional cut-point method ##########
### and the (6-state-normal-)HMM based method ##################
################################################################
\donttest{		
traditional_cut_off_point_method <- cut_off_point_method(x=x, 
 	 cut_points = c(200,500,1000), 
 	 names_activity_ranges = c("SED","LIG","MOD","VIG"), 
 	 bout_lengths = c(1,1,2,4,5,10,11,20,21,60,61,260),
 	 plotting = 1)

HMM_based_cut_off_point_method <- HMM_based_method(x=x,
	 max_scaled_x = 200, 
 	 cut_points = c(200,500,1000), 
 	 min_m = 6, 
 	 max_m = 6,
 	 BW_limit_accuracy = 0.5, 
 	 BW_max_iter = 10,
 	 names_activity_ranges = c("SED","LIG","MOD","VIG"), 
 	 distribution_class = "norm", 
 	 bout_lengths = c(1,1,2,4,5,10,11,20,21,60,61,260),
 	 plotting = 1)
}
}
