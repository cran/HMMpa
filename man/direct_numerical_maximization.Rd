\name{direct_numerical_maximization}
\alias{direct_numerical_maximization}

\title{
Estimation by Directly Maximizing the log-Likelihood 
}
\description{
Estimates the parameters of a (stationary) discrete-time hidden Markov model by directly maximizing the log-likelihood of the model using the \link[stats]{nlm}-function. See MacDonald & Zucchini (2009, Paragraph 3) for further details.
}
\usage{
direct_numerical_maximization(x, m, delta, gamma, 
     distribution_class, distribution_theta, 
     DNM_limit_accuracy = 0.001, DNM_max_iter = 50, 
     DNM_print = 2)
}

\arguments{
  \item{x}{%%
a vector object containing the time-series of observations that are assumed to be realizations of the (hidden Markov state dependent) observation process of the model.}
  \item{m}{a (finite) number of states in the hidden Markov chain.}
  \item{delta}{%%
a vector object containing starting values for the marginal probability distribution of the \code{m} states of the Markov chain at the time point \code{t=1}. This implementation of the algorithm uses the stationary distribution as delta.}
  \item{gamma}{%%
a matrix (\code{nrow=ncol=m}) containing starting values for the transition matrix of the hidden Markov chain.}
  \item{distribution_class}{%%
a single character string object with the abbreviated name of the \code{m} observation distributions of the Markov dependent observation process. The following distributions are supported by this algorithm: Poisson (\code{pois}); generalized Poisson (\code{genpois}); normal (\code{norm}, discrete log-Likelihood not applicable by this algorithm).}
  \item{distribution_theta}{%%
a list object containing starting values for the parameters of the \code{m} observation distributions that are dependent on the hidden Markov state.}
  \item{DNM_limit_accuracy}{%%
a single numerical value representing the convergence criterion of the direct numerical maximization algorithm using the \link[stats]{nlm}-function. Default value is \code{0.001}.}
  \item{DNM_max_iter}{%%
a single numerical value representing the maximum number of iterations of the direct numerical maximization using the \link[stats]{nlm}-function. Default value is \code{50}.}
  \item{DNM_print}{%%
a single numerical value to determine the level of printing of the \code{nlm}-function.  See \code{nlm}-function for further informations. The value \code{0} suppresses, that no printing will be outputted. Default value is \code{2} for full printing.}
}

\value{
\code{direct_numerical_maximization } returns a list containing the estimated parameters of the hidden Markov model and other components.

  \item{x}{input time-series of observations.}
  \item{m}{input number of hidden states in the Markov chain.}
  \item{logL}{%%
a numerical value representing the logarithmized likelihood calculated by the \code{\link{forward_backward_algorithm}}.}
  \item{AIC}{%%
a numerical value representing Akaike's information criterion for the hidden Markov model with estimated parameters.}
  \item{BIC}{%%
a numerical value representing the Bayesian information criterion for the hidden Markov model with estimated parameters.}
  \item{delta}{%%
a vector object containing the estimates for the marginal probability distribution of the \code{m} states of the Markov chain at time-point point \code{t=1}.}
  \item{gamma}{%%
a matrix containing the estimates for the transition matrix of the hidden Markov chain.}
  \item{distribution_theta}{%%
a list object containing estimates for the parameters of the \code{m} observation distributions that are dependent on the hidden Markov state.}
  \item{distribution_class}{input distribution class.}
}
\references{
MacDonald, I. L., Zucchini, W. (2009) \emph{Hidden Markov Models for Time Series: An Introduction Using R}, Boca Raton: Chapman & Hall.
}
\author{
The basic algorithm of a Poisson-HMM is provided by MacDonald & Zucchini (2009, Paragraph A.1). Extension and implementation by Vitali Witowski (2013).
}

\seealso{
\code{\link{HMM_based_method}}, 
\code{\link{HMM_training}}, 
\code{\link{Baum_Welch_algorithm}}, 
\code{\link{forward_backward_algorithm}},

\code{\link{initial_parameter_training}}
}

\examples{

################################################################
### Fictitious observations ####################################
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


### Assummptions (number of states, probability vector, 
### transition matrix, and distribution parameters)
    m <-4
delta <- c(0.25,0.25,0.25,0.25)
gamma <- 0.7 * diag(m) + rep(0.3 / m)
distribution_class <- "pois"
distribution_theta <- list(lambda = c(4,9,17,25))

### Estimation of a HMM using the method of 
### direct numerical maximization
\donttest{
trained_HMM_with_m_hidden_states <- 
		direct_numerical_maximization(x = x, 
      m = m, 
      delta = delta, 
      gamma = gamma, 
      distribution_class = distribution_class,
      DNM_max_iter=100,
      distribution_theta = distribution_theta)

print(trained_HMM_with_m_hidden_states)
}
}
\keyword{ ts }

