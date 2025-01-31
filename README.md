# Analysing Accelerometer Data Using Hidden Markov Models (HMMpa)

`HMMpa`is an R-package providing function to analyses accelerometer data 
(known as a time-series of (impulse)-counts) to quantify length and 
intensity of physical activity using hidden Markov models. 
It also contains the traditional cut-off point method.

 Usually, so called **activity ranges** are used to classify an activity as 
 'sedentary', 'moderate' and so on. Activity ranges are separated 
 by certain thresholds (**cut-off points**). The choice of these cut-off points 
 depends on different components like the subjects' age or the type of accelerometer device. 
 
 Cut-off point values and defined activity ranges are important input values of the 
 following analyzing tools provided by this package:
 
 1. **Cut-off point method** (_assigns an activity range to a count given its 
 total magintude_). This traditional approach assigns an activity range to each count of 
 the time-series independently of each other given its total magnitude.
 
 2. **HMM-based method** (_assigns an activity range to a count given its 
 underlying PA-level_). This approach uses a stochastic model (the hidden Markov model 
 or HMM) to identify the (Markov dependent) time-series of physical activity states 
 underlying the given time-series of accelerometer counts. In contrast to the cut-off 
 point method, this approach assigns activity ranges to the estimated PA-levels 
 corresponding to the hidden activity states, and not directly to the accelerometer counts.
 
 The new procedure for analyzing accelerometer data can be roughly described as follows:
 
 a) First, a hidden Markov model (HMM) is trained to estimate the number `m` of 
 hidden physical activity states and the model specific parameters 
 (`delta, gamma, distribution_theta`). 
 
 b) Then, a user-sepcified decoding algorithm decodes the trained HMM to classify each 
 accelerometer count into the `m` hidden physical activity states.  
 
 c) Finally, the estimated distribution mean values (PA-levels) 
 corresponding to the hidden physical activity states are extracted and the accelerometer 
 counts are assigned by the total magnitudes of their corresponding PA-levels to given 
 physical activity ranges (e.g. 'sedentary', 'light', 'moderate' and 'vigorous') by the 
 traditional cut-off point method.

To see examples, type (after installation)
```R
?HMMpa
```

### Installation 
To install, type in R

```R
install.packages("HMMpa")
# or
devtools::install_github("bips-hb/HMMpa")
```

### Reference
Witowski V, Foraita R, Pitsiladis Y, Pigeot I, Wirsik N (2014). 
_Using Hidden Markov Models to Improve Quantifying Physical Activity in Accelerometer Data – A Simulation Study_
PLOS One. doi:10.1371/journal.pone.0114089