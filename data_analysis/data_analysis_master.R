# flow of the complete data analyses
# necessary packages
require(numDeriv)
require(MASS)
require(rootSolve)
require(Rcpp)
require(tidyverse)
require(ggplot2)
require(scales)
require(patchwork)

# MLE
souce("MLestimation.R") # output MLE estimates using synthetic Taiwan data

# harm risk factor analysis
# input: rda file ouput by MLestimation.R or MLE_result.rda; 
#		 and synthetic_Taiwan.rda
# output: synthetic_harmriskfactor.rda
source("harmriskfactor.R")

# landmark time analysis
# input: rda file ouput by MLestimation.R or MLE_result.rda; 
#		 and synthetic_Taiwan.rda
# output: synthetic_landmarktime.rda
source("landmarktime.R")

# profile analysis
# input: rda file ouput by MLestimation.R or MLE_result.rda; 
#		 and synthetic_Taiwan.rda
# output: synthetic_gsprofile.rda
source("profile_analysis.R")

# create figures summarize results
# input: rda file ouput by MLestimation.R or MLE_result.rda,
# 		 synthetic_Taiwan.rda, synthetic_harmriskfactor.rda,
# 		 synthetic_landmarktime.rda, and synthetic_gsprofile.rda
# output: figures 2, 3, 4, and table 3
source("illustration.R")
