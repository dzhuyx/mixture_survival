rm(list = ls())
library(rootSolve)
# real data, LHID-2010
# load("Before PSM-for Zheyu.RData")
load("synthetic_Taiwan.rda") # load synthetic data
source("parmix_function_poly.R")
# load variable list 
varlist_long <- c("gs2", "gender", "HT", "DM", "AF", "before_STROKE", "dd",                  
                  "SocioeconomicStatus2", "SocioeconomicStatus1", 
                  "area4", "area2", "area3", "area5",
                  "Teaching2", "age_young", "age_old")
varlist_short <- c("gs2", "gender", "HT", "DM", "before_STROKE",
                   "area4", "Teaching4", "Teaching3",
                   "Teaching0", "age_young", "age_old")
varlist_shortint <- c("gs2_age_old")

t0 <- 1
tupper <- 365
VARS19 <- unique(c(varlist_long, varlist_short, varlist_shortint, "age_middle"))

data$EventDays1 <- data$EventDays
data$stroke1 <- 1
data$stroke1[data$EventDays > 365] <- 0
data$EventDays1[data$EventDays > 365] <- 365

data$EventDays <- data$EventDays1
data$stroke <- data$stroke1
data <- data[, c("ID", "stroke", "EventDays", VARS19)]

gs1 <- data$gs2 == 0
gs2 <- data$gs2 == 1

n1 <- sum(gs1)
n2 <- sum(gs2)
n <- nrow(data)

X0 <- cbind(1, matrix(unlist(data[, varlist_long]), nrow = n))
X1 <- cbind(1, matrix(unlist(data[, varlist_short]), nrow = n),
            matrix(unlist(data[, varlist_shortint]), nrow = n))
X2 <- X1
p0 <- ncol(X0)
p1 <- ncol(X1)
p2 <- ncol(X2)

# choose from the Option 1 or 2, run one. Default is Option 1

## Option 1: 
## load real data analysis MLE result
load("MLE_result.rda") 
fit0 <- fit
## Option 1 ends

# ## Option 2:
# ## load sythetic data analysis MLE result (result obtained from running MLestimation.R)
# load("synthetic_MLE_result.rda")
# ## Option 2 ends

AH_fun <- function(theta1, theta2, t) {
   out <- theta1 / (1 - theta2) * ((t+t0)^(1-theta2) - t0^(1-theta2))
   return(out)
}
AAH_fun <- function(theta0, theta1, theta2, t) {
   out <- theta1 / (1 - theta2) * ((t+t0)^(1-theta2) - t0^(1-theta2)) / 
      (theta0 * t)
   return(out)
}
AAHS_fun <- function(theta0, theta1, theta2, t) {
   out <- theta1 / (1 - theta2) * ((t+t0)^(1-theta2) - t0^(1-theta2)) / 
      (1 - theta0 * t)
   return(out)
}

AH_observed <- function(par, X1, X2, t) {
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   
   theta1 <- exp(X1 %*% par[1:p1+p0])
   theta2 <- exp(X2 %*% par[1:p2+p0+p1])
   return(mean(AH_fun(theta1, theta2, t)))
}

AAH_observed <- function(par, X0, X1, X2, t) {
   p0 <- ncol(X0)
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   theta0 <- exp(X0 %*% par[1:p0])
   theta1 <- exp(X1 %*% par[1:p1+p0])
   theta2 <- exp(X2 %*% par[1:p2+p0+p1])
   
   baseline <- mean(theta0 * t)
   return(mean(AH_fun(theta1, theta2, t)) / baseline)
}

AAHS_observed <- function(par, X0, X1, X2, t) {
   p0 <- ncol(X0)
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   theta0 <- exp(X0 %*% par[1:p0])
   theta1 <- exp(X1 %*% par[1:p1+p0])
   theta2 <- exp(X2 %*% par[1:p2+p0+p1])
   
   baseline <- mean(theta0 * t)
   return(mean(AH_fun(theta1, theta2, t)) / (1 - baseline))
}

# profile analysis: observed
ah_observed_gs1 <- rep(NA, tupper)
ah_logse_observed_gs1 <- rep(NA, tupper)
ah_observed_gs2 <- rep(NA, tupper)
ah_logse_observed_gs2 <- rep(NA, tupper)

aah_observed_gs1 <- rep(NA, tupper)
aah_logse_observed_gs1 <- rep(NA, tupper)
aah_observed_gs2 <- rep(NA, tupper)
aah_logse_observed_gs2 <- rep(NA, tupper)

aahs_observed_gs1 <- rep(NA, tupper)
aahs_logse_observed_gs1 <- rep(NA, tupper)
aahs_observed_gs2 <- rep(NA, tupper)
aahs_logse_observed_gs2 <- rep(NA, tupper)

# profile analysis: standard
ah_standard_gs1 <- rep(NA, tupper)
ah_logse_standard_gs1 <- rep(NA, tupper)
ah_standard_gs2 <- rep(NA, tupper)
ah_logse_standard_gs2 <- rep(NA, tupper)

aah_standard_gs1 <- rep(NA, tupper)
aah_logse_standard_gs1 <- rep(NA, tupper)
aah_standard_gs2 <- rep(NA, tupper)
aah_logse_standard_gs2 <- rep(NA, tupper)

aahs_standard_gs1 <- rep(NA, tupper)
aahs_logse_standard_gs1 <- rep(NA, tupper)
aahs_standard_gs2 <- rep(NA, tupper)
aahs_logse_standard_gs2 <- rep(NA, tupper)

# profile analysis: adjusted
ah_adjusted_gs1 <- rep(NA, tupper)
ah_logse_adjusted_gs1 <- rep(NA, tupper)
ah_adjusted_gs2 <- rep(NA, tupper)
ah_logse_adjusted_gs2 <- rep(NA, tupper)

aah_adjusted_gs1 <- rep(NA, tupper)
aah_logse_adjusted_gs1 <- rep(NA, tupper)
aah_adjusted_gs2 <- rep(NA, tupper)
aah_logse_adjusted_gs2 <- rep(NA, tupper)

aahs_adjusted_gs1 <- rep(NA, tupper)
aahs_logse_adjusted_gs1 <- rep(NA, tupper)
aahs_adjusted_gs2 <- rep(NA, tupper)
aahs_logse_adjusted_gs2 <- rep(NA, tupper)

for (t in 1:tupper) {
   X0 <- cbind(1, matrix(unlist(data[, varlist_long]), nrow = n))
   X1 <- cbind(1, matrix(unlist(data[, varlist_short]), nrow = n),
               matrix(unlist(data[, varlist_shortint]), nrow = n))
   X2 <- cbind(1, gs2)
   
   X0gs1 <- cbind(1, 0, matrix(unlist(data[, varlist_long[-1]]), nrow = n))
   X0gs2 <- cbind(1, 1, matrix(unlist(data[, varlist_long[-1]]), nrow = n))
   
   X1gs1 <- cbind(1, 0, matrix(unlist(data[, varlist_short[-1]]), nrow = n),
                  0 * matrix(unlist(data[, varlist_shortint]), nrow = n))
   X1gs2 <- cbind(1, 1, matrix(unlist(data[, varlist_short[-1]]), nrow = n),
                  1 * matrix(unlist(data[, varlist_shortint]), nrow = n))
   
   X2gs1 <- cbind(1, rep(0, n))
   X2gs2 <- cbind(1, rep(1, n))
   
   print(t)
   # observed performance based on AH, AAH and AAHS
   # observe AH
   ah_observed_gs1[t] <- AH_observed(fit0$par, X1[gs1, ], X2[gs1, , drop = F], t)
   tempfun <- function(par) {
      return(AH_observed(par, X1[gs1, ], X2[gs1, , drop = F], t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_observed_gs1[t] <- sqrt(gradtemp %*% 
                                       solve(fit0$hessian) %*% 
                                       t(gradtemp) / ah_observed_gs1[t]^2)
   
   ah_observed_gs2[t] <- AH_observed(fit0$par, X1[gs2, ], X2[gs2, , drop = F], t)
   tempfun <- function(par) {
      return(AH_observed(par, X1[gs2, ], X2[gs2, , drop = F], t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_observed_gs2[t] <- sqrt(gradtemp %*% 
                                       solve(fit0$hessian) %*% 
                                       t(gradtemp) / ah_observed_gs2[t]^2)
   
   
   
   # observed AAH
   aah_observed_gs1[t] <- AAH_observed(fit0$par, X0[gs1, ], X1[gs1, ], 
                                       X2[gs1, , drop = F], t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0[gs1, ], X1[gs1, ], 
                          X2[gs1,  , drop = F], t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_observed_gs1[t] <- sqrt(gradtemp %*% 
                                        solve(fit0$hessian) %*% 
                                        t(gradtemp) / aah_observed_gs1[t]^2)
   
   aah_observed_gs2[t] <- AAH_observed(fit0$par, X0[gs2, ], X1[gs2, ], 
                                       X2[gs2, , drop = F], t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0[gs2, ], X1[gs2, ], 
                          X2[gs2, , drop = F], t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_observed_gs2[t] <- sqrt(gradtemp %*% 
                                        solve(fit0$hessian) %*% 
                                        t(gradtemp) / aah_observed_gs2[t]^2)
   
   
   # standard performance based on AAH and AAHS
   # standard AH
   ah_standard_gs1[t] <- (AH_observed(fit0$par, X1 = X1[gs1, ],
                                      X2 = X2[gs1, , drop = F], t) + 
                             AH_observed(fit0$par, X1 = X1gs2[gs1, ], 
                                         X2 = X2gs2[gs1, , drop = F], t)) / 2
   tempfun <- function(par) {
      return((AH_observed(par, X1 = X1[gs1, ], X2 = X2[gs1, , drop = F], t) + 
                 AH_observed(par, X1 = X1gs2[gs1, ], 
                             X2 = X2gs2[gs1, , drop = F], t)) / 2)
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_standard_gs1[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                    %*% t(gradtemp) / ah_standard_gs1[t]^2)
   
   
   ah_standard_gs2[t] <- (AH_observed(fit0$par, X1 = X1[gs2, ], 
                                      X2 = X2[gs2, , drop = F], t) + 
                             AH_observed(fit0$par, X1 = X1gs1[gs2, ], 
                                         X2 = X2gs1[gs2, , drop = F], t)) / 2
   tempfun <- function(par) {
      return((AH_observed(par, X1 = X1[gs2, ], X2 = X2[gs2, , drop = F], t) + 
                 AH_observed(par, X1 = X1gs1[gs2, ], 
                             X2 = X2gs1[gs2, , drop = F], t)) / 2)
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_standard_gs2[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                    %*% t(gradtemp) / ah_standard_gs2[t]^2)
   
   
   # standard AAH
   aah_standard_gs1[t] <- AAH_observed(fit0$par, X0 = rbind(X0[gs1, ], X0[gs1, ]), 
                                       X1 = rbind(X1[gs1, ], X1gs2[gs1, ]), 
                                       X2 = rbind(X2[gs1, , drop = F], X2gs2[gs1, , drop = F]), t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0 = rbind(X0[gs1, ], X0[gs1, ]), 
                          X1 = rbind(X1[gs1, ], X1gs2[gs1, ]), 
                          X2 = rbind(X2[gs1, , drop = F], X2gs2[gs1, , drop = F]), t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_standard_gs1[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                     %*% t(gradtemp) / aah_standard_gs1[t]^2)
   
   
   aah_standard_gs2[t] <- AAH_observed(fit0$par, X0 = rbind(X0[gs2, ], X0[gs2, ]), 
                                       X1 = rbind(X1[gs2, ], X1gs1[gs2, ]), 
                                       X2 = rbind(X2[gs2, , drop = F], X2gs1[gs2, , drop = F]), t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0 = rbind(X0[gs2, ], X0[gs2, ]), 
                          X1 = rbind(X1[gs2, ], X1gs1[gs2, ]), 
                          X2 = rbind(X2[gs2, , drop = F], X2gs1[gs2, , drop = F]), t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_standard_gs2[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                     %*% t(gradtemp) / aah_standard_gs2[t]^2)
   
   
   
   # adjusted performance
   # adjusted AH
   ah_adjusted_gs1[t] <- AH_observed(fit0$par, X1 = X1gs1, X2 = X2gs1, t)
   tempfun <- function(par) {
      return(AH_observed(par, X1 = X1gs1, X2 = X2gs1, t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_adjusted_gs1[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                    %*% t(gradtemp) / ah_adjusted_gs1[t]^2)
   
   ah_adjusted_gs2[t] <- AH_observed(fit0$par, X1 = X1gs2, X2 = X2gs2, t)
   tempfun <- function(par) {
      return(AH_observed(par, X1 = X1gs2, X2 = X2gs2, t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   ah_logse_adjusted_gs2[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                    %*% t(gradtemp) / ah_adjusted_gs2[t]^2)
   
   # adjusted AAH
   aah_adjusted_gs1[t] <- AAH_observed(fit0$par, X0 = X0, X1 = X1gs1,
                                       X2 = X2gs1, t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0 = X0, X1 = X1gs1, X2 = X2gs1, t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_adjusted_gs1[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                     %*% t(gradtemp) / aah_adjusted_gs1[t]^2)
   
   aah_adjusted_gs2[t] <- AAH_observed(fit0$par, X0 = X0, X1 = X1gs2,
                                       X2 = X2gs2, t)
   tempfun <- function(par) {
      return(AAH_observed(par, X0 = X0, X1 = X1gs2,
                          X2 = X2gs2, t))
   }
   gradtemp <- gradient(tempfun, fit0$par)
   aah_logse_adjusted_gs2[t] <- sqrt(gradtemp %*% solve(fit0$hessian) 
                                     %*% t(gradtemp) / aah_adjusted_gs2[t]^2)
   
}
save.image(file = "synthetic_gsprofile.rda")

