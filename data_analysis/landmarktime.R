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

# without any covariates
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



# general score function for tau
tau_scorefun <- function(tvec, X, harm_fun, phi, par) {
   # tvec: time points to major harm level, a vector
   # X: list of covariate matrices, length = q
   # harm_fun: harm measure function, function of t, X, par
   # phi: major harm level function, function of X and par, returns a vector
   # par: model parameter estimation
   # structure of par should correpond to that of X
   
   # length(tvec) == length(phi output)
   q <- length(tvec)
   major_harm <- phi(X, par) # return q-vector
   harm <- harm_fun(tvec, X, par) # return q-vector
   return(harm - major_harm)
}

# general function for tau estimation and inference
tau_estfun <- function(X, harm_fun, phi, par0, hessian0,
                       tupper = 365, initval = c(7, 30), maxit = 1000) {
   q <- length(initval)
   
   tempfun <- function(t) {
      return(sum(tau_scorefun(t, X, harm_fun, phi, par0)^2))
   }
   fit_tau <- optim(par = initval, fn = tempfun,
                    control = list(maxit = maxit))
   
   tempfun <- function(t) {
      return(tau_scorefun(t, X, harm_fun, phi, par0))
   }
   grad_tau <- jacobian(tempfun, fit_tau$par)
   
   tempfun <- function(par) {
      return(tau_scorefun(fit_tau$par, X, harm_fun, phi, par))
   }
   grad_par <- jacobian(tempfun, par0)
   
   sigma2_tau <- solve(grad_tau) %*% (grad_par) %*% 
      solve(hessian0) %*% t(grad_par) %*% t(solve(grad_tau))
   return(list(est = fit_tau$par,
               sigma = sigma2_tau,
               fit_tau = fit_tau))
}

AH_fun <- function(t, X, par) {
   X1 <- X[[2]]
   X2 <- X[[3]]
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   
   theta1 <- exp(X1 %*% par[1:p1+p0])
   theta2 <- exp(X2 %*% par[1:p2+p0+p1])
   
   out <- theta1 / (1 - theta2) * ((t+t0)^(1-theta2) - t0^(1-theta2)) 
   # t0 = 1, as defined earlier in this script
   return(mean(out))
}

RAH_fun <- function(t, X, par) {
   X0 <- X[[1]]
   X1 <- X[[2]]
   X2 <- X[[3]]
   p0 <- ncol(X0)
   p1 <- ncol(X1)
   p2 <- ncol(X2)
   
   theta0 <- exp(X0 %*% par[1:p0])
   theta1 <- exp(X1 %*% par[1:p1+p0])
   theta2 <- exp(X2 %*% par[1:p2+p0+p1])
   
   baseline <- mean(theta0 * t)
   return(mean(AH_fun(t, X, par)) / baseline)
}

# estimate various time to major harm level in observed population
fit_tau_ah <- tau_estfun(X = list(list(X0, X1, X2), 
                                  list(X0, X1, X2)), 
                         harm_fun = function(tvec, X, par) {
                            return(c(AH_fun(tvec[1], X[[1]], par),
                                     AH_fun(tvec[2], X[[2]], par)))
                         }, 
                         phi = function(X, par) {
                            return(c(0.8 * AH_fun(tupper, X[[1]], par), 
                                     0.5 * AH_fun(tupper, X[[2]], par)))
                         }, 
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian)



fit_tau_rah <- tau_estfun(X = list(list(X0, X1, X2),
                                   list(X0, X1, X2)), 
                          harm_fun = function(tvec, X, par) {
                             return(c(RAH_fun(tvec[1], X[[1]], par),
                                      RAH_fun(tvec[2], X[[2]], par)))
                          }, 
                          phi = function(X, par) {
                             return(c(4, 1))
                          }, 
                          par0 = fit0$par, 
                          hessian0 = fit0$hessian)



# compute measure for subgroups defined according to 
# all 11 binary variables
library(tidyverse)
mycreatecovmat <- function(var, val) {
   result <- list(cbind(1, matrix(unlist(data[data[, var] %in% val, varlist_long]), 
                                  ncol = length(varlist_long))),
                  cbind(1, matrix(unlist(data[data[, var] %in% val, c(varlist_short, varlist_shortint)]), 
                                  ncol = length(varlist_short) + length(varlist_shortint))),
                  cbind(1, matrix(unlist(data[data[, var] %in% val, c(varlist_short, varlist_shortint)]), 
                                  ncol = length(varlist_short) + length(varlist_shortint))))
   return(result)
}
# age groups
X_age <- list(mycreatecovmat("age_young", 1),
              mycreatecovmat("age_middle", 1),
              mycreatecovmat("age_old", 1))
fit_ah_age <- tau_estfun(X = X_age,
                         harm_fun = function(tvec, X, par) {
                            return(c(AH_fun(tvec[1], X[[1]], par),
                                     AH_fun(tvec[2], X[[2]], par),
                                     AH_fun(tvec[3], X[[3]], par),
                                     AH_fun(tvec[4], X[[1]], par),
                                     AH_fun(tvec[5], X[[2]], par),
                                     AH_fun(tvec[6], X[[3]], par)))
                         },
                         phi = function(X, par) {
                            return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                     0.8 * AH_fun(tupper, X[[2]], par),
                                     0.8 * AH_fun(tupper, X[[3]], par),
                                     0.5 * AH_fun(tupper, X[[1]], par),
                                     0.5 * AH_fun(tupper, X[[2]], par),
                                     0.5 * AH_fun(tupper, X[[3]], par)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(30, 30, 30, 7, 7, 7),
                         maxit = 10000)
fit_ah_age$fit_tau$convergence

fit_rah_age <- tau_estfun(X = X_age,
                          harm_fun = function(tvec, X, par) {
                             return(c(RAH_fun(tvec[1], X[[1]], par),
                                      RAH_fun(tvec[2], X[[2]], par),
                                      RAH_fun(tvec[3], X[[3]], par),
                                      RAH_fun(tvec[4], X[[1]], par),
                                      RAH_fun(tvec[5], X[[2]], par),
                                      RAH_fun(tvec[6], X[[3]], par)))
                          },
                          phi = function(X, par) {
                             return(c(rep(4, 3), rep(1, 3)))
                          },
                          par0 = fit0$par, 
                          hessian0 = fit0$hessian, 
                          initval = c(90, 90, 90, 300, 300, 300),
                          maxit = 10000)
fit_rah_age$fit_tau$convergence

# gender groups
X_gender <- list(mycreatecovmat("gender", 0),
                 mycreatecovmat("gender", 1))
fit_ah_gender <- tau_estfun(X = X_gender,
                            harm_fun = function(tvec, X, par) {
                               return(c(AH_fun(tvec[1], X[[1]], par),
                                        AH_fun(tvec[2], X[[2]], par),
                                        AH_fun(tvec[3], X[[1]], par),
                                        AH_fun(tvec[4], X[[2]], par)))
                            },
                            phi = function(X, par) {
                               return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                        0.8 * AH_fun(tupper, X[[2]], par),
                                        0.5 * AH_fun(tupper, X[[1]], par),
                                        0.5 * AH_fun(tupper, X[[2]], par)))
                            },
                            par0 = fit0$par, 
                            hessian0 = fit0$hessian, 
                            initval = c(30, 30, 7, 7),
                            maxit = 10000)
fit_ah_gender$fit_tau$convergence

fit_rah_gender <- tau_estfun(X = X_gender,
                             harm_fun = function(tvec, X, par) {
                                return(c(RAH_fun(tvec[1], X[[1]], par),
                                         RAH_fun(tvec[2], X[[2]], par),
                                         RAH_fun(tvec[3], X[[1]], par),
                                         RAH_fun(tvec[4], X[[2]], par)))
                             },
                             phi = function(X, par) {
                                return(c(rep(4, 2), rep(1, 2)))
                             },
                             par0 = fit0$par, 
                             hessian0 = fit0$hessian, 
                             initval = c(90, 90, 300, 300),
                             maxit = 10000)
fit_rah_gender$fit_tau$convergence

# DM
X_DM <- list(mycreatecovmat("DM", 0),
             mycreatecovmat("DM", 1))
fit_ah_DM <- tau_estfun(X = X_DM,
                        harm_fun = function(tvec, X, par) {
                           return(c(AH_fun(tvec[1], X[[1]], par),
                                    AH_fun(tvec[2], X[[2]], par),
                                    AH_fun(tvec[3], X[[1]], par),
                                    AH_fun(tvec[4], X[[2]], par)))
                        },
                        phi = function(X, par) {
                           return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                    0.8 * AH_fun(tupper, X[[2]], par),
                                    0.5 * AH_fun(tupper, X[[1]], par),
                                    0.5 * AH_fun(tupper, X[[2]], par)))
                        },
                        par0 = fit0$par, 
                        hessian0 = fit0$hessian, 
                        initval = c(30, 30, 7, 7),
                        maxit = 10000)
fit_ah_DM$fit_tau$convergence

fit_rah_DM <- tau_estfun(X = X_DM,
                         harm_fun = function(tvec, X, par) {
                            return(c(RAH_fun(tvec[1], X[[1]], par),
                                     RAH_fun(tvec[2], X[[2]], par),
                                     RAH_fun(tvec[3], X[[1]], par),
                                     RAH_fun(tvec[4], X[[2]], par)))
                         },
                         phi = function(X, par) {
                            return(c(rep(4, 2), rep(1, 2)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(90, 90, 300, 300),
                         maxit = 10000)
fit_rah_DM$fit_tau$convergence

# before_STROKE
X_stroke <- list(mycreatecovmat("before_STROKE", 0),
                 mycreatecovmat("before_STROKE", 1))
fit_ah_stroke <- tau_estfun(X = X_stroke,
                            harm_fun = function(tvec, X, par) {
                               return(c(AH_fun(tvec[1], X[[1]], par),
                                        AH_fun(tvec[2], X[[2]], par),
                                        AH_fun(tvec[3], X[[1]], par),
                                        AH_fun(tvec[4], X[[2]], par)))
                            },
                            phi = function(X, par) {
                               return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                        0.8 * AH_fun(tupper, X[[2]], par),
                                        0.5 * AH_fun(tupper, X[[1]], par),
                                        0.5 * AH_fun(tupper, X[[2]], par)))
                            },
                            par0 = fit0$par, 
                            hessian0 = fit0$hessian, 
                            initval = c(30, 30, 7, 7),
                            maxit = 10000)
fit_ah_stroke$fit_tau$convergence

fit_rah_stroke <- tau_estfun(X = X_stroke,
                             harm_fun = function(tvec, X, par) {
                                return(c(RAH_fun(tvec[1], X[[1]], par),
                                         RAH_fun(tvec[2], X[[2]], par),
                                         RAH_fun(tvec[3], X[[1]], par),
                                         RAH_fun(tvec[4], X[[2]], par)))
                             },
                             phi = function(X, par) {
                                return(c(rep(4, 2), rep(1, 2)))
                             },
                             par0 = fit0$par, 
                             hessian0 = fit0$hessian, 
                             initval = c(90, 90, 300, 300),
                             maxit = 10000)
fit_rah_stroke$fit_tau$convergence

# HT
X_HT <- list(mycreatecovmat("HT", 0),
             mycreatecovmat("HT", 1))
fit_ah_HT <- tau_estfun(X = X_HT,
                        harm_fun = function(tvec, X, par) {
                           return(c(AH_fun(tvec[1], X[[1]], par),
                                    AH_fun(tvec[2], X[[2]], par),
                                    AH_fun(tvec[3], X[[1]], par),
                                    AH_fun(tvec[4], X[[2]], par)))
                        },
                        phi = function(X, par) {
                           return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                    0.8 * AH_fun(tupper, X[[2]], par),
                                    0.5 * AH_fun(tupper, X[[1]], par),
                                    0.5 * AH_fun(tupper, X[[2]], par)))
                        },
                        par0 = fit0$par, 
                        hessian0 = fit0$hessian, 
                        initval = c(30, 30, 7, 7),
                        maxit = 10000)
fit_ah_HT$fit_tau$convergence

fit_rah_HT <- tau_estfun(X = X_HT,
                         harm_fun = function(tvec, X, par) {
                            return(c(RAH_fun(tvec[1], X[[1]], par),
                                     RAH_fun(tvec[2], X[[2]], par),
                                     RAH_fun(tvec[3], X[[1]], par),
                                     RAH_fun(tvec[4], X[[2]], par)))
                         },
                         phi = function(X, par) {
                            return(c(rep(4, 2), rep(1, 2)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(90, 90, 300, 300),
                         maxit = 10000)
fit_rah_HT$fit_tau$convergence

# gs2
X_gs2 <- list(mycreatecovmat("gs2", 0),
              mycreatecovmat("gs2", 1))
fit_ah_gs2 <- tau_estfun(X = X_gs2,
                         harm_fun = function(tvec, X, par) {
                            return(c(AH_fun(tvec[1], X[[1]], par),
                                     AH_fun(tvec[2], X[[2]], par),
                                     AH_fun(tvec[3], X[[1]], par),
                                     AH_fun(tvec[4], X[[2]], par)))
                         },
                         phi = function(X, par) {
                            return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                     0.8 * AH_fun(tupper, X[[2]], par),
                                     0.5 * AH_fun(tupper, X[[1]], par),
                                     0.5 * AH_fun(tupper, X[[2]], par)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(30, 30, 7, 7),
                         maxit = 10000)
fit_ah_gs2$fit_tau$convergence

fit_rah_gs2 <- tau_estfun(X = X_gs2,
                          harm_fun = function(tvec, X, par) {
                             return(c(RAH_fun(tvec[1], X[[1]], par),
                                      RAH_fun(tvec[2], X[[2]], par),
                                      RAH_fun(tvec[3], X[[1]], par),
                                      RAH_fun(tvec[4], X[[2]], par)))
                          },
                          phi = function(X, par) {
                             return(c(rep(4, 2), rep(1, 2)))
                          },
                          par0 = fit0$par, 
                          hessian0 = fit0$hessian, 
                          initval = c(90, 90, 300, 300),
                          maxit = 10000)
fit_rah_gs2$fit_tau$convergence


# AF
X_AF <- list(mycreatecovmat("AF", 0),
             mycreatecovmat("AF", 1))
fit_ah_AF <- tau_estfun(X = X_AF,
                        harm_fun = function(tvec, X, par) {
                           return(c(AH_fun(tvec[1], X[[1]], par),
                                    AH_fun(tvec[2], X[[2]], par),
                                    AH_fun(tvec[3], X[[1]], par),
                                    AH_fun(tvec[4], X[[2]], par)))
                        },
                        phi = function(X, par) {
                           return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                    0.8 * AH_fun(tupper, X[[2]], par),
                                    0.5 * AH_fun(tupper, X[[1]], par),
                                    0.5 * AH_fun(tupper, X[[2]], par)))
                        },
                        par0 = fit0$par, 
                        hessian0 = fit0$hessian, 
                        initval = c(30, 30, 7, 7),
                        maxit = 10000)
fit_ah_AF$fit_tau$convergence

fit_rah_AF <- tau_estfun(X = X_AF,
                         harm_fun = function(tvec, X, par) {
                            return(c(RAH_fun(tvec[1], X[[1]], par),
                                     RAH_fun(tvec[2], X[[2]], par),
                                     RAH_fun(tvec[3], X[[1]], par),
                                     RAH_fun(tvec[4], X[[2]], par)))
                         },
                         phi = function(X, par) {
                            return(c(rep(4, 2), rep(1, 2)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(90, 90, 300, 300),
                         maxit = 10000)
fit_rah_AF$fit_tau$convergence

# dd (hospitalization)
X_dd <- list(mycreatecovmat("dd", 0),
             mycreatecovmat("dd", 1))
fit_ah_dd <- tau_estfun(X = X_dd,
                        harm_fun = function(tvec, X, par) {
                           return(c(AH_fun(tvec[1], X[[1]], par),
                                    AH_fun(tvec[2], X[[2]], par),
                                    AH_fun(tvec[3], X[[1]], par),
                                    AH_fun(tvec[4], X[[2]], par)))
                        },
                        phi = function(X, par) {
                           return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                    0.8 * AH_fun(tupper, X[[2]], par),
                                    0.5 * AH_fun(tupper, X[[1]], par),
                                    0.5 * AH_fun(tupper, X[[2]], par)))
                        },
                        par0 = fit0$par, 
                        hessian0 = fit0$hessian, 
                        initval = c(30, 30, 7, 7),
                        maxit = 10000)
fit_ah_dd$fit_tau$convergence

fit_rah_dd <- tau_estfun(X = X_dd,
                         harm_fun = function(tvec, X, par) {
                            return(c(RAH_fun(tvec[1], X[[1]], par),
                                     RAH_fun(tvec[2], X[[2]], par),
                                     RAH_fun(tvec[3], X[[1]], par),
                                     RAH_fun(tvec[4], X[[2]], par)))
                         },
                         phi = function(X, par) {
                            return(c(rep(4, 2), rep(1, 2)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(90, 90, 300, 300),
                         maxit = 10000)
fit_rah_dd$fit_tau$convergence

# SES
data$SocioeconomicStatus0 <- as.numeric(data$SocioeconomicStatus1 != 1 &
                                           data$SocioeconomicStatus2 != 1)
X_SES <- list(mycreatecovmat("SocioeconomicStatus0", 1),
              mycreatecovmat("SocioeconomicStatus1", 1),
              mycreatecovmat("SocioeconomicStatus2", 1))
fit_ah_SES <- tau_estfun(X = X_SES,
                         harm_fun = function(tvec, X, par) {
                            return(c(AH_fun(tvec[1], X[[1]], par),
                                     AH_fun(tvec[2], X[[2]], par),
                                     AH_fun(tvec[3], X[[3]], par),
                                     AH_fun(tvec[4], X[[1]], par),
                                     AH_fun(tvec[5], X[[2]], par),
                                     AH_fun(tvec[6], X[[3]], par)))
                         },
                         phi = function(X, par) {
                            return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                     0.8 * AH_fun(tupper, X[[2]], par),
                                     0.8 * AH_fun(tupper, X[[3]], par),
                                     0.5 * AH_fun(tupper, X[[1]], par),
                                     0.5 * AH_fun(tupper, X[[2]], par),
                                     0.5 * AH_fun(tupper, X[[3]], par)))
                         },
                         par0 = fit0$par, 
                         hessian0 = fit0$hessian, 
                         initval = c(30, 30, 30, 7, 7, 7),
                         maxit = 10000)
fit_ah_SES$fit_tau$convergence

fit_rah_SES <- tau_estfun(X = X_SES,
                          harm_fun = function(tvec, X, par) {
                             return(c(RAH_fun(tvec[1], X[[1]], par),
                                      RAH_fun(tvec[2], X[[2]], par),
                                      RAH_fun(tvec[3], X[[3]], par),
                                      RAH_fun(tvec[4], X[[1]], par),
                                      RAH_fun(tvec[5], X[[2]], par),
                                      RAH_fun(tvec[6], X[[3]], par)))
                          },
                          phi = function(X, par) {
                             return(c(rep(4, 3), rep(1, 3)))
                          },
                          par0 = fit0$par, 
                          hessian0 = fit0$hessian, 
                          initval = c(90, 90, 90, 300, 300, 300),
                          maxit = 10000)
fit_rah_SES$fit_tau$convergence

# area
data$area16 <- as.numeric(data$area2 != 1 &
                             data$area3 != 1 &
                             data$area4 != 1 & 
                             data$area5 != 1)
X_area <- list(mycreatecovmat("area2", 1),
               mycreatecovmat("area3", 1),
               mycreatecovmat("area4", 1),
               mycreatecovmat("area5", 1),
               mycreatecovmat("area16", 1))
fit_ah_area <- tau_estfun(X = X_area,
                          harm_fun = function(tvec, X, par) {
                             return(c(AH_fun(tvec[1], X[[1]], par),
                                      AH_fun(tvec[2], X[[2]], par),
                                      AH_fun(tvec[3], X[[3]], par),
                                      AH_fun(tvec[4], X[[4]], par),
                                      AH_fun(tvec[5], X[[5]], par),
                                      AH_fun(tvec[6], X[[1]], par),
                                      AH_fun(tvec[7], X[[2]], par),
                                      AH_fun(tvec[8], X[[3]], par),
                                      AH_fun(tvec[9], X[[4]], par),
                                      AH_fun(tvec[10], X[[5]], par)))
                          },
                          phi = function(X, par) {
                             return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                      0.8 * AH_fun(tupper, X[[2]], par),
                                      0.8 * AH_fun(tupper, X[[3]], par),
                                      0.8 * AH_fun(tupper, X[[4]], par),
                                      0.8 * AH_fun(tupper, X[[5]], par),
                                      0.5 * AH_fun(tupper, X[[1]], par),
                                      0.5 * AH_fun(tupper, X[[2]], par),
                                      0.5 * AH_fun(tupper, X[[3]], par),
                                      0.5 * AH_fun(tupper, X[[4]], par),
                                      0.5 * AH_fun(tupper, X[[5]], par)))
                          },
                          par0 = fit0$par, 
                          hessian0 = fit0$hessian, 
                          initval = c(rep(30, 5), rep(7, 5)),
                          maxit = 10000)
fit_ah_area$fit_tau$convergence

fit_rah_area <- tau_estfun(X = X_area,
                           harm_fun = function(tvec, X, par) {
                              return(c(RAH_fun(tvec[1], X[[1]], par),
                                       RAH_fun(tvec[2], X[[2]], par),
                                       RAH_fun(tvec[3], X[[3]], par),
                                       RAH_fun(tvec[4], X[[4]], par), 
                                       RAH_fun(tvec[5], X[[5]], par),
                                       RAH_fun(tvec[6], X[[1]], par),
                                       RAH_fun(tvec[7], X[[2]], par),
                                       RAH_fun(tvec[8], X[[3]], par),
                                       RAH_fun(tvec[9], X[[4]], par), 
                                       RAH_fun(tvec[10], X[[5]], par)))
                           },
                           phi = function(X, par) {
                              return(c(rep(4, 5), rep(1, 5)))
                           },
                           par0 = fit0$par, 
                           hessian0 = fit0$hessian, 
                           initval = c(rep(90, 5), rep(300, 5)),
                           maxit = 10000)
fit_rah_area$fit_tau$convergence

# Teaching
data$Teaching1 <- as.numeric(data$Teaching0 != 1 &
                                data$Teaching2 != 1 &
                                data$Teaching3 != 1 &
                                data$Teaching4 != 1)
X_teaching <- list(mycreatecovmat("Teaching0", 1),
                   mycreatecovmat("Teaching2", 1),
                   mycreatecovmat("Teaching3", 1),
                   mycreatecovmat("Teaching4", 1),
                   mycreatecovmat("Teaching1", 1))
fit_ah_teaching <- tau_estfun(X = X_teaching,
                              harm_fun = function(tvec, X, par) {
                                 return(c(AH_fun(tvec[1], X[[1]], par),
                                          AH_fun(tvec[2], X[[2]], par),
                                          AH_fun(tvec[3], X[[3]], par),
                                          AH_fun(tvec[4], X[[4]], par),
                                          AH_fun(tvec[5], X[[5]], par),
                                          AH_fun(tvec[6], X[[1]], par),
                                          AH_fun(tvec[7], X[[2]], par),
                                          AH_fun(tvec[8], X[[3]], par),
                                          AH_fun(tvec[9], X[[4]], par),
                                          AH_fun(tvec[10], X[[5]], par)))
                              },
                              phi = function(X, par) {
                                 return(c(0.8 * AH_fun(tupper, X[[1]], par),
                                          0.8 * AH_fun(tupper, X[[2]], par),
                                          0.8 * AH_fun(tupper, X[[3]], par),
                                          0.8 * AH_fun(tupper, X[[4]], par),
                                          0.8 * AH_fun(tupper, X[[5]], par),
                                          0.5 * AH_fun(tupper, X[[1]], par),
                                          0.5 * AH_fun(tupper, X[[2]], par),
                                          0.5 * AH_fun(tupper, X[[3]], par),
                                          0.5 * AH_fun(tupper, X[[4]], par),
                                          0.5 * AH_fun(tupper, X[[5]], par)))
                              },
                              par0 = fit0$par, 
                              hessian0 = fit0$hessian, 
                              initval = c(rep(30, 5), rep(7, 5)),
                              maxit = 10000)
fit_ah_teaching$fit_tau$convergence

fit_rah_teaching <- tau_estfun(X = X_teaching,
                               harm_fun = function(tvec, X, par) {
                                  return(c(RAH_fun(tvec[1], X[[1]], par),
                                           RAH_fun(tvec[2], X[[2]], par),
                                           RAH_fun(tvec[3], X[[3]], par),
                                           RAH_fun(tvec[4], X[[4]], par), 
                                           RAH_fun(tvec[5], X[[5]], par),
                                           RAH_fun(tvec[6], X[[1]], par),
                                           RAH_fun(tvec[7], X[[2]], par),
                                           RAH_fun(tvec[8], X[[3]], par),
                                           RAH_fun(tvec[9], X[[4]], par), 
                                           RAH_fun(tvec[10], X[[5]], par)))
                               },
                               phi = function(X, par) {
                                  return(c(rep(4, 5), rep(1, 5)))
                               },
                               par0 = fit0$par, 
                               hessian0 = fit0$hessian, 
                               initval = c(rep(90, 5), rep(300, 5)),
                               maxit = 10000)
fit_rah_teaching$fit_tau$convergence



save.image(file = "synthetic_landmarktime.rda")

