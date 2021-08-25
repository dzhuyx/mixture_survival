rm(list = ls())

# setwd to root
setwd("/Users/daisyzhu/Documents/Research Projects/SPADE/manuscript/code_submission/mixture_survival")
# setwd
setwd("data_analysis/")
library(rootSolve)
# real data, LHID-2010
# load("Before PSM-for Zheyu.RData")
# synthetic Taiwan data
load("synthetic_Taiwan.rda") # change path as fit
# load relevant functions
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
VARS19 <- unique(c(varlist_long, varlist_short, varlist_shortint))
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

optim_fun <- function(theta) {
      out <- -llh_fun(Y = data$EventDays, d = data$stroke, 
                      X0 = X0, X1 = X1, X2 = X2, 
                      alpha = theta[1:p0], 
                      beta1 = theta[1:p1+p0], 
                      beta2 = theta[1:p2 + p0 + p1],
                      t0 = t0)
      return(out)
}
theta0 <- c(-14, rep(0.5, p0-1), -8, rep(0, p1 - 1), 4, rep(0, p2-1))

set.seed(37)
nstart <- 3000
start <- rbind(mvrnorm(nstart, mu = theta0, Sigma = diag(rep(2, p0+p1+p2))), theta0)

# the following loop can be run as is
# or can be modified to be run as parallel jobs on a computing cluster
# the following implementation is mainly written as demostration
result <- list()
llh <- NULL
for (k in 1:3000) {
      x <- start[k, ]
      fit <- try(optim(x, fn = optim_fun, hessian = T, 
                       control = list(maxit = 500000, reltol = 1e-80, trace = 1)), silent = T)
      result[[k]] <- fit
      llh <- c(llh, fit$value)
}
llh[llh == 999999999] <- NA
fit0 <- result[[which(llh == min(llh, na.rm = T))]]
save(list = c("fit0"), file = "synthetic_MLE_result.rda")

# choose from the Option 1 or 2, run one.

## Option 1: 
## load real data analysis MLE result
load("MLE_result.rda") 
fit0 <- fit
## Option 1 ends

## Option 2:
## load sythetic data analysis MLE result (result from MLestimation.R)
load("synthetic_MLE_result.rda")
## Option 2 ends

est <- fit0$par
se <- sqrt(diag(solve(fit0$hessian)))
pval <- (est + qnorm(0.975) * se) * (est - qnorm(0.975) * se) > 0

# goodness-of-fit tests
gof_mle_poly(data$EventDays, data$stroke,
             X0 = X0, X1 = X1, X2 = X2, fit0 = fit0,
             avec = c(0, quantile(data$EventDays[data$stroke == 1], 1:4/4)),
             t0 = t0)
# all estimation done