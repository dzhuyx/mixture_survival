# create synthetic Taiwan data to illustrate data analyses
# data is generated using summary statistics and estimated coefficients

rm(list = ls())
# setwd to root
setwd("/Users/daisyzhu/Documents/Research Projects/SPADE/manuscript/code_submission/mixture_survival")
# setwd
setwd("synthetic_data/")
# read specification
load("Taiwan_X_spec.rda")
load("MLE_result.rda")
fit0 <- fit

# set seed
seed <- 47 # change this seed to generate a different synthetic dataset
set.seed(seed)

# generate covariates according to X_freq
data <- do.call(cbind, lapply(X_freq, function(x) {
      return(rbinom(n, 1, x))
}))
data <- data.frame(data)
colnames(data) <- unique(c(varlist_long, varlist_short))
data$add <- data$gs2 * data$age_old
colnames(data)[ncol(data)] <- varlist_shortint

# create X matrices, code from data analysis
X0 <- cbind(1, matrix(unlist(data[, varlist_long]), nrow = n))
X1 <- cbind(1, matrix(unlist(data[, varlist_short]), nrow = n),
            matrix(unlist(data[, varlist_shortint]), nrow = n))
X2 <- X1 
p0 <- ncol(X0)
p1 <- ncol(X1)
p2 <- ncol(X2)

# use fit0$par to generate event data under "correct model"
t0 <- 1
tau <- 365
theta0 <- c(exp(X0 %*% fit0$par[1:p0]))
theta1 <- c(exp(X1 %*% fit0$par[1:p1 + p0]))
theta2 <- c(exp(X2 %*% fit0$par[1:p2 + p0 + p1]))
dist_fun <- function(y, theta0, theta1, theta2) {
      result <- theta0 * y + 
            theta1 * ((y+t0)^(1-theta2)/(1-theta2) - 
                            t0^(1-theta2)/(1-theta2))
      return(result)
}
# generate uniform
U <- runif(n)
Y <- sapply(1:n, function(i) {
      tempfun <- function(y) {
            return(dist_fun(y, theta0[i], theta1[i], theta2[i])-U[i])
      }
      if (tempfun(tau) < 0) {
            result <- tau + 1 # this is not the real event date
                              # only set up so that it will be censored
      } else {
            fit <- uniroot(tempfun, c(0, 365))
            result <- fit$root
      }
      return(result)
})

data$EventDays <- Y
data$stroke <- 1 # placeholder
data$ID <- 1:nrow(data)
# add age_middle variable
data$age_middle <- as.numeric(data$age_young == 0 & data$age_old == 0)

save(list = "data", file = "synthetic_Taiwan.rda")