rm(list = ls())
library(gdata)
source("parmix_functions.R")

# preset:
#     1. scenario_index (correspond to scenario_index in table1.csv, 1-8);
#     2. n (sample size = 1500, 5000, 15000);
#     3. seed (to index each replication, 1-1000)

nstart <- 50 # number of initial values to try in optimization
set.seed(seed)

table1 <- read.xls(xls = "table1.xlsx")
alpha0 <- c(table1$alpha0[table1$scenario_index == scenario_index],
            table1$alpha1[table1$scenario_index == scenario_index],
            table1$alpha2[table1$scenario_index == scenario_index])
beta10 <- table1$beta1[table1$scenario_index == scenario_index]
beta20 <- table1$beta2[table1$scenario_index == scenario_index]

X0 <- cbind(1, rbinom(n, 1, 0.5), rnorm(n))
X1 <- matrix(1, nrow = n, ncol = 1)
X2 <- matrix(1, nrow = n, ncol = 1)
censoring <- runif(n, 200, 250)

p0 <- ncol(X0)
tempt <- runif(n)
time <- sapply(1:n, function(i) {
      theta0 <- exp(X0[i, , drop = F] %*% alpha0)
      theta1 <- exp(X1[i, , drop = F] %*% beta10)
      theta2 <- expit(X2[i, , drop = F] %*% beta20)
      funi <- function(t) {
            out <- theta0 * t + theta1 / log(theta2) * (theta2^t - 1) - tempt[i]
            return(out)
      }
      if (funi(0) *funi(censoring[i]) > 0) {
            out <- censoring[i]
      } else {
            fit <- uniroot(funi, interval = c(0, censoring[i]))
            out <- fit$root
      }
      return(out)
})

# obtain MLE
Y <- sapply(1:n, function(i) return(min(time[i], censoring[i])))
d <- sapply(1:n, function(i) return(as.numeric(time[i] < censoring[i])))

# estimation
optim_fun <- function(theta) {
      out <- -llh_fun(Y = Y, d = d, X0 = X0, X1 = X1, X2 = X2, 
                      alpha = theta[1:3], beta1 = theta[4], beta2 = theta[5])
      return(out)
}

theta0 <- c(alpha0, beta10, beta20)

start <- rbind(mvrnorm(nstart, mu = theta0, Sigma = diag(c(2, 2, 2, 2, 2))), theta0)
fit <- list()
for (i in 1:nstart) {
      print(i)
      x <- start[i, ]
      fit[[i]] <- try(optim(x, fn = optim_fun, hessian = T, control = list(maxit = 50000,   reltol = 1e-8)), silent = T)
      if (!is.null(as.list(fit[[i]])$par)) {
            print(paste(fit[[i]]$convergence, fit[[i]]$value))
      } 
}
llh <- do.call("c", lapply(fit, function(x) {
      if (is.null(as.list(x)$par)) {
            out <- -1
      } else if (x$convergence == 0 & x$value != min(99999, 99 * n)) {
            out <- x$value
      } else {
            out <- -1
      }
      return(out)
}))

fit0 <- fit[[which(llh == min(llh[llh > 0]))[1]]]
p <- length(fit0$par)

alpha0 <- fit0$par[1:3]
beta10 <- fit0$par[4]
beta20 <- fit0$par[5]


theta00 <- c(exp(X0 %*% alpha0))
theta10 <- c(exp(X1 %*% beta10))
theta20 <- c(expit(X2 %*% beta20))

# mle with no misdiagnosis group
optim_fun <- function(theta) {
      out <- -llh_fun_nomis(Y = Y, d = d, X0 = X0, alpha = theta)
      return(out)
}
theta0 <- alpha0
start <- rbind(mvrnorm(nstart, mu = theta0, Sigma = diag(c(2, 2, 2))), theta0)
fit <- list()
for (i in 1:nstart) {
      print(i)
      x <- start[i, ]
      fit[[i]] <- try(optim(x, fn = optim_fun, hessian = T, control = list(maxit = 50000,   reltol = 1e-8)), silent = T)
      if (!is.null(as.list(fit[[i]])$par)) {
            print(paste(fit[[i]]$convergence, fit[[i]]$value))
      }
}
llh <- do.call("c", lapply(fit, function(x) {
      if (is.null(as.list(x)$par)) {
            out <- -1
      } else if (x$convergence == 0 & x$value != min(99999, 99 * n)) {
            out <- x$value
      } else {
            out <- -1
      }
      return(out)
}))
fit_nomis <- fit[[which(llh == min(llh[llh > 0]))[1]]]

# compute invSigma for testing
avec <- c(0, quantile(Y[d==1], c(1:5)/6))

gof0 <- gof_mle(Y, d, X0, X1, X2, fit0, avec)
gof_nomis <- gof_mle_nomis(Y = Y, d = d, X0 = X0,
                           fit0 = fit_nomis, avec = avec)


# estimate attributable harm and relative attributable harm on day 30
# setup harm measure functions
harm_measure_fun <- function(t, X, harm_fun, par) {
      # t is a scalar
      # X: list of covariate matrices
      # harm_fun is a function of t, X[[i]], and par
      q <- length(X)
      
      result <- NULL
      for (i in 1:q) {
            result <- c(result, harm_fun(t, X[[i]], par))
      }
      return(result)
}

harm_measure_est_fun <- function(t, X, harm_fun, par0, hessian0) {
      est <- harm_measure_fun(t, X, harm_fun, par0)
      
      tempfun <- function(par) {
            return(harm_measure_fun(t, X, harm_fun, par))
      }
      gradpar <- gradient(tempfun, par0)
      sigma2 <- (gradpar) %*% solve(hessian0) %*% t(gradpar)
      return(list(est = est, sigma2 = sigma2))
}

expit <- function(x) {
      return(exp(x) / (1 + exp(x)))
}

AH_fun <- function(t, X, par) {
      X1 <- X[[2]]
      X2 <- X[[3]]
      p1 <- ncol(X1)
      p2 <- ncol(X2)
      
      theta1 <- exp(X1 %*% par[1:p1+p0])
      theta2 <- expit(X2 %*% par[1:p2+p0+p1])
      
      out <- theta1 / log(theta2) * (theta2^t - 1)
      # t0 = 1, as defined earlier in this script
      return(mean(out))
}

logRAH_fun <- function(t, X, par) {
      X0 <- X[[1]]
      X1 <- X[[2]]
      X2 <- X[[3]]
      p0 <- ncol(X0)
      p1 <- ncol(X1)
      p2 <- ncol(X2)
      
      theta0 <- exp(X0 %*% par[1:p0])
      theta1 <- exp(X1 %*% par[1:p1+p0])
      theta2 <- expit(X2 %*% par[1:p2+p0+p1])
      
      baseline <- mean(theta0 * t)
      return(log(mean(AH_fun(t, X, par)) / baseline))
}

ah_30 <- harm_measure_est_fun(t = 30, 
                              X = list(list(matrix(X0, nrow = n), 
                                            matrix(X1, nrow = n), 
                                            matrix(X2, nrow = n))),
                              harm_fun = AH_fun, 
                              par0 = fit0$par, 
                              hessian0 = fit0$hessian)
lograh_30 <- harm_measure_est_fun(t = 30, 
                                  X = list(list(matrix(X0, nrow = n), 
                                                matrix(X1, nrow = n), 
                                                matrix(X2, nrow = n))),
                                  harm_fun = logRAH_fun, 
                                  par0 = fit0$par, 
                                  hessian0 = fit0$hessian)

# obtain true value
ah0_30 <- harm_measure_est_fun(t = 30, 
                               X = list(list(matrix(X0, nrow = n), 
                                             matrix(X1, nrow = n), 
                                             matrix(X2, nrow = n))),
                               harm_fun = AH_fun, 
                               par0 = theta0, 
                               hessian0 = fit0$hessian)

lograh0_30 <- harm_measure_est_fun(t = 30, 
                                   X = list(list(matrix(X0, nrow = n), 
                                                 matrix(X1, nrow = n), 
                                                 matrix(X2, nrow = n))),
                                   harm_fun = logRAH_fun, 
                                   par0 = theta0, 
                                   hessian0 = fit0$hessian)

# harm-pivotal time point for 80% of one-year attributable harm
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

tupper <- 250
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

fit_tau0_ah <- tau_estfun(X = list(list(X0, X1, X2), 
                                   list(X0, X1, X2)), 
                          harm_fun = function(tvec, X, par) {
                                return(c(AH_fun(tvec[1], X[[1]], par),
                                         AH_fun(tvec[2], X[[2]], par)))
                          }, 
                          phi = function(X, par) {
                                return(c(0.8 * AH_fun(tupper, X[[1]], par), 
                                         0.5 * AH_fun(tupper, X[[2]], par)))
                          }, 
                          par0 = theta0, 
                          hessian0 = fit0$hessian)

main.dir <- getwd()
sub.dir <- paste0("sim_gof_s", scenario_index, "_n", n)
if(file.exists(sub.dir)) {
      setwd(file.path(main.dir, sub.dir))
} else {
      dir.create(file.path(main.dir, sub.dir))
      setwd(file.path(main.dir, sub.dir))
}
save.image(file = paste0("sim_gof_s", scenario_index, "_n", n, "_seed", seed, ".rda"))
setwd(file.path(main.dir))