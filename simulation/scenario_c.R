rm(list = ls())
library(gdata)
source("parmix_functions.R")

# preset:
#     1. scenario_index (correspond to scenario_index in table1.csv, 1c-4c);
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
            out <- theta0 * t + theta1 * pbeta(t/60, shape1 = 2, shape2 = 5) * 60 - tempt[i]
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