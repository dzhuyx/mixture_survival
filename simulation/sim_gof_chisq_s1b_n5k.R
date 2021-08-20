# , lib.loc="/users/yzhu/Rlibs"
source("parmix_functions.R")

n <- 5000
scenario_index <- 1
nstart <- 50

seed <- as.numeric(Sys.getenv("SGE_TASK_ID"))
set.seed(seed)

X0 <- cbind(1, rbinom(n, 1, 0.5), rnorm(n))
alpha0 <- c(-8, -0.1, 0.1)
X1 <- matrix(1, nrow = n, ncol = 1)
beta10 <- -8
X2 <- matrix(1, nrow = n, ncol = 1)
beta20 <- 2
censoring <- runif(n, 200, 250)


tempt <- runif(n)
time <- sapply(1:n, function(i) {
      theta0 <- exp(X0[i, , drop = F] %*% alpha0)
      funi <- function(t) {
            out <- theta0 * t - tempt[i]
            return(out)
      }
      if (funi(0) * funi(censoring[i]) > 0) {
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
mean(d)

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
# theta0 <- fit0$par


alpha0 <- fit0$par[1:3]
beta10 <- fit0$par[4]
beta20 <- fit0$par[5]


theta00 <- c(exp(X0 %*% alpha0))
theta10 <- c(exp(X1 %*% beta10))
theta20 <- c(expit(X2 %*% beta20))
ttt <- theta00*Y + theta10 * (theta20^Y -1)/log(theta20)

# mle with no misdiagnosis group, gs = 1
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
sub.dir <- paste0("sim_gof_chisq_s", scenario_index, "b_n", n)
if(file.exists(sub.dir)) {
      setwd(file.path(main.dir, sub.dir))
} else {
      dir.create(file.path(main.dir, sub.dir))
      setwd(file.path(main.dir, sub.dir))
}
save.image(file = paste0("sim_gof_chisq_s", scenario_index, "b_n", n, "_seed", seed, ".rda"))
setwd(file.path(main.dir))
