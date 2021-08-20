# rerun survival mixture model 
# include short-term risk factors (varlist_short) in decay rate
library(rootSolve, lib.loc="/users/yzhu/Rlibs")
load("Before PSM-for Zheyu.RData")
source("parmix_function_poly.R")
# load("varlist_revised.rda")  
# load variable list
varlist_long <- c("gs2", "gender", "HT", "DM", "AF", "before_STROKE", "dd",                  
                  "SocioeconomicStatus2", "SocioeconomicStatus1", 
                  "area4", "area2", "area3", "area5",
                  "Teaching2", "age_young", "age_old")
varlist_short <- c("gs2", "gender", "HT", "DM", "before_STROKE",
                   "area4", "Teaching4", "Teaching3",
                   "Teaching0", "age_young", "age_old")

varlist_shortint <- c("gs2_age_old")

# k <- 37
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

t0 <- 1
VARS19 <- c("age", "gender", "SocioeconomicStatus", "UrbanStatus", 
            "area", "HT", "DM", "HYPERLIPIDEMIA", "MI", "HF",
            "AF", "CARDIAC_DYSRHYTHMIA", "before_STROKE", 
            "COPD", "CANCER", "CKD", "cardiologist_care", 
            "dd", "Teaching")
# data$EventDays <- data$EventDays + rnorm(nrow(data), 0, 0.2)
data$EventDays1 <- data$EventDays
data$stroke1 <- 1
data$stroke1[data$EventDays > 365] <- 0
data$EventDays1[data$EventDays > 365] <- 365
# data$age <- log(data$age)
# data$age <- scale(data$age)

data$EventDays <- data$EventDays1
data$stroke <- data$stroke1
data <- data[, c("ID", "GeneralistSpecialist", "stroke", "EventDays", VARS19)]

# plot(survfit(Surv(EventDays, stroke) ~ GeneralistSpecialist, data = data),
#      xlim = c(0, 365), ylim = c(0.99, 1), lty = c(1, 2), lwd = 3,
#      xlab = "Time (days)", ylab = "Stroke-free proportion")
# legend(10, 0.995, c("Special care", "General care"), lty = c(1, 2), bty = "n")

# data <- data[data$EventDays > 3, ]



# # # K-M plot
# Y <- data$EventDays
# d <- data$stroke
# plot(survfit(Surv(Y, d) ~ 1), ylim = c(0.99,1), xlim = c(0, 365))


cat_to_binary <- function(IDvec, x, name) {
      temp <- unique(x)
      out <- data.frame(ID = IDvec)
      colnames(out)[1] <- "ID"
      for (j in 1:(length(temp) - 1)) {
            out <- cbind(out, as.numeric(x == temp[j]))
            colnames(out)[j+1] <- paste0(name, temp[j])
      }
      return(out)
}

var_cat <- c("SocioeconomicStatus", "UrbanStatus", "area", "Teaching")
for (var in var_cat) {
      data <- merge(data, cat_to_binary(data$ID, data[, var], var), by = "ID")
      data[, var] <- NULL
}

data$age_young <- as.numeric(data$age <= 40)
data$age_middle <- as.numeric(data$age > 40 & data$age <= 60)
data$age_old <- as.numeric(data$age > 60)
data$gs2 <- as.numeric(data$GeneralistSpecialist == 2)

# create interaction variables
for (var in varlist_short) {
      temp <- as.numeric(data$GeneralistSpecialist == 2) * data[, var]
      data <- cbind(data, temp)
      colnames(data)[ncol(data)] <- paste0("gs2_", var)
}


#data$age <- scale(data$age)


gs1 <- data$GeneralistSpecialist == 1
gs2 <- data$GeneralistSpecialist == 2

# set.seed(37)
# data <- data[c(sample(which(gs1), 20000, replace = F),
#                sample(which(gs2), 20000, replace = F)), ]
gs1 <- data$GeneralistSpecialist == 1
gs2 <- data$GeneralistSpecialist == 2
n1 <- sum(gs1)
n2 <- sum(gs2)
n <- nrow(data)

# without any covariates
# p <- length(varlist)
X0 <- cbind(1, matrix(unlist(data[, varlist_long]), nrow = n))
X1 <- cbind(1, matrix(unlist(data[, varlist_short]), nrow = n),
            matrix(unlist(data[, varlist_shortint]), nrow = n))
# X2 <- matrix(1, nrow = n, ncol = 1)
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


# load("est.rda")
# temphess <- hessian(optim_fun, est)

# # DEoptim
# fit_DE <- DEoptim(optim_fun, 
#                   lower = c(-20, rep(-2, p), -20, 0, -2, -2), 
#                   upper = c(-10, rep(-1, p), -10, 0, 2, 2),
#                   control = DEoptim.control(trace = T))

set.seed(37)
start <- rbind(mvrnorm(6000, mu = theta0, Sigma = diag(rep(2, p0+p1+p2))), theta0)

x <- start[k, ]
fit <- try(optim(x, fn = optim_fun, hessian = T, 
                 control = list(maxit = 500000, reltol = 1e-80, trace = 1)), silent = T)


if (fit$value != 999999999) {
      main.dir <- getwd()
      sub.dir <- paste0("taiwan_parmix_05192021_poly")
      if(file.exists(sub.dir)) {
            setwd(file.path(main.dir, sub.dir))
      } else {
            dir.create(file.path(main.dir, sub.dir))
            setwd(file.path(main.dir, sub.dir))
      }
      save(list = c("fit"), file = paste0("taiwan_parmix_05192021_poly_start", k, ".rda"))
      setwd(file.path(main.dir))
}


# # summarize result from cluster
# setwd("parmix/taiwan_parmix_05192021_poly")
# filenames <- list.files(path = getwd())
# llh <- do.call("c", lapply(filenames, function(x) {
#    load(x)
#    return(fit$value)
# }))

####################################
# after summarizing from cluster
load("taiwan_parmix_05192021_poly_start2189.rda")
fit0 <- fit
est <- fit0$par
est

# load data and optim_fun first before running this line
# temphess <- hessian(optim_fun, est)

se <- sqrt(diag(solve(fit0$hessian)))
pval <- (est + qnorm(0.975) * se) * (est - qnorm(0.975) * se) > 0
# n <- nrow(data)
# 
# 
# long_order <- c(15, 16, 2:7, 9, 8, 10:13, 1, 14)
# short_order <- c(10, 11, 2:6, 1, 7:9)
# total_order <- c(1, long_order + 1,
#                  length(long_order) + 2,
#                  short_order + length(long_oder) + 2,
#                  30:32)
# round(est, 2)[total_order]
# round(se, 2)[total_order]
# round(est + qnorm(0.975) * se, 2)[total_order]
# round(est - qnorm(0.975) * se, 2)[total_order]
# sapply(round(pnorm(est/se), 2), function(x) {
#       if (x > 0.5) {
#             return(1 - x)
#       } else {
#             return(x)
#       }
# })[total_order] * 2



alpha <- est[1:p0]
beta1 <- est[1:p1+p0]
beta2 <- est[1:p2+p0+p1]

# n <- nrow(data)
# X0 <- matrix(unlist(X0), nrow = n)
# X1 <- matrix(unlist(X1), nrow = n)
# X2 <- matrix(unlist(X2), nrow = n)

theta0 <- exp(X0 %*% alpha)
theta1 <- exp(X1 %*% beta1)
theta2 <- exp(X2 %*% beta2)

theta3 <- rep(t0, nrow(X0))

surv3 <- NULL
dens3 <- NULL
for (t in 0:364) {
      surv3 <- c(surv3, mean(sapply(1:n, function(i) {
            return(1 - (theta0[i] * t + theta1[i] / (1-theta2[i]) * ((t+theta3[i])^(1-theta2[i]) - theta3[i]^(1-theta2[i]))))
      })))
      dens3 <- c(dens3, mean(sapply(1:n, function(i) {
            return(theta0[i] + theta1[i] * (t + theta3[i]) ^ (-theta2[i]))
      })))
}

# densfun <- function(t) {
#       dens3 <- mean(sapply(1:n, function(i) {
#             return(theta0[i] + theta1[i] * (t + theta3[i]) ^ (-theta2[i]))
#       }))
#       return(dens3)
# }


Y <- data$EventDays
d <- data$stroke

library(survival)
fit_km <- survfit(Surv(data$EventDays, data$stroke) ~ 1)
fit_km_young <- survfit(Surv(data$EventDays, data$stroke) ~ 1, subset = data$age_young==1)
jpeg(filename = "temp.jpeg", width = 1000, height = 700)
plot(c(0, fit_km$time), c(1, fit_km$surv), ylim = c(0.99,1), xlim = c(0, 365), type = "l")
lines(0:364, surv3, col = "red")
dev.off()


gof_mle_poly(data$EventDays, data$stroke,
             X0 = X0, X1 = X1, X2 = X2, fit0 = fit0,
             avec = c(0, quantile(data$EventDays[data$stroke == 1], 1:4/4)),
             t0 = t0)
# gof result: chisqstat = 7.597443, p-value = 0.1074884