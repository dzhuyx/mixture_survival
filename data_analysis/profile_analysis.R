# profile analyses
# for new model run on 05/19/2021
# run taiwan_parmix_finalv7_harmgs.R instead
rm(list = ls())
library(rootSolve)
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

# k <- 3
# k <- as.numeric(Sys.getenv("SGE_TASK_ID"))

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
X2 <- X1 # new analyses updated on 05/19/2021
# X2 <- cbind(1, gs2) # old analyses covariate matrix
p0 <- ncol(X0)
p1 <- ncol(X1)
p2 <- ncol(X2)


# load analysis result rda, updated in 05/2021
load("taiwan_parmix_05192021_poly_start2189.rda")
fit0 <- fit
tupper <- 365

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

logRAH_fun <- function(t, X, par) {
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
      return(log(mean(AH_fun(t, X, par)) / baseline))
}

# create 6 types of covariate matrix list
mycreatecovmat <- function(var, val, subvar = NULL, subval = NULL) {
      # var and val are scalars
      # subvar is a vector, subval is a list of the same length
      # each element of subval is a vector of length n = sample size
      X0 <- cbind(1, matrix(unlist(data[data[, var] %in% val, varlist_long]), 
                            ncol = length(varlist_long)))
      X1 <- cbind(1, matrix(unlist(data[data[, var] %in% val, c(varlist_short, varlist_shortint)]), 
                            ncol = length(varlist_short) + length(varlist_shortint)))
      if (!is.null(subvar)) {
            varlist2 <- c(varlist_short, varlist_shortint)
            for (i in 1:length(subvar)) {
                  if (subvar[i] %in% varlist_long) {
                        X0[, which(varlist_long == subvar[i]) + 1] <- subval[[i]][data[, var] %in% val]
                  }
                  if (subvar[i] %in% varlist2) {
                        X1[, which(varlist2 == subvar[i]) + 1] <- subval[[i]][data[, var] %in% val]
                  }
            }
      }
      X2 <- X1
      result <- list(X0, X1, X2)
      return(result)
}
# observed populations
X_gs1treatgs1 <- mycreatecovmat(var = "gs2", val = 0)
X_gs2treatgs2 <- mycreatecovmat(var = "gs2", val = 1)
# standard populations
X_gs2treatgs1 <- mycreatecovmat(var = "gs2", val = 0, 
                                subvar = c("gs2", "gs2_age_old"),
                                subval = list(rep(1, nrow(data)),
                                              rep(1, nrow(data)) * data$age_old))
X_gs1treatgs2 <- mycreatecovmat(var = "gs2", val = 1,
                                subvar = c("gs2", "gs2_age_old"),
                                subval = list(rep(0, nrow(data)), 
                                              rep(0, nrow(data)) * data$age_old))
# adjusted populations
X_gs2treatall <- mycreatecovmat(var = "gs2", val = c(0, 1),
                                subvar = c("gs2", "gs2_age_old"),
                                subval = list(rep(1, nrow(data)),
                                              rep(1, nrow(data)) * data$age_old))
X_gs1treatall <- mycreatecovmat(var = "gs2", val = c(0, 1),
                                subvar = c("gs2", "gs2_age_old"),
                                subval = list(rep(0, nrow(data)), 
                                              rep(0, nrow(data)) * data$age_old))


# attributable harm curves for six populations
harm_curve_ah_profile <- do.call(rbind, lapply(1:tupper, function(t) {
      temp <- harm_measure_est_fun(t = t, 
                                   X = list(X_gs1treatgs1,
                                            X_gs2treatgs2,
                                            X_gs1treatgs2, 
                                            X_gs2treatgs1,
                                            X_gs1treatall,
                                            X_gs2treatall),
                                   harm_fun = AH_fun, 
                                   par0 = fit0$par, 
                                   hessian0 = fit0$hessian)
      # values that can be obtained directly
      ah <- temp$est
      se <- sqrt(diag(temp$sigma2))
      
      # gs1standard, average of gs1treatgs1 and gs2treatgs1
      ah_gs1standard <- (temp$est[1] + temp$est[4]) / 2
      se_gs1standard <- sqrt(t(c(1/2, 1/2)) %*% temp$sigma2[c(1, 4), c(1, 4)] %*%
            c(1/2, 1/2))
      
      # gs2standard, average of gs1treatgs2 and gs2treatgs2
      ah_gs2standard <- (temp$est[2] + temp$est[3]) / 2
      se_gs2standard <- sqrt(t(c(1/2, 1/2)) %*% temp$sigma2[c(2, 3), c(2, 3)] %*%
            c(1/2, 1/2))
      
      # add standard performance measures to vector
      ah <- c(ah, ah_gs1standard, ah_gs2standard)
      se <- c(se, se_gs1standard, se_gs2standard)
      
      result <- data.frame(t = rep(t, 24),
                           type = rep(c("gs1treatgs1", "gs2treatgs2",
                                    "gs1treatgs2", "gs2treatgs1",
                                    "gs1treatall", "gs2treatall", 
                                    "gs1standard", "gs2standard"), 3),
                           ah = c(ah, 
                                  ah - qnorm(0.975) * se,
                                  ah + qnorm(0.975) * se),
                           type = rep(c("est", "lower", "upper"), each = 8))
}))

# relative attributable harm functions for six populations
harm_curve_rah_profile <- do.call(rbind, lapply(1:tupper, function(t) {
      temp <- harm_measure_est_fun(t = t, 
                                   X = list(X_gs1treatgs1,
                                            X_gs2treatgs2,
                                            X_gs1treatgs2, 
                                            X_gs2treatgs1,
                                            X_gs1treatall,
                                            X_gs2treatall),
                                   harm_fun = logRAH_fun, 
                                   par0 = fit0$par, 
                                   hessian0 = fit0$hessian)
      lograh <- temp$est
      logse <- sqrt(diag(temp$sigma2)) # se of lograh
      
      # gs1standard, "average" of gs1treatgs1 and gs2treatgs1
      lograh_gs1standard <- log((exp(lograh[1]) + exp(lograh[4])))
      tempgrad <- c(exp(lograh[1]), exp(lograh[4])) / 
            (exp(lograh[1]) + exp(lograh[4]))
      logse_gs1standard <- sqrt(t(tempgrad) %*% temp$sigma2[c(1, 4), c(1, 4)] %*%
                                      tempgrad)
      
      # gs2standard, "average" of gs1treatgs2 and gs2treatgs2
      lograh_gs2standard <- log((exp(lograh[2]) + exp(lograh[3])))
      tempgrad <- c(exp(lograh[2]), exp(lograh[3])) / 
            (exp(lograh[2]) + exp(lograh[3]))
      logse_gs2standard <- sqrt(t(tempgrad) %*% temp$sigma2[c(2, 3), c(2, 3)] %*%
                                      tempgrad)
      
      # add standard performance measures to vector
      lograh <- c(lograh, lograh_gs1standard, lograh_gs2standard)
      logse <- c(logse, logse_gs1standard, logse_gs2standard)
      
      result <- data.frame(t = rep(t, 24),
                           type = rep(c("gs1treatgs1", "gs2treatgs2",
                                    "gs1treatgs2", "gs2treatgs1",
                                    "gs1treatall", "gs2treatall",
                                    "gs1standard", "gs2standard"), 3),
                           rah = c(exp(lograh),
                                   exp(lograh - qnorm(0.975) * logse),
                                   exp(lograh + qnorm(0.975) * logse)),
                           type = rep(c("est", "lower", "upper"), each = 8))
}))

save.image(file = "taiwan_parmix_05192021_gsprofile.rda")