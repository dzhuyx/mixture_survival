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

# colMeans of c(varlist_long, varlist_short, varlist_shortint)
X_freq <- colMeans(data[, unique(c(varlist_long, varlist_short))])
save(list = c("varlist_long", "varlist_short", "varlist_shortint",
              "X_freq", "n"), file = "Taiwan_X_spec.rda")


# load analysis result rda, updated in 05/2021
load("taiwan_parmix_05192021_poly_start2189.rda") # this seed is selected according to llh
fit0 <- fit
save(list = "fit0", file = "Taiwan_par_spec.rda")
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


# calculate subgroup harm measures at 30-day
mycreatecovmat <- function(var, val) {
      result <- list(cbind(1, matrix(unlist(data[data[, var] %in% val, varlist_long]), 
                                     ncol = length(varlist_long))),
                     cbind(1, matrix(unlist(data[data[, var] %in% val, c(varlist_short, varlist_shortint)]), 
                                     ncol = length(varlist_short) + length(varlist_shortint))),
                     cbind(1, matrix(unlist(data[data[, var] %in% val, c(varlist_short, varlist_shortint)]), 
                                     ncol = length(varlist_short) + length(varlist_shortint))))
      return(result)
}
# age
X_age <- list(mycreatecovmat("age_young", 1),
              mycreatecovmat("age_middle", 1),
              mycreatecovmat("age_old", 1))
harm_ah_age <- harm_measure_est_fun(t = 30, 
                                    X = X_age,
                                    harm_fun = AH_fun, 
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_age <- harm_measure_est_fun(t = 30, 
                                     X = X_age,
                                     harm_fun = logRAH_fun, 
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

# gender groups
X_gender <- list(mycreatecovmat("gender", 0),
                 mycreatecovmat("gender", 1))
harm_ah_gender <- harm_measure_est_fun(t = 30,
                                       X = X_gender,
                                       harm_fun = AH_fun,
                                       par0 = fit0$par, 
                                       hessian0 = fit0$hessian)
harm_rah_gender <- harm_measure_est_fun(t = 30,
                                        X = X_gender,
                                        harm_fun = logRAH_fun,
                                        par0 = fit0$par, 
                                        hessian0 = fit0$hessian)

# DM
X_DM <- list(mycreatecovmat("DM", 0),
             mycreatecovmat("DM", 1))
harm_ah_DM <- harm_measure_est_fun(t = 30,
                                   X = X_DM,
                                   harm_fun = AH_fun,
                                   par0 = fit0$par, 
                                   hessian0 = fit0$hessian)
harm_rah_DM <- harm_measure_est_fun(t = 30,
                                    X = X_DM,
                                    harm_fun = logRAH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)

# before_STROKE
X_stroke <- list(mycreatecovmat("before_STROKE", 0),
                 mycreatecovmat("before_STROKE", 1))
harm_ah_stroke <- harm_measure_est_fun(t = 30,
                                       X = X_stroke,
                                       harm_fun = AH_fun,
                                       par0 = fit0$par, 
                                       hessian0 = fit0$hessian)
harm_rah_stroke <- harm_measure_est_fun(t = 30,
                                        X = X_stroke,
                                        harm_fun = logRAH_fun,
                                        par0 = fit0$par, 
                                        hessian0 = fit0$hessian)

# HT
X_HT <- list(mycreatecovmat("HT", 0),
             mycreatecovmat("HT", 1))
harm_ah_HT <- harm_measure_est_fun(t = 30,
                                   X = X_HT,
                                   harm_fun = AH_fun,
                                   par0 = fit0$par, 
                                   hessian0 = fit0$hessian)
harm_rah_HT <- harm_measure_est_fun(t = 30,
                                    X = X_HT,
                                    harm_fun = logRAH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
# gs2
X_gs2 <- list(mycreatecovmat("gs2", 0),
              mycreatecovmat("gs2", 1))
harm_ah_gs2 <- harm_measure_est_fun(t = 30,
                                    X = X_gs2,
                                    harm_fun = AH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_gs2 <- harm_measure_est_fun(t = 30,
                                     X = X_gs2,
                                     harm_fun = logRAH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

# AF
X_AF <- list(mycreatecovmat("AF", 0),
             mycreatecovmat("AF", 1))
harm_ah_AF <- harm_measure_est_fun(t = 30,
                                    X = X_AF,
                                    harm_fun = AH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_AF <- harm_measure_est_fun(t = 30,
                                     X = X_AF,
                                     harm_fun = logRAH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

# dd (hospitalization)
X_dd <- list(mycreatecovmat("dd", 0),
             mycreatecovmat("dd", 1))
harm_ah_dd <- harm_measure_est_fun(t = 30,
                                    X = X_dd,
                                    harm_fun = AH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_dd <- harm_measure_est_fun(t = 30,
                                     X = X_dd,
                                     harm_fun = logRAH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

# SES
data$SocioeconomicStatus0 <- as.numeric(data$SocioeconomicStatus1 != 1 &
                                              data$SocioeconomicStatus2 != 1)
X_SES <- list(mycreatecovmat("SocioeconomicStatus0", 1),
              mycreatecovmat("SocioeconomicStatus1", 1),
              mycreatecovmat("SocioeconomicStatus2", 1))
harm_ah_SES <- harm_measure_est_fun(t = 30,
                                    X = X_SES,
                                    harm_fun = AH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_SES <- harm_measure_est_fun(t = 30,
                                     X = X_SES,
                                     harm_fun = logRAH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

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
harm_ah_area <- harm_measure_est_fun(t = 30,
                                    X = X_area,
                                    harm_fun = AH_fun,
                                    par0 = fit0$par, 
                                    hessian0 = fit0$hessian)
harm_rah_area <- harm_measure_est_fun(t = 30,
                                     X = X_area,
                                     harm_fun = logRAH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)

# hospital type
data$Teaching1 <- as.numeric(data$Teaching0 != 1 &
                                   data$Teaching2 != 1 &
                                   data$Teaching3 != 1 &
                                   data$Teaching4 != 1)
X_teaching <- list(mycreatecovmat("Teaching0", 1),
                   mycreatecovmat("Teaching2", 1),
                   mycreatecovmat("Teaching3", 1),
                   mycreatecovmat("Teaching4", 1),
                   mycreatecovmat("Teaching1", 1))
harm_ah_teaching <- harm_measure_est_fun(t = 30,
                                     X = X_teaching,
                                     harm_fun = AH_fun,
                                     par0 = fit0$par, 
                                     hessian0 = fit0$hessian)
harm_rah_teaching <- harm_measure_est_fun(t = 30,
                                      X = X_teaching,
                                      harm_fun = logRAH_fun,
                                      par0 = fit0$par, 
                                      hessian0 = fit0$hessian)




# compute harm curves for age subgroups
harm_curve_ah_age <- do.call(rbind, lapply(1:tupper, function(t) {
      temp <- harm_measure_est_fun(t = t,
                                   X = X_age, 
                                   harm_fun = AH_fun, 
                                   par0 = fit0$par,
                                   hessian0 = fit0$hessian)
      ah <-  temp$est
      se <- sqrt(diag(temp$sigma2))
      result <- data.frame(t = t, 
                           age = c("young", "middle", "old"),
                           ah = ah, se = se,
                           lower = ah - qnorm(0.975) * se,
                           upper = ah + qnorm(0.975) * se)
      return(result)
}))

harm_curve_rah_age <- do.call(rbind, lapply(1:tupper, function(t) {
      temp <- harm_measure_est_fun(t = t,
                                   X = X_age, 
                                   harm_fun = logRAH_fun, 
                                   par0 = fit0$par,
                                   hessian0 = fit0$hessian)
      lograh <-  temp$est
      se <- sqrt(diag(temp$sigma2))
      result <- data.frame(t = t,
                           age = c("young", "middle", "old"),
                           rah = exp(lograh), 
                           logse = se,
                           lower = exp(lograh - qnorm(0.975) * se),
                           upper = exp(lograh + qnorm(0.975) * se))
      return(result)
}))

save.image(file = "taiwan_parmix_05192021_harmriskfactor.rda")

# create figures