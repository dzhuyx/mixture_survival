library(rootSolve)
# run harm_measure_revised.R before running this script
source("harm_function.R")
tupper <- 365
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

k <- 37


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
X2 <- cbind(1, gs2)
p0 <- ncol(X0)
p1 <- ncol(X1)
p2 <- ncol(X2)

# load analysis result rda, updated in 05/2021
load("taiwan_parmix_05192021_poly_start2189.rda")
fit0 <- fit
# # make plots for coefficient estimates and CIs
# q_coef_long <- ggplot() + geom_point(aes(x = ))

# harm risk factor analysis
tupper <- 365
# need to set X0, X1, X2, p0, p1, p2 correctly before running the following codes

AH_fun <- function(theta1, theta2, t) {
      out <- theta1 / (1 - theta2) * ((t+t0)^(1-theta2) - t0^(1-theta2))
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


# # overall summary
# library(survival)
# fit <- survfit(Surv(EventDays, stroke) ~ 1, data = data)
# temp <- sapply(1:365, function(t) {
#    return(AH_observed(est, X1, X2, t))
# })
# max(temp)
# max(temp) * 0.8
# max(which(temp <= max(temp) * 0.8)) # 39
# (max(temp) - temp[39]) * 10000
# 
# max(temp) * 0.9
# max(which(temp <= max(temp) * 0.9)) # 104
# (max(temp) - temp[104]) * 10000
# 
# (min(fit$surv[fit$time <= 39]) - min(fit$surv[fit$time <= 365])) * 10000
# (min(fit$surv[fit$time <= 104]) - min(fit$surv[fit$time <= 365])) * 10000


# subgroups analysis and plots
my4cols <- c("#E7B800", "#2E9FDF", "#FC4E07", "#66cd00")
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
my2cols <- c("#2E9FDF", "#FC4E07")
est <- fit0$par
se <- sqrt(diag(solve(fit0$hessian)))
est_lower <- est - qnorm(0.975) * se
est_upper <- est + qnorm(0.975) * se
signif <- (est + qnorm(0.975) * se) * (est - qnorm(0.975) * se) > 0

AH_observed_subgroup <- function(par, X1, X2, tupper, group_ind) {
      out <- sapply(1:tupper, function(t) {
            return(AH_observed(par = par, X1 = X1[group_ind, ], 
                               X2 = X2[group_ind, , drop = F], t))
      })
      return(out)
}

AAH_observed_subgroup <- function(par, X0, X1, X2, tupper, group_ind) {
      out <- sapply(1:tupper, function(t) {
            return(AAH_observed(par = par, X0 = X0[group_ind, ],
                                X1 = X1[group_ind, ], X2 = X2[group_ind, , drop = F], t))
      })
      return(out)
}

c(varlist_short, varlist_shortint)[signif[c(19:30)]]
c(varlist_short, varlist_shortint)[signif[c(32:43)]]

# attributable harm
result <- data.frame()
for (gs in 1) {
      # age <= 40
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_young == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "upper"))
      
      
      
      # age 40-60
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_middle == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "upper"))
      
      
      # age > 60
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_old == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "upper"))
      
      # gender female
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$gender == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "upper"))
      
      
      # gender male
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$gender == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "upper"))
      
      
      
      
      # DM = 0
      group_ind <- which(data$GeneralistSpecialist %in% c(1, 2) & data$DM == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "upper"))
      
      # DM = 1
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$DM == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "upper"))
      
      # before_stroke = 0
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$before_STROKE == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "upper"))
      
      
      # before_stroke = 1
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$before_STROKE == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "upper"))
      
      # gs = 2
      
      group_ind <- which(data$GeneralistSpecialist == 2)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "upper"))
      
      # gs = 1
      group_ind <- which(data$GeneralistSpecialist == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "upper"))
      
      
      # HT = 0
      group_ind <- which(data$HT == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "upper"))
      
      # HT = 1
      group_ind <- which(data$HT == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "upper"))
      
      # area4 = 1
      group_ind <- which(data$area4 == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "upper"))
      
      # area4 = 0
      group_ind <- which(data$area4 == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "upper"))
      
      # Teaching = 0
      group_ind <- which(data$Teaching0 == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "upper"))
      
      # Teaching = 1, 2
      group_ind <- which(data$Teaching0 == 0 & data$Teaching3 == 0 & data$Teaching4 == 0)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "upper"))
      
      # Teaching = 3
      group_ind <- which(data$Teaching3 == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "upper"))
      
      
      # Teaching = 4
      group_ind <- which(data$Teaching4 == 1)
      tempfun <- function(par) {
            return(AH_observed_subgroup(par, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "upper"))
      
      
      
}
result_AH <- result

# relative attributable harm
result <- data.frame()
for (gs in 1) {
      # age <= 40
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_young == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_ageyoung"), tupper),
                                         est = "upper"))
      
      
      
      # age 40-60
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_middle == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_agemiddle"), tupper),
                                         est = "upper"))
      
      
      # age > 60
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$age_old == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_ageold"), tupper),
                                         est = "upper"))
      
      # gender female
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$gender == 0)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_genderfemale"), tupper),
                                         est = "upper"))
      
      
      # gender male
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$gender == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gendermale"), tupper),
                                         est = "upper"))
      
      
      
      
      # DM = 0
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$DM == 0)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_DM0"), tupper),
                                         est = "upper"))
      
      # DM = 1
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$DM == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_DM1"), tupper),
                                         est = "upper"))
      
      # before_stroke = 0
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$before_STROKE == 0)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_beforestroke0"), tupper),
                                         est = "upper"))
      
      
      # before_stroke = 1
      group_ind <- which(data$GeneralistSpecialist %in% c(1,2) & data$before_STROKE == 1)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_beforestroke1"), tupper),
                                         est = "upper"))
      
      
      # gs = 2
      group_ind <- which(data$GeneralistSpecialist == 2)
      tempfun <- function(par) {
            return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
            out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                        %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gs2"), tupper),
                                         est = "upper"))
      
      
      # gs = 1
      group_ind <- which(data$GeneralistSpecialist == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_gs1"), tupper),
                                         est = "upper"))
      
      # HT = 0
      group_ind <- which(data$HT == 0)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_HT0"), tupper),
                                         est = "upper"))
      
      
      # HT = 1
      group_ind <- which(data$HT == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_HT1"), tupper),
                                         est = "upper"))
      
      # area4 = 0
      group_ind <- which(data$area4 == 0)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_area40"), tupper),
                                         est = "upper"))
      
      # area4 = 1
      group_ind <- which(data$area4 == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_area41"), tupper),
                                         est = "upper"))
      
      # Teaching = 0
      group_ind <- which(data$Teaching0 == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching0"), tupper),
                                         est = "upper"))
      
      # Teaching = 1, 2
      group_ind <- which(data$Teaching0 == 0 & data$Teaching3 == 0 & data$Teaching4 == 0)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching12"), tupper),
                                         est = "upper"))
      
      
      # Teaching = 3
      group_ind <- which(data$Teaching3 == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching3"), tupper),
                                         est = "upper"))
      
      
      # Teaching = 4
      group_ind <- which(data$Teaching4 == 1)
      tempfun <- function(par) {
         return(AAH_observed_subgroup(par, X0, X1, X2, tupper, group_ind))
      }
      temp_est <- tempfun(fit0$par)
      gradtemp <- gradient(tempfun, fit0$par)
      temp_se <- sapply(1:tupper, function(t) {
         out <- sqrt(t(gradtemp[t, ]) %*% solve(fit0$hessian) 
                     %*% (gradtemp[t, ]) / temp_est[t]^2)
      })
      temp_lower <- exp(log(temp_est) + qnorm(0.025) * temp_se)
      temp_upper <- exp(log(temp_est) + qnorm(0.975) * temp_se)
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_est,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "est"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_lower,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "lower"))
      result <- rbind(result, data.frame(t = 1:tupper, ah = temp_upper,
                                         group = rep(paste0("gs", gs, "_Teaching4"), tupper),
                                         est = "upper"))
      
      
}
result_AAH <- result


save.image(file = "taiwan_parmix_05192021_poly_harmriskfactor.rda")

library(beepr)
beep()

temp <- levels(factor(result_AH$group))
result_AH$group <- factor(result_AH$group, levels = c(temp[3], temp[1], temp[2], 
                                                      temp[4:19]))
result_AAH$group <- factor(result_AAH$group, levels = c(temp[3], temp[1], temp[2], 
                                                        temp[4:19]))


library(ggplot2)
# make revised plots
group_name <- "gs1_age"
colorname <- "Age"
labelnames <- c("<= 40", "40 - 60", "> 60")



group_vec <- unique(result_AH$group)[grepl(group_name, unique(result_AH$group))]
dat_plot <- result_AH[result_AH$group %in% group_vec, ]
dat_plot$ah <- dat_plot$ah * 10000
q1 <- ggplot(data = dat_plot) + 
      geom_line(aes(x = t, y = ah, 
                    group = interaction(group, est),
                    color = group, linetype = est), size = 1) + 
      theme_bw() + 
      scale_color_manual(name = colorname, labels = labelnames, values = my3cols) + 
      scale_linetype_discrete(guide = "none") + 
      xlab("Time (days)") + ylab("Harmed cases per 10,000 patients") + 
      xlim(1, 365) + ylim(0, 100) 

group_vec <- unique(result_AAH$group)[grepl(group_name, unique(result_AAH$group))]
dat_plot <- result_AAH[result_AAH$group %in% group_vec, ]
dat_plot$ah <- log(dat_plot$ah)
q2 <- ggplot(data = dat_plot) + 
      geom_line(aes(x = t, y = ah, 
                    group = interaction(group, est),
                    color = group, linetype = est), size = 1) + 
      theme_bw() + 
      scale_color_manual(name = colorname, labels = labelnames, values = my3cols) + 
      scale_linetype_discrete(guide = "none") + 
      xlab("Time (days)") + ylab("Relative attribuable harm (log)") + 
      xlim(1, 365) 

library(patchwork)
(((q1 ) / (q2 )) & theme(legend.position = "bottom")) + 
      plot_layout(guides = "collect")


library(scales)
# make additional histograms
# AH at 30 days for gs1
dat_plot <- result_AH[result_AH$t == 30 & grepl("gs1", result_AH$group), ]
dat_plot <- do.call(rbind, lapply(split(dat_plot, dat_plot$group), function(x) {
      out <- x[1, ]
      out$ci_upper <- x$ah[x$est == "upper"] - x$ah[x$est == "est"]
      out$ci_lower <- x$ah[x$est == "est"] - x$ah[x$est == "lower"]
      return(out)
}))
dat_plot$group <- as.character(dat_plot$group)
dat_plot$val <- sapply(dat_plot$group, function(x) {
      x <- strsplit(x, "_")[[1]][2]
      if (grepl("gs1", x)) {
         out <- 0
      } else if (grepl("gs2", x)) {
         out <- 1
      } else if (grepl("Teaching", x)) {
         if (grepl("0", x)) {
            out <- 0
         } else if (grepl("12", x)) {
            out <- 1
         } else if (grepl("3", x)) {
            out <- 2
         } else {
            out <- 3
         }
      } else if (grepl("0", x) | grepl("young", x) | grepl("female", x)) {
            out <- 0
      } else if (grepl("1", x) | grepl("mid", x) | grepl("gendermale", x)) {
            out <- 1
      } else {
            out <- 2
      }
      return(out)
})

dat_plot$var <- sapply(dat_plot$group, function(x) {
      if (grepl("age", x)) {
            out <- "Age"
      } else if (grepl("gender", x)) {
            out <- "Gender male"
      } else if (grepl("DM", x)) {
            out <- "Diabetes mellitus"
      } else if (grepl("stroke", x)) {
            out <- "Stroke history"
      } else if (grepl("HT", x)){
            out <- "Hypertension"
      } else if(grepl("SocioeconomicStatus2", x)) {
            out <- "Moderate SES"
      } else if (grepl("area4", x)) {
            out <- "Eastern Taiwan"
      } else if (grepl("Teaching", x)) {
         out <- "Hospital Type"
      } else {
         out <- "Generalist Hospital"
      }
      
      return(out)
})

dat_plot$var <- factor(dat_plot$var, levels = c("Age", "Gender male", 
                                                "Hypertension", 
                                                "Diabetes mellitus", 
                                                "Stroke history", "Moderate SES", 
                                                "Eastern Taiwan",
                                                "Hospital Type",
                                                "Generalist Hospital"))




q5 <- ggplot(data = dat_plot, 
             aes(x = var, y = ah * 10000, fill = factor(val))) + 
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) + 
      geom_errorbar(aes(ymin = (ah - ci_lower) * 10000,
                        ymax = (ah + ci_upper) * 10000,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) + 
      scale_fill_manual(name = "Group", values = my4cols, 
                        labels = c("<=40 / No / Unclassified", 
                                   "40-60 / Yes / Medical Center or Regional Hospital", 
                                   "> 60 / Local Hospital",
                                   "Private Practice")) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 50), oob = rescale_none) +
      xlab("") + ylab("Harmed cases per 10,000 patients")




# AAH at 30 days for gs1
dat_plot <- result_AAH[result_AAH$t == 30 & grepl("gs1", result_AAH$group), ]
dat_plot <- do.call(rbind, lapply(split(dat_plot, dat_plot$group), function(x) {
      out <- x[1, ]
      out$ci_upper <- x$ah[x$est == "upper"] - x$ah[x$est == "est"]
      out$ci_lower <- x$ah[x$est == "est"] - x$ah[x$est == "lower"]
      return(out)
}))
dat_plot$group <- as.character(dat_plot$group)
dat_plot$val <- sapply(dat_plot$group, function(x) {
   x <- strsplit(x, "_")[[1]][2]
   if (grepl("gs1", x)) {
      out <- 0
   } else if (grepl("gs2", x)) {
      out <- 1
   } else if (grepl("Teaching", x)) {
      if (grepl("0", x)) {
         out <- 0
      } else if (grepl("12", x)) {
         out <- 1
      } else if (grepl("3", x)) {
         out <- 2
      } else {
         out <- 3
      }
   } else if (grepl("0", x) | grepl("young", x) | grepl("female", x)) {
      out <- 0
   } else if (grepl("1", x) | grepl("mid", x) | grepl("gendermale", x)) {
      out <- 1
   } else {
      out <- 2
   }
   return(out)
})

dat_plot$var <- sapply(dat_plot$group, function(x) {
   if (grepl("age", x)) {
      out <- "Age"
   } else if (grepl("gender", x)) {
      out <- "Gender male"
   } else if (grepl("DM", x)) {
      out <- "Diabetes mellitus"
   } else if (grepl("stroke", x)) {
      out <- "Stroke history"
   } else if (grepl("HT", x)){
      out <- "Hypertension"
   } else if(grepl("SocioeconomicStatus2", x)) {
      out <- "Moderate SES"
   } else if (grepl("area4", x)) {
      out <- "Eastern Taiwan"
   } else if (grepl("Teaching", x)) {
      out <- "Hospital Type"
   } else {
      out <- "Generalist Hospital"
   }
   
   return(out)
})

dat_plot$var <- factor(dat_plot$var, levels = c("Age", "Gender male", 
                                                "Hypertension", 
                                                "Diabetes mellitus", 
                                                "Stroke history", "Moderate SES", 
                                                "Eastern Taiwan",
                                                "Hospital Type",
                                                "Generalist Hospital"))



q7 <- ggplot(data = dat_plot, 
             aes(x = var, y = (ah), fill = factor(val))) + 
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) + 
      geom_errorbar(aes(ymin = (ah - ci_lower),
                        ymax = (ah + ci_upper),
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) + 
      scale_fill_manual(name = "Group", values = my4cols, 
                        labels = c("<=40 / No / Unclassified", 
                                   "40-60 / Yes / Medical Center or Regional Hospital", 
                                   "> 60 / Local Hospital",
                                   "Private Practice")) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 20), oob = rescale_none) +
      xlab("") + ylab("Relative attributable harm")

((((q5) | (q7)) & theme(legend.position = "bottom")) + 
      plot_layout(guides = "collect")) # 1600 * 500

((((q1 ) | (q2 )) & theme(legend.position = "bottom")) + 
   plot_layout(guides = "collect"))
