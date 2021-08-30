rm(list = ls())

# calculate expected number of events and harmed cases
# calculation through simulation, numerical variation exists
set.seed(37)
expit <- function(x) {
   return(exp(x) / (1 + exp(x)))
}
logit <- function(x) {
   return(log(x / (1-x)))
}
n <- 1000000
X0 <- cbind(1, rbinom(n, 1, 0.5), rnorm(n))
X1 <- matrix(1, nrow = n, ncol = 1)
X2 <- matrix(1, nrow = n, ncol = 1)
censoring <- runif(n, 200, 250)
table1 <- read.xls(xls = "table1.xlsx")

dat_gof <- data.frame()
# add `observed` (number of events) and `elevated` (number of harmed cases)
# for scenarios 1-8 and 1b, 3b
for (scenario_index in 1:8) {
   alpha0 <- c(table1$alpha0[table1$scenario_index == scenario_index],
               table1$alpha1[table1$scenario_index == scenario_index],
               table1$alpha2[table1$scenario_index == scenario_index])
   beta10 <- table1$beta1[table1$scenario_index == scenario_index]
   beta20 <- table1$beta2[table1$scenario_index == scenario_index]
   
   # number of total events
   density_all <- sapply(1:n, function(i) {
      theta0 <- exp(X0[i, , drop = F] %*% alpha0)
      theta1 <- exp(X1[i, , drop = F] %*% beta10)
      theta2 <- expit(X2[i, , drop = F] %*% beta20)
      funi <- function(t) {
         out <- theta0 * t + theta1 / log(theta2) * (theta2^t - 1)
         return(out)
      }
      return(funi(censoring[i]))
   })
   
   # baseline_rate
   density_baseline <- sapply(1:n, function(i) {
      theta0 <- exp(X0[i, , drop = F] %*% alpha0)
      theta1 <- 0
      theta2 <- expit(X2[i, , drop = F] %*% beta20)
      funi <- function(t) {
         out <- theta0 * t
         return(out)
      }
      return(funi(censoring[i]))
   })
   
   dat_gof <- rbind(dat_gof, data.frame(scenario = rep(scenario_index, 3),
                                        sample_size = c(1500, 5000, 15000),
                                        observed = mean(density_all) * 
                                           c(1500, 5000, 15000),
                                        elevated = (mean(density_all) - mean(density_baseline)) *
                                           c(1500, 5000, 15000)))
   if (scenario_index %in% c(1, 3)) {
      dat_gof <- rbind(dat_gof, data.frame(scenario = rep(paste0(scenario_index, "b"), 3),
                                           sample_size = c(1500, 5000, 15000),
                                           observed = mean(density_baseline) * 
                                              c(1500, 5000, 15000),
                                           elevated = rep(NA, 3)))
   }
}
# add numbers for 1c-4c
for (scenario_index in paste0(1:4, "c")) {
   alpha0 <- c(table1$alpha0[table1$scenario_index == scenario_index],
               table1$alpha1[table1$scenario_index == scenario_index],
               table1$alpha2[table1$scenario_index == scenario_index])
   beta10 <- table1$beta1[table1$scenario_index == scenario_index]
   beta20 <- table1$beta2[table1$scenario_index == scenario_index]
   
   # number of total events
   density_all <- sapply(1:n, function(i) {
      theta0 <- exp(X0[i, , drop = F] %*% alpha0)
      theta1 <- exp(X1[i, , drop = F] %*% beta10)
      theta2 <- expit(X2[i, , drop = F] %*% beta20)
      funi <- function(t) {
         out <- theta0 * t + theta1 * pbeta(t/60, shape1 = 2, shape2 = 5) * 60
         return(out)
      }
      return(funi(censoring[i]))
   })
   
   dat_gof <- rbind(dat_gof, data.frame(scenario = rep(scenario_index, 3),
                                        sample_size = c(1500, 5000, 15000),
                                        observed = mean(density_all) * 
                                           c(1500, 5000, 15000),
                                        elevated = rep(NA, 3)))
}




# summarize simulation gof results
main.dir <- getwd()
# summarize gof for scenarios 1-8
out_a <- data.frame()
for (scenario in as.character(1:8)) {
      for (sample_size in c(1500, 5000, 15000)) {
            print(paste0("scenario:", scenario,",\ sample size:", sample_size))
            setwd(paste0(main.dir, "/sim_gof_s", scenario, "_n", sample_size))
            filenames <- list.files(path = getwd())
            temp <- do.call(rbind, lapply(filenames, function(x) {
                  load(x)
                  return(c(gof0$pval, gof_nomis$pval)) # type-I error and/or power
                                                       # could contain NAs
            }))
            out_a <- rbind(out_a, 
                           data.frame(scenario = scenario,
                                      sample_size = sample_size,
                                      error = sum(temp[which(!is.na(temp[, 1]))[1:1000], 1] <= 0.05) / 1000,
                                      power = sum(temp[which(!is.na(temp[, 2]))[1:1000], 2] <= 0.05) / 1000))
      }
}
setwd(main.dir)

# summarize gof for scenarios 1b and 3b
out_b <- data.frame()
for (scenario in c("1b", "3b")) {
   for (sample_size in c(1500, 5000, 15000)) {
      print(paste0("scenario:", scenario,",\ sample size:", sample_size))
      setwd(paste0(main.dir, "/sim_gof_s", scenario, "_n", sample_size))
      filenames <- list.files(path = getwd())
      temp <- do.call(rbind, lapply(filenames, function(x) {
         load(x)
         return(c(gof0$pval, gof_nomis$pval)) # type-I error and/or power
         # could contain NAs
      }))
      out_b <- rbind(out_b, 
                     data.frame(scenario = scenario,
                                sample_size = sample_size,
                                error = sum(temp[which(!is.na(temp[, 1]))[1:1000], 1] <= 0.05) / 1000,
                                power = sum(temp[which(!is.na(temp[, 2]))[1:1000], 2] <= 0.05) / 1000))
   }
}
setwd(main.dir)

# summarize gof for 1c-4c
out_c <- data.frame()
for (scenario in paste0(1:4, "c")) {
   for (sample_size in c(1500, 5000, 15000)) {
      print(paste0("scenario:", scenario,",\ sample size:", sample_size))
      setwd(paste0(main.dir, "/sim_gof_s", scenario, "_n", sample_size))
      filenames <- list.files(path = getwd())
      temp <- do.call(rbind, lapply(filenames, function(x) {
         load(x)
         return(c(gof0$pval, gof_nomis$pval)) # type-I error and/or power
         # could contain NAs
      }))
      out_c <- rbind(out_c, 
                     data.frame(scenario = scenario,
                                sample_size = sample_size,
                                error = sum(temp[which(!is.na(temp[, 1]))[1:1000], 1] <= 0.05) / 1000,
                                power = sum(temp[which(!is.na(temp[, 2]))[1:1000], 2] <= 0.05) / 1000))
   }
}
setwd(main.dir)

# merge with descriptive numbers
dat_gof <- merge(dat_gof, rbind(out_a, out_b, out_c), 
                 by = c("scenario", "sample_size"),
                 all.x = T, all.y = T)
dat_gof$scenario <- factor(dat_gof$scenario, 
                           levels = c(as.character(1:8), "1b", "3b",
                                      paste0(1:4, "c")))
dat_gof <- dat_gof[order(dat_gof$scenario), ]
dat_gof$scenario <- as.character(dat_gof$scenario)


# create Figure 1
library(gdata)
library(ggplot2)

dat_gof <- read.xls(xls = "gof_simulation_setting.xlsx", sheet = 3)
dat_gof$logobserved <- log(dat_gof$observed)
dat_gof$logelevated <- log(dat_gof$elevated)

dat_gof$label <- sapply(1:nrow(dat_gof), function(i) {
   return(paste0(dat_gof$scenario[i]))
})

q3 <- ggplot(data = dat_gof[25:30, ], aes(x = logobserved, y = type.I.error, label = label)) + 
   geom_text(hjust = 1, vjust = -0.5, size = 3) +
   geom_point(aes(x = logobserved, y = type.I.error, 
                  group = factor(sample.size), 
                  shape = factor(sample.size), 
                  color = factor(sample.size)),
              size = 6, alpha = 0.7) +
   theme_bw() + 
   geom_smooth(aes(x = logobserved, y = type.I.error,
                   group = factor(sample.size), 
                   color = factor(sample.size)),
               method = "loess", se = F, color = "black") + 
   geom_hline(yintercept = 0.05, linetype = "dotdash", color = "gray") +
   xlab("Log number of observed cases") +
   ylab("Type I error") +
   scale_color_manual(name = "Sample size", values = my3cols, guide = "none") + 
   scale_shape_discrete(name = "Sample size", guide = "none") + 
   ggtitle("(c) Testing existence of misdiagnosis-related harm when there is no harm.")


fit <- lm(log((power-0.001) / (1.001-power)) ~ logobserved * factor(sample.size), data = dat_gof[31:42, ])
dat_gof$predicted_power[31:42] <- (1.001 * exp(predict(fit)) + 0.01) /  (1 + exp(predict(fit)))

q4 <- ggplot(data = dat_gof[31:42, ], aes(x = logobserved, y = power, label = label)) + 
   geom_text(hjust = 1, vjust = -0.5, size = 3)+
   geom_point(aes(x = logobserved, y = power, 
                  group = factor(sample.size), 
                  shape = factor(sample.size), 
                  color = factor(sample.size)),
              size = 6, alpha = 0.7) + 
   theme_bw() + 
   geom_smooth(aes(x = logobserved, y = predicted_power,
                   group = factor(sample.size), 
                   linetype = factor(sample.size),
                   color = factor(sample.size)),
               method = "loess", se = F) + 
   xlab("Log number of observed cases") + ylab("Power") + 
   scale_color_manual(name = "Sample size", values = my3cols) + 
   scale_shape_discrete(name = "Sample size") + 
   scale_linetype_discrete(guide = "none") +
   ggtitle("(d) Testing existence of misdiagnosis-related harm when there is harm.")


dat_gof <- dat_gof[1:24, ]
fit <- lm(log((power-0.001) / (1.001-power)) ~ logelevated * factor(sample_size), data = dat_gof[1:24, ])
dat_gof$predicted_power <- (1.001 * exp(predict(fit)) + 0.01) /  (1 + exp(predict(fit)))

fit <- lm(log((error-0.05) / (0.995 - error)) ~ logobserved * factor(sample_size), 
          data = dat_gof[dat_gof$error> 0.05, ][1:24, ])
dat_gof$predicted_error <- (0.995 * exp(predict(fit, newdata = dat_gof)) + 0.05) /  (1 + exp(predict(fit, newdata = dat_gof)))


q1 <- ggplot(data = dat_gof[1:24, ], aes(x = logobserved, y = error, label = label)) + 
   geom_text(hjust = 1, vjust = -0.5, size = 3) +
   geom_point(aes(x = logobserved, y = error, 
                  group = factor(sample_size), 
                  shape = factor(sample_size), 
                  color = factor(sample_size)),
              size = 6, alpha = 0.7) +
   theme_bw() + 
   geom_smooth(aes(x = logobserved, y = predicted_error,
                   group = factor(sample_size), 
                   linetype = factor(sample_size), 
                   color = factor(sample_size)),
               method = "loess", se = F) + 
   geom_hline(yintercept = 0.05, linetype = "dotdash", color = "gray") +
   xlab("Log number of observed cases") +
   ylab("Type I error") + 
   scale_color_manual(name = "Sample size", values = my3cols) + 
   scale_shape_discrete(name = "Sample size") + 
   scale_linetype_discrete(guide = "none") +
   ggtitle("(a) Testing parametric assumption when model is correctly specified.")

q2 <- ggplot(data = dat_gof[1:24, ], aes(x = logelevated, y = power, label = label)) + 
   geom_text(hjust = 1, vjust = -0.5, size = 3)+
   geom_point(aes(x = logelevated, y = power, 
                  group = factor(sample_size), 
                  shape = factor(sample_size), 
                  color = factor(sample_size)),
              size = 6, alpha = 0.7) + 
   theme_bw() + 
   geom_smooth(aes(x = logelevated, y = predicted_power,
                   group = factor(sample_size), 
                   linetype = factor(sample_size),
                   color = factor(sample_size)),
               method = "loess", se = F) + 
   xlab("Log number of harmed cases") +
   ylab("Power") + labs(shape = "Sample size") + labs(linetype = "Sample size") +
   scale_color_manual(name = "Sample size", values = my3cols) + 
   scale_shape_discrete(name = "Sample size") + 
   scale_linetype_discrete(guide = "none")+
   ggtitle("(b) Testing parametric assumption when model is misspecified.")

library(patchwork)
jpeg(filename = "parmix_sim_gof.jpeg", width = 1300, height = 800)
(((q1+q2) / (q3+q4)) & theme(legend.position = "bottom")) + 
   plot_layout(guides = "collect")
dev.off()

# create Table 2
out <- NULL
scenario_par0 <- rbind(c(-8, -0.1, 0.1, -8, 2),
                       c(-8, -0.1, 0.1, -6, 2),
                       c(-10, -0.1, 0.1, -8, 2),
                       c(-10, -0.1, 0.1, -6, 2), 
                       c(-8, -0.1, 0.1, -8, 1),
                       c(-8, -0.1, 0.1, -6, 1),
                       c(-10, -0.1, 0.1, -8, 1),
                       c(-10, -0.1, 0.1, -6, 1))
for (scenario in 1:8) {
   for (sample_size in c(1500, 5000, 15000)) {
      print(paste0("scenario:", scenario,",\ sample size:", sample_size))
      est0 <- scenario_par0[scenario, ]
      p <- length(est0)
      setwd(paste0(main.dir, "/sim_gof_s", scenario, "_n", sample_size))
      filenames <- list.files(path = getwd())
      
      if (length(filenames) <= 1000) {
         stop("not enough replications")
      } 
      temp <- do.call(rbind, lapply(filenames, function(x) {
         load(x)
         if (sum(d) != 0) {
            # coefficient estimates
            est <- fit0$par
            se <- sqrt(diag(solve(fit0$hessian)))
            cp <- (est - qnorm(0.975) * se <= est0) & (est + qnorm(0.975) * se >= est0)
            # harm measures
            hm <- c(ah_30$est, lograh_30$est, fit_tau_ah$est[1])
            hm0 <- c(ah0_30$est, lograh0_30$est, fit_tau0_ah$est[1])
            se_hm <- sqrt(c(ah_30$sigma2, lograh_30$sigma2, fit_tau_ah$sigma[1, 1]))
            cp_hm <- (hm - qnorm(0.975) * se_hm <= hm0) & (hm + qnorm(0.975) * se_hm >= hm0)
            
            return(c(est - est0, hm - hm0, 
                     se, se_hm, 
                     cp, cp_hm, fit0$convergence))
         } else {
            print(x)
            return(NULL)
         }
      }))
      temp <- temp[complete.cases(temp) & temp[, 3*(p+3)+1] == 0, ][1:1000, ]
      out <- rbind(out, c(colMeans(temp[, 1:(p + 3)]), 
                          colMeans(temp[, 1:(p+3)+p+3]),
                          sqrt(diag(cov(temp[, 1:(p+3)]))),
                          colMeans(temp[, 1:(p+3)+2*(p+3)])))
      print(c(colMeans(temp[, 1:(p+3)]), 
              colMeans(temp[, 1:(p+3)+(p+3)]),
              sqrt(diag(cov(temp[, 1:(p+3)]))),
              colMeans(temp[, 1:(p+3)+2*(p+3)])))
   }
}
setwd(main.dir)
write.csv(out, file = "Table2_raw.csv")

# reformat Table 2 to the structure shown in manuscript
library(gdata)
# reformat simulation results in `sim_mle_s1to8_n1500to15k.csv`
# to one appropriate manuscript tables
sim_result <- read.csv("Table2_raw.csv")
sim_result <- data.frame(sim_result)
colnames(sim_result) <- c("row", paste0("bias", 1:8),
                          paste0("mse", 1:8),
                          paste0("ese", 1:8), 
                          paste0("cp", 1:8))
sim_result$scenario <- rep(1:8, each = 3)
sim_result$sample_size <- rep(c(1500, 5000, 15000), 8)

dat <- sim_result[, c(34, 35, 2:33)] # rearrange columns
dat[, 1:24+2] <- apply(dat[, 1:24+2], c(1, 2), function(x) {
   temp <- strsplit(formatC(x, format = "e", 
                            flag = "0", digits = 2), "e")[[1]]
   if (as.numeric(temp[2]) == 0) {
      out <- temp[1]
   } else {
      out <- paste0(temp[1], "\\times 10^{", 
                    as.numeric(temp[2]), "}")
   }
   return(out)
})

dat[, 25:32+2] <- apply(dat[, 25:32+2], c(1, 2), function(x) {
   out <- formatC(x, format = "f", flag = "0", digits = 2)
   return(out)
})

# output a table for each sample size
library(reshape2)
dat_bias <- melt(dat[, c(1:10)], 
                 id.vars = c("scenario", "sample_size"),
                 value.name = "bias")
dat_bias$par <- sapply(dat_bias$variable, function(x) {
   return(paste0("parameter", 
                 strsplit(as.character(x), "bias")[[1]][2]))
})
dat_bias$variable <- NULL
dat_bias <- dat_bias[, c("scenario", "sample_size",
                         "par", "bias")]

dat_mse <- melt(dat[, c(1:2, 11:18)],
                id.vars = c("scenario", "sample_size"),
                value.name = "mse")
dat_mse$par <- sapply(dat_mse$variable, function(x) {
   return(paste0("parameter", 
                 strsplit(as.character(x), "mse")[[1]][2]))
})
dat_mse$variable <- NULL

dat_ese <- melt(dat[, c(1:2, 18:26)],
                id.vars = c("scenario", "sample_size"),
                value.name = "ese")
dat_ese$par <- sapply(dat_ese$variable, function(x) {
   return(paste0("parameter", 
                 strsplit(as.character(x), "ese")[[1]][2]))
})
dat_ese$variable <- NULL

dat_cp <- melt(dat[, c(1:2, 27:34)],
               id.vars = c("scenario", "sample_size"),
               value.name = "cp")
dat_cp$par <- sapply(dat_cp$variable, function(x) {
   return(paste0("parameter", 
                 strsplit(as.character(x), "cp")[[1]][2]))
})
dat_cp$variable <- NULL

# merge bias, mse, ese, and cp
dat_out <- merge(dat_bias, dat_mse, 
                 by = c("scenario", "sample_size", "par"))
dat_out <- merge(dat_out, dat_ese,
                 by = c("scenario", "sample_size", "par"))
dat_out <- merge(dat_out, dat_cp,
                 by = c("scenario", "sample_size", "par"))

# load simulation setting numbers 
# add expected total events and harm-attributable events to table
dat_setting <- read.xls(xls = "gof_simulation_setting.xlsx",
                        sheet = 1)
dat_setting <- dat_setting[1:8, c("s", "sim_1.5k", "sim_5k", "sim_15k",
                                  "power_1.5k", "power_5k", "power_15k")]
colnames(dat_setting)[1] <- "scenario"

dat_setting_total <- melt(dat_setting[, 1:4], 
                          id.vars = "scenario", 
                          value.name = "total")
dat_setting_total$sample_size <- NA
dat_setting_total$sample_size[dat_setting_total$variable == "sim_1.5k"] <- 1500
dat_setting_total$sample_size[dat_setting_total$variable == "sim_5k"] <- 5000
dat_setting_total$sample_size[dat_setting_total$variable == "sim_15k"] <- 15000

dat_setting_harm <- melt(dat_setting[, c(1, 5:7)],
                         id.vars = "scenario",
                         value.name = "harm")
dat_setting_harm$sample_size <- NA
dat_setting_harm$sample_size[dat_setting_harm$variable == "power_1.5k"] <- 1500
dat_setting_harm$sample_size[dat_setting_harm$variable == "power_5k"] <- 5000
dat_setting_harm$sample_size[dat_setting_harm$variable == "power_15k"] <- 15000

# merge setting data with simulation results
dat_out <- merge(dat_out, dat_setting_total,
                 by = c("scenario", "sample_size"))
dat_out <- merge(dat_out, dat_setting_harm, 
                 by = c("scenario", "sample_size"))

dat_out$total <- sapply(dat_out$total, function(x) {
   # transform total and harm counts to number per 10,000 persons
   return(formatC(x/1500*10000, format = "f", digits = 1))
})
dat_out$harm <- sapply(dat_out$harm, function(x) {
   return(formatC(x/1500*10000, format = "f", digits = 1))
})
# arrange rows and output
out <- cbind(dat_out[dat_out$sample_size == 1500, c(1, 9, 11, 3:7)],
             dat_out[dat_out$sample_size == 5000, c(4:7)],
             dat_out[dat_out$sample_size == 15000, c(4:7)])

write.csv(out, file = "Table2.csv")
