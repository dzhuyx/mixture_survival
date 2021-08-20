# summarize gof simulation results
out <- NULL

my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
my2cols <- c("#2E9FDF", "#FC4E07")

for (scenario in 2) {
      for (sample_size in c(1500, 5000, 15000)) {
            print(paste0("scenario:", scenario,",\ sample size:", sample_size))
            setwd(paste0("/users/yzhu/parmix_simulation/sim_gof_chisq_s", scenario, "b_n", sample_size))
            filenames <- list.files(path = getwd())
            temp <- do.call(rbind, lapply(filenames, function(x) {
                  load(x)
                  return(c(gof0$pval, gof_nomis$pval))
            }))
            out <- rbind(out, c(sum(temp[which(!is.na(temp[, 1]))[1:1000], 1] <= 0.05) / 1000,
                                sum(temp[which(!is.na(temp[, 2]))[1:1000], 2] <= 0.05) / 1000))
            print(c(sum(temp[which(!is.na(temp[, 1]))[1:1000], 1] <= 0.05) / 1000,
                    sum(temp[which(!is.na(temp[, 2]))[1:1000], 2] <= 0.05) / 1000))
      }
}
setwd("/users/yzhu/parmix_simulation")
write.csv(out, file = "sim_gof_chisq_s1tos4_n1500to15k.csv")

library(gdata)
library(ggplot2)
# make gof simulation plots
dat_gof <- read.xls(xls = "gof_simulation_setting.xlsx", sheet = 3)
dat_gof$logobserved <- log(dat_gof$observed)
dat_gof$logelevated <- log(dat_gof$elevated)
#dat_gof <- dat_gof[1:24, ]
dat_gof$label <- sapply(1:nrow(dat_gof), function(i) {
      return(paste0(dat_gof$scenario[i]))
})

#dat_gof <- dat_gof[1:24, ]
fit <- lm(log((power-0.001) / (1.001-power)) ~ logelevated * factor(sample.size), data = dat_gof[1:24, ])
dat_gof$predicted_power <- (1.001 * exp(predict(fit)) + 0.01) /  (1 + exp(predict(fit)))

fit <- lm(log((type.I.error-0.05) / (0.995 - type.I.error)) ~ logobserved * factor(sample.size), 
          data = dat_gof[dat_gof$type.I.error> 0.05, ][1:24, ])
dat_gof$predicted_error <- (0.995 * exp(predict(fit, newdata = dat_gof)) + 0.05) /  (1 + exp(predict(fit, newdata = dat_gof)))


q1 <- ggplot(data = dat_gof[1:24, ], aes(x = logobserved, y = type.I.error, label = label)) + 
      geom_text(hjust = 1, vjust = -0.5, size = 3) +
      geom_point(aes(x = logobserved, y = type.I.error, 
                    group = factor(sample.size), 
                    shape = factor(sample.size), 
                    color = factor(sample.size)),
                 size = 6, alpha = 0.7) +
      theme_bw() + 
      geom_smooth(aes(x = logobserved, y = predicted_error,
                      group = factor(sample.size), 
                      linetype = factor(sample.size), 
                      color = factor(sample.size)),
                  method = "loess", se = F) + 
      geom_hline(yintercept = 0.05, linetype = "dotdash", color = "gray") +
      xlab("Log number of observed cases") +
      ylab("Type I error") + 
   scale_color_manual(name = "Sample size", values = my3cols) + 
   scale_shape_discrete(name = "Sample size") + 
   scale_linetype_discrete(guide = "none") +
   ggtitle("(a) Testing parametric assumption when model is correctly specified.")
      
q1


q2 <- ggplot(data = dat_gof[1:24, ], aes(x = logelevated, y = power, label = label)) + 
      geom_text(hjust = 1, vjust = -0.5, size = 3)+
      geom_point(aes(x = logelevated, y = power, 
                    group = factor(sample.size), 
                    shape = factor(sample.size), 
                    color = factor(sample.size)),
                 size = 6, alpha = 0.7) + 
      theme_bw() + 
      geom_smooth(aes(x = logelevated, y = predicted_power,
                    group = factor(sample.size), 
                    linetype = factor(sample.size),
                    color = factor(sample.size)),
                  method = "loess", se = F) + 
      xlab("Log number of harmed cases") +
      ylab("Power") + labs(shape = "Sample size") + labs(linetype = "Sample size") +
   scale_color_manual(name = "Sample size", values = my3cols) + 
   scale_shape_discrete(name = "Sample size") + 
   scale_linetype_discrete(guide = "none")+
   ggtitle("(b) Testing parametric assumption when model is misspecified.")
q2




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
q3


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
q4

# patchwork q1-4
library(patchwork)
(((q1+q2) / (q3+q4)) & theme(legend.position = "bottom")) + 
   plot_layout(guides = "collect")

# used for summarizing latest simulation results 05192021
# summarize MLE simulation results
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
            setwd(paste0("/users/yzhu/parmix/sim_gof_05192021_s", scenario, "_n", sample_size))
            filenames <- list.files(path = getwd())
            if ("sim_mle_s1to8_n1500to15k_05192021.csv" %in% filenames) {
                filenames <- setdiff(filenames, "sim_mle_s1to8_n1500to15k_05192021.csv")
            }

            # remove the error and out files
            if (grepl("error", filenames)) {
                filenames <- filenames[-which(grepl("error", filenames))]
            }
            if (grepl("out", filenames)) {
                filenames <- filenames[-which(grepl("out", filenames))]
            }

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
setwd("/users/yzhu/parmix")
write.csv(out, file = "sim_mle_s1to8_n1500to15k_05192021.csv")