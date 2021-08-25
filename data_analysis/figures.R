# figures for data analysis

# figure for harm risk factors
rm(list = ls())
library(ggplot2)
library(scales)
library(patchwork)

my2cols <- c("#E7B800", "#2E9FDF")
my3cols <- c(my2cols, "#FC4E07")
my4cols <- c(my3cols, "#66cd00")
my5cols <- c(my4cols, "#AC92EB")
my6cols <- c(my5cols, "#AFC1CC")
load("synthetic_harmriskfactor.rda")
# extract ah data
harm_ah_riskfactor <- list(age = harm_ah_age, 
                           gender = harm_ah_gender,
                           HT = harm_ah_HT,
                           DM = harm_ah_DM,
                           AF = harm_ah_AF,
                           stroke = harm_ah_stroke,
                           dd = harm_ah_dd,
                           SES = harm_ah_SES,
                           area = harm_ah_area,
                           teaching = harm_ah_teaching,
                           gs2 = harm_ah_gs2)
varlist_harmriskfactor <- c("age", "gender", "HT", "DM", "AF",
                            "stroke", "dd", "SES", "area", 
                            "teaching", "gs2")
ah <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(harm_ah_riskfactor[[var]]$est)
}))
se <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(sqrt(diag(harm_ah_riskfactor[[var]]$sigma2)))
}))
ci_lower <- ah - qnorm(0.975) * se
ci_upper <- ah + qnorm(0.975) * se

# extract rah data
harm_rah_riskfactor <- list(age = harm_rah_age,
                            gender = harm_rah_gender,
                            HT = harm_rah_HT,
                            DM = harm_rah_DM,
                            AF = harm_rah_AF,
                            stroke = harm_rah_stroke,
                            dd = harm_rah_dd,
                            SES = harm_rah_SES, 
                            area = harm_rah_area,
                            teaching = harm_rah_teaching,
                            gs2 = harm_rah_gs2)
lograh <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(harm_rah_riskfactor[[var]]$est)
}))
logse <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(sqrt(diag(harm_rah_riskfactor[[var]]$sigma2)))
}))
rah <- exp(lograh)
rah_ci_lower <- exp(lograh - qnorm(0.975) * logse)
rah_ci_upper <- exp(lograh + qnorm(0.975) * logse)

dat_plot <- data.frame(ah = ah, ah_se = se,
                       ah_ci_lower = ci_lower, ah_ci_upper = ci_upper,
                       rah = rah, lograh_se = logse,
                       rah_ci_lower = rah_ci_lower, rah_ci_upper = rah_ci_upper,
                       var = c(rep("Age", 3),
                               rep("Gender male", 2),
                               rep("Hypertension", 2),
                               rep("Diabetes mellitus", 2),
                               rep("Atrial fibrillation", 2),
                               rep("Stroke history", 2), 
                               rep("Hospitalization", 2),
                               rep("Socioeconomic status", 3),
                               rep("Area", 5),
                               rep("Hospital type", 5),
                               rep("Generalist hospital", 2)),
                       val = c(1:3, rep(1:2, 6), 1:3,
                               1:5, 1:5, 1:2))
dat_plot$var <- factor(dat_plot$var, 
                       levels = c("Age", "Gender male", "Hypertension",
                                  "Diabetes mellitus", "Atrial fibrillation",
                                  "Stroke history", "Hospitalization",
                                  "Socioeconomic status", "Area",
                                  "Hospital type", "Generalist hospital"))

# perform inter-group heterogeneity hypothesis testing
pval_ah <- NULL
pval_rah <- NULL
for (i in 1:length(harm_ah_riskfactor)) {
      temp_ah <- harm_ah_riskfactor[[i]]
      temp_rah <- harm_rah_riskfactor[[i]]
      
      if (length(temp_ah$est) == 2) {
            # pval for attributable harm
            tstat <- diff(temp_ah$est)
            tstat_se <- t(c(1, -1)) %*% temp_ah$sigma2 %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah <- c(pval_ah, pval)
            # pval for relative attributable harm
            tstat <- diff(temp_rah$est)
            tstat_se <- t(c(1, -1)) %*% temp_rah$sigma2 %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah <- c(pval_rah, pval)
      } else {
            p <- length(temp_ah$est)
            tempmat <- cbind(0, diag(1, p-1)) - 
                  cbind(diag(1, p-1), 0)
            # pval for attributable harm
            tempdiff <- tempmat %*% temp_ah$est
            tempcov <- tempmat %*% temp_ah$sigma2 %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah <- c(pval_ah, pval)
            # pval for relative attributable harm
            tempdiff <- tempmat %*% temp_rah$est
            tempcov <- tempmat %*% temp_rah$sigma2 %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah <- c(pval_rah, pval)
      }
}

group_names <- levels(dat_plot$var)
group_names_ah <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_ah[i], ")"))
})
group_names_rah <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_rah[i], ")"))
})

# create figures
# figures for thirty-day harm
q1_ah <- ggplot(data = dat_plot,
                aes(x = var, y = ah * 10000, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = ah_ci_lower * 10000,
                        ymax = ah_ci_upper * 10000,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_ah) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 120), oob = rescale_none) +
      xlab("") + ylab("Harmed cases per 10,000 patients")

q1_rah <- ggplot(data = dat_plot,
                aes(x = var, y = rah, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = rah_ci_lower,
                        ymax = rah_ci_upper,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_rah) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 45), oob = rescale_none) +
      xlab("") + ylab("Relative attributable harm")

((((q1_ah) / (q1_rah)) & theme(legend.position = "bottom")) + 
            plot_layout(guides = "collect")) # 1600 * 500

# figures for time-dynamic harm in different age groups
dat_plot <- rbind(data.frame(time = rep(harm_curve_ah_age$t, 3),
                             ah = c(harm_curve_ah_age$ah,
                                    harm_curve_ah_age$lower,
                                    harm_curve_ah_age$upper),
                             rah = log(c(harm_curve_rah_age$rah,
                                         harm_curve_rah_age$lower,
                                         harm_curve_rah_age$upper)),
                             type = c(rep(c("est", "lower", "upper"), 
                                          each = 365*3)),
                             age = rep(harm_curve_ah_age$age, 3)))
dat_plot$type <- factor(dat_plot$type)
dat_plot$age <- factor(dat_plot$age)
colorname <- "Age"
labelnames <- c("<= 40", "40 - 60", "> 60")
q1_ahcurve_age <- ggplot(data = dat_plot) + 
      geom_line(aes(x = time, y = ah * 10000,
                    group = interaction(type, age),
                    color = age, linetype = type), size = 1) + 
      theme_bw() + 
      scale_color_manual(name = colorname, labels = labelnames, values = my3cols) + 
      scale_linetype_discrete(guide = "none") + 
      xlab("Time (days)") + ylab("Harmed cases per 10,000 patients") + 
      xlim(1, 365) + ylim(0, 120) 
q1_ahcurve_age

q1_rahcurve_age <- ggplot(data = dat_plot) + 
      geom_line(aes(x = time, y = rah, 
                    group = interaction(age, type),
                    color = age, linetype = type), size = 1) + 
      theme_bw() + 
      scale_color_manual(name = colorname, labels = labelnames, values = my3cols) + 
      scale_linetype_discrete(guide = "none") + 
      xlab("Time (days)") + ylab("Relative attribuable harm (log)") + 
      xlim(1, 365) 
q1_rahcurve_age

((((q1_ahcurve_age) | (q1_rahcurve_age)) & theme(legend.position = "bottom")) + 
            plot_layout(guides = "collect")) # 1600 * 500



# figure for harm-pivotal time points 
rm(list = ls())
my2cols <- c("#E7B800", "#2E9FDF")
my3cols <- c(my2cols, "#FC4E07")
my4cols <- c(my3cols, "#66cd00")
my5cols <- c(my4cols, "#AC92EB")
my6cols <- c(my5cols, "#AFC1CC")
load("synthetic_landmarktime.rda")
tau_list_ah <- list(age = fit_ah_age,
                    gender = fit_ah_gender,
                    HT = fit_ah_HT,
                    DM = fit_ah_DM,
                    AF = fit_ah_AF,
                    stroke = fit_ah_stroke,
                    dd = fit_ah_dd,
                    SES = fit_ah_SES,
                    area = fit_ah_area,
                    teaching = fit_ah_teaching,
                    gs2 = fit_ah_gs2)
# fit_ah_x results are for 0.8 and 0.5 of one-year harm
varlist_harmriskfactor <- c("age", "gender", "HT", "DM", "AF",
                            "stroke", "dd", "SES", "area", 
                            "teaching", "gs2")
# extract ah information
tau_ah <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(tau_list_ah[[var]]$est)
}))
tause_ah <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(sqrt(diag(tau_list_ah[[var]]$sigma)))
}))
ah_ci_lower <- tau_ah - qnorm(0.975) * tause_ah
ah_ci_upper <- tau_ah + qnorm(0.975) * tause_ah
ah_ci_lower <- sapply(1:length(ah_ci_lower), function(i) {
      return(max(ah_ci_lower[i], 0))
})
# extract rah information
tau_list_rah <- list(age = fit_rah_age,
                    gender = fit_rah_gender,
                    HT = fit_rah_HT,
                    DM = fit_rah_DM,
                    AF = fit_rah_AF,
                    stroke = fit_rah_stroke,
                    dd = fit_rah_dd,
                    SES = fit_rah_SES,
                    area = fit_rah_area,
                    teaching = fit_rah_teaching,
                    gs2 = fit_rah_gs2)
tau_rah <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(tau_list_rah[[var]]$est)
}))
tause_rah <- do.call("c", lapply(varlist_harmriskfactor, function(var) {
      return(sqrt(diag(tau_list_rah[[var]]$sigma)))
}))
rah_ci_lower <- tau_rah - qnorm(0.975) * tause_rah
rah_ci_upper <- tau_rah + qnorm(0.975) * tause_rah
rah_ci_lower <- sapply(1:length(rah_ci_lower), function(i) {
      return(max(rah_ci_lower[i], 0))
})
# dat_plot_all
dat_plot_all <- data.frame(tau_ah = tau_ah, tause_ah = tause_ah,
                           ah_ci_lower = ah_ci_lower,
                           ah_ci_upper = ah_ci_upper,
                           tau_rah = tau_rah, tause_rah = tause_rah,
                           rah_ci_lower = rah_ci_lower,
                           rah_ci_upper = rah_ci_upper,
                           var = c(rep(rep("Age", 3), 2),
                                   rep(rep("Gender male", 2), 2),
                                   rep(rep("Hypertension", 2), 2),
                                   rep(rep("Diabetes mellitus", 2), 2),
                                   rep(rep("Atrial fibrillation", 2), 2),
                                   rep(rep("Stroke history", 2), 2), 
                                   rep(rep("Hospitalization", 2), 2),
                                   rep(rep("Socioeconomic status", 3), 2),
                                   rep(rep("Area", 5), 2),
                                   rep(rep("Hospital type", 5), 2),
                                   rep(rep("Generalist hospital", 2), 2)),
                           val = c(rep(1:3, 2), rep(rep(1:2, 6), 2), 
                                   rep(1:3, 2), rep(1:5, 2), 
                                   rep(1:5, 2), rep(1:2, 2)),
                           harm = c(rep(c(0.8, 0.5), each = 3),
                                    rep(rep(c(0.8, 0.5), each = 2), 6),
                                    rep(c(0.8, 0.5), each = 3),
                                    rep(c(0.8, 0.5), each = 5),
                                    rep(c(0.8, 0.5), each = 5),
                                    rep(c(0.8, 0.5), each = 2)))
dat_plot_all$var <- factor(dat_plot_all$var, 
                       levels = c("Age", "Gender male", "Hypertension",
                                  "Diabetes mellitus", "Atrial fibrillation",
                                  "Stroke history", "Hospitalization",
                                  "Socioeconomic status", "Area",
                                  "Hospital type", "Generalist hospital"))
# hypothesis testing, in results 80 then 50
pval_ah50 <- NULL
pval_ah80 <- NULL
pval_rah50 <- NULL
pval_rah80 <- NULL
for (i in 1:length(tau_list_ah)) {
      temp_ah <- tau_list_ah[[i]]
      temp_rah <- tau_list_rah[[i]]
      
      p <- length(temp_ah$est) / 2
      if (p == 2) {
            # pval for ah80
            tstat <- diff(temp_ah$est[1:p])
            tstat_se <- t(c(1, -1)) %*% temp_ah$sigma[1:p, 1:p] %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah80 <- c(pval_ah80, pval)
            
            # pval for ah50
            tstat <- diff(temp_ah$est[1:p+p])
            tstat_se <- t(c(1, -1)) %*% temp_ah$sigma[1:p+p, 1:p+p] %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah50 <- c(pval_ah50, pval)
            
            # pval for rah80
            tstat <- diff(temp_rah$est[1:p])
            tstat_se <- t(c(1, -1)) %*% temp_rah$sigma[1:p, 1:p] %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah80 <- c(pval_rah80, pval)
            
            # pval for rah50
            tstat <- diff(temp_rah$est[1:p+p])
            tstat_se <- t(c(1, -1)) %*% temp_rah$sigma[1:p+p, 1:p+p] %*% c(1, -1)
            pval <- 1 - pnorm(abs(tstat), 0, tstat_se)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah50 <- c(pval_rah50, pval)
      } else {
            tempmat <- cbind(0, diag(1, p-1)) - 
                  cbind(diag(1, p-1), 0)
            
            # pval for ah80
            tempdiff <- tempmat %*% temp_ah$est[1:p]
            tempcov <- tempmat %*% temp_ah$sigma[1:p, 1:p] %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah80 <- c(pval_ah80, pval)
            
            # pval for ah50
            tempdiff <- tempmat %*% temp_ah$est[1:p+p]
            tempcov <- tempmat %*% temp_ah$sigma[1:p+p, 1:p+p] %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_ah50 <- c(pval_ah50, pval)
            
            # pval for rah80
            tempdiff <- tempmat %*% temp_rah$est[1:p]
            tempcov <- tempmat %*% temp_rah$sigma[1:p, 1:p] %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah80 <- c(pval_rah80, pval)
            
            # pval for rah50
            tempdiff <- tempmat %*% temp_rah$est[1:p+p]
            tempcov <- tempmat %*% temp_rah$sigma[1:p+p, 1:p+p] %*% t(tempmat)
            tstat <- t(tempdiff) %*% solve(tempcov) %*% tempdiff
            pval <- 1-pchisq(tstat, df = p-1)
            pval <- ifelse(pval < 0.01, "< 0.01", round(pval, 2))
            pval_rah50 <- c(pval_rah50, pval)
      }
}

group_names <- levels(dat_plot_all$var)
group_names_ah50 <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_ah50[i], ")"))
})
group_names_ah80 <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_ah80[i], ")"))
})
group_names_rah50 <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_rah50[i], ")"))
})
group_names_rah80 <- sapply(1:length(group_names), function(i) {
      return(paste0(group_names[i], "\n(", pval_rah80[i], ")"))
})
dat_plot50 <- dat_plot_all[dat_plot_all$harm == 0.5, ]
dat_plot80 <- dat_plot_all[dat_plot_all$harm == 0.8, ]
# make figures
# attributable harm
q2_ah50 <- ggplot(data = dat_plot50,
                  aes(x = var, y = tau_ah, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = ah_ci_lower,
                        ymax = ah_ci_upper,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_ah50) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 100), oob = rescale_none) +
      xlab("") + ylab("Number of days to\n50% of one-year attributable harm")

q2_ah80 <- ggplot(data = dat_plot80,
                  aes(x = var, y = tau_ah, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = ah_ci_lower,
                        ymax = ah_ci_upper,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_ah80) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 300), oob = rescale_none) +
      xlab("") + ylab("Number of days to\n80% of one-year attributable harm")

((((q2_ah50) / (q2_ah80)) & theme(legend.position = "bottom")) + 
            plot_layout(guides = "collect")) # 1600 * 500

# relative attributable harm
q2_rah50 <- ggplot(data = dat_plot50,
                  aes(x = var, y = tau_rah, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = rah_ci_lower,
                        ymax = rah_ci_upper,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_rah50) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 2000), oob = rescale_none) +
      xlab("") + ylab("Number of days to\nan relative attributable harm of 1")

q2_rah80 <- ggplot(data = dat_plot80,
                  aes(x = var, y = tau_rah, fill = factor(val))) +
      geom_bar(position = position_dodge(), stat = "identity",
               colour = "black", size = .3) +
      geom_errorbar(aes(ymin = rah_ci_lower,
                        ymax = rah_ci_upper,
                        group = val),
                    width = 0.3, position = position_dodge(0.9)) +
      scale_fill_manual(name = "Group", values = my5cols, 
                        labels = c("<=40 (age) /\nLow (Socioeconomic status) /\nNorthern (area) /\nUnclassified (hospital type) /\nNo (others)", 
                                   "40-60 (age) /\nMedium (Socioeconomic status) /\nMiddle (area) /\nRegional (hospital type) /\nYes (others)", 
                                   "> 60 (age) /\nHigh (Socioeconomic status) /\nSouthern (area) /\nLocal (hospital type)",
                                   "Gaoping (area) /\nPrivate practice (hospital type)",
                                   "Taipei or Eastern (area) /\nMedical center (hospital type)")) +
      scale_x_discrete(labels = group_names_rah80) +
      theme_bw() + theme(axis.text.x = element_text(angle = -20, hjust = 0.5)) + 
      scale_y_continuous(limits = c(0, 300), oob = rescale_none) +
      xlab("") + ylab("Number of days to\nan relative attributable harm of 4")


(((q2_ah50) / (q2_ah80) / 
        (q2_rah80) / (q2_rah50) & theme(legend.position = "bottom")) + 
            plot_layout(guides = "collect")) # 1600 * 500


# figures for profile analysis
rm(list = ls())
load("taiwan_parmix_05192021_gsprofile.rda")

# gs1 observed vs standard
dat_plot <- harm_curve_ah_profile[harm_curve_ah_profile$type %in% 
                                        c(), ]



# coefficient estimates
rm(list = ls())
load("taiwan_parmix_05192021_poly_start2189.rda")
varlist_long <- c("gs2", "gender", "HT", "DM", "AF", "before_STROKE", "dd",                  
                  "SocioeconomicStatus2", "SocioeconomicStatus1", 
                  "area4", "area2", "area3", "area5",
                  "Teaching2", "age_young", "age_old")
varlist_short <- c("gs2", "gender", "HT", "DM", "before_STROKE",
                   "area4", "Teaching4", "Teaching3",
                   "Teaching0", "age_young", "age_old")
varlist_shortint <- c("gs2_age_old")
varlist_all <- c("intercept", unique(c(varlist_long, varlist_short, varlist_shortint)))

est <- fit$par
se <- sqrt(diag(solve(fit$hessian)))
lower <- est - qnorm(0.975) * se
upper <- est + qnorm(0.975) * se


beta_est <- data.frame(var = c("intercept", varlist_long),
                        beta = est[1:(length(varlist_long) + 1)],
                       beta_lower = lower[1:(length(varlist_long) + 1)],
                       beta_upper = upper[1:(length(varlist_long) + 1)])
gamma1_est <- data.frame(var = c("intercept", varlist_short, varlist_shortint),
                         gamma1 = est[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                             length(varlist_long) + 1],
                         gamma1_lower = lower[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                                  length(varlist_long) + 1],
                         gamma1_upper = upper[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                                  length(varlist_long) + 1])
gamma2_est <- data.frame(var = c("intercept", varlist_short, varlist_shortint),
                         gamma2 = est[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                             length(varlist_long) + 1 + length(varlist_short) + 
                                             length(varlist_shortint) + 1],
                         gamma2_lower = lower[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                                  length(varlist_long) + 1 + length(varlist_short) + 
                                                  length(varlist_shortint) + 1],
                         gamma2_upper = upper[1:(length(varlist_short) + length(varlist_shortint) + 1) +
                                                  length(varlist_long) + 1 + length(varlist_short) + 
                                                  length(varlist_shortint) + 1])
coef_est <- merge(beta_est, gamma1_est, by = "var", all = T)
coef_est <- merge(coef_est, gamma2_est, by = "var", all = T)
# order rows
coef_est$var <- factor(coef_est$var, 
                       levels = c("intercept", "age_young", "age_old",
                                  "gender", "DM", "HT", "AF", "before_STROKE", "dd",
                                  "SocioeconomicStatus1", "SocioeconomicStatus2",
                                  "area2", "area3", "area4", "area5",
                                  "Teaching0", "Teaching2", "Teaching3", "Teaching4",
                                  "gs2", "gs2_age_old"))
coef_est <- coef_est[order(coef_est$var), ]
result_est <- t(apply(coef_est, 1, function(x) {
      x <- unlist(x)
      x[-1] <- round(as.numeric(x[-1]), 2)
      beta <- paste0(x[2], " (", x[3], ", ", x[4], ")")
      gamma1 <- paste0(x[5], " (", x[6], ", ", x[7], ")")
      gamma2 <- paste0(x[8], " (", x[9], ", ", x[10], ")")
      return(c(beta, gamma1, gamma2))
}))
result_est <- data.frame(result_est)
write.csv(cbind(coef_est[, 1], result_est), 
          file = "taiwan_parmix_est.csv")

