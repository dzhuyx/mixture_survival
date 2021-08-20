library(gdata)

# reformat simulation results in `sim_mle_s1to8_n1500to15k.csv`
# to one appropriate manuscript tables
sim_result <- read.csv("sim_mle_s1to8_n1500to15k_05192021.csv")
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


colnames(dat_out) 
# arrange rows and output
out <- cbind(dat_out[dat_out$sample_size == 1500, c(1, 9, 11, 3:7)],
             dat_out[dat_out$sample_size == 5000, c(4:7)],
             dat_out[dat_out$sample_size == 15000, c(4:7)])

out

write.csv(out, file = "parmix_simulation_mle_table_05192021.csv")
