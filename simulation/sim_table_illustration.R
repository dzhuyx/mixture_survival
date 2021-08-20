library(ggplot2)
my2cols <- c("#00BFC4", "#F8777d")

# scenario parameters for illustrating figures
scenario_par0 <- rbind(c(10, 10, 2), c(10, 20, 2),
                       c(1, 10, 2), c(1, 20, 2),
                       c(10, 10, 1), c(10, 20, 1),
                       c(1, 10, 1), c(1, 20, 1))
expit <- function(x) {
      return(exp(x) / (1 + exp(x)))
}

myexprate <- function(theta, t) {
      out <- theta[1] + theta[2] * expit(theta[3])^t
      return(out)
}

risk_mat <- do.call(rbind, lapply(1:8, function(i) {
      x <- scenario_par0[i, ]
      temp <- sapply(0:100, function(t) {
            return(myexprate(x, t))
      })
      short <- max(which(temp >= min(temp) + 0.1))
      col <- rep(NA, 101)
      # for slow and fast decay 
      if (i %in% 1:4) {
            col[1:short] <- rep(my2cols[1], short)
      } else {
            col[1:short] <- rep(my2cols[2], short)
      }
      # high and low baseline risk
      if (i %in% c(1, 2, 5, 6)) {
            col[short:100+1] <- rep(my2cols[2], 100-short+1)
      } else {
            col[short:100+1] <- rep(my2cols[1], 100-short+1)
      }
      
      return(data.frame(time = 0:100, 
                        scenario = rep(i, 101),
                        risk = temp,
                        col = col))
}))

# y: 0-30, x: 0:100
for (scenario in 1:8) {
      q <- ggplot(data = risk_mat[risk_mat$scenario == scenario, ]) + 
            geom_line(aes(x = time, y = risk, col = col), size = 1) + 
            geom_hline(yintercept = max(risk_mat$risk[risk_mat$scenario == scenario]),
                       col = ifelse(scenario %in% c(1, 3, 5, 7),
                                    my2cols[1], my2cols[2]),
                       linetype = "dashed", size = 1) +
            theme(text = element_text(size=20),
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.position = "none",
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            xlab("") + ylab("") +
            scale_color_manual(values = unique(risk_mat$col[risk_mat$scenario == scenario])) + 
            xlim(c(0, 80)) + ylim(c(0, 30))
      
      jpeg(filename = paste0("table_s", scenario, ".jpeg"),
           width = 120, height = 100)
      plot(q)
      dev.off()
}

      
