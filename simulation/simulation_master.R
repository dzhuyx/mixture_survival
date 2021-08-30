# run simulation
# can also implement this on a computing cluster 
# using parallel jobs

for (scenario_index in as.character(1:8)) {
      for (n in c(1500, 5000, 15000)) {
            for (seed in 1:1000) {
                  source("scenario_a.R")
            }
      }
}

for (scenario_index in c("1b", "3b")) {
      for (n in c(1500, 5000, 15000)) {
            for (seed in 1:1000) {
                  source("scenario_b.R")
            }
      }
}

for (scenario_index in c("1c", "2c", "3c", "4c")) {
      for (n in c(1500, 5000, 15000)) {
            for (seed in 1:1000) {
                  source("scenario_c.R")
            }
      }
}

# all simulation results are stored in subfolder sim_gof_s$scenario_index$_n$n$
# under name sim_gof_s$scenario_index$_n$n$_seed$seed$.rda

# summarize simulation results, create Figure 1 and Table 2
source("sim_summary.R")

# create illustrative figures in Table 2
source("Table2_figures.R")
