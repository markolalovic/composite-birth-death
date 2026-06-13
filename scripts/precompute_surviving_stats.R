# scripts/precompute_surviving_stats.R
#
# Precompute sufficient statistics from surviving trajectories.
#
# Baseline parameter, used for Figures 4--7:
#   theta0 = c(b1 = 1.01, b2 = 3.70, mu = 1)
#
# Output for T = 100, 200, ..., 1000:
#   data/simulated/stats_baseline_T{T}_M200.rds
#
# Null parameter, used for Figure 8:
#   theta_null = c(b1 = 2.875, b2 = 0, mu = 1)
#
# Output:
#   data/simulated/stats_null_T1000_M1000.rds
#
# Each file contains:
#   stats    : list of sufficient-statistics objects
#   n_tries  : number of simulated paths required to obtain each survivor
#
# The baseline seed rule is seed = 1000 * T + m, matching the estimator
# scripts used for Figures 4--7.

library(parallel)

source("R/01_rates.R")
source("R/02_simulation.R")
source("R/03_path_statistics.R")

N <- 100
I0 <- 10

theta_baseline <- c(b1 = 1.01, b2 = 3.70, mu = 1.0)
theta_null     <- c(b1 = 2.875, b2 = 0.00, mu = 1.0)

ncores <- detectCores()

dir.create("data/simulated", recursive = TRUE, showWarnings = FALSE)


simulate_stats_one <- function(m, T_max, theta, seed) {
  set.seed(seed)

  tries <- 0

  repeat {
    tries <- tries + 1
    traj <- simulate_path(theta, N = N, I0 = I0, time_max = T_max)
    if (tail(traj$k, 1) > 0) break
  }

  list(
    stats = compute_sufficient_stats(traj, N),
    n_tries = tries
  )
}


simulate_stats_set <- function(theta, T_max, M, outfile, seed_fun) {
  cat("\n=== Simulating", outfile, "===\n")
  t0 <- proc.time()

  sims <- mclapply(
    1:M,
    function(m) simulate_stats_one(
      m = m,
      T_max = T_max,
      theta = theta,
      seed = seed_fun(T_max, m)
    ),
    mc.cores = ncores
  )

  obj <- list(
    N = N,
    I0 = I0,
    T_max = T_max,
    M = M,
    theta = theta,
    stats = lapply(sims, `[[`, "stats"),
    n_tries = sapply(sims, `[[`, "n_tries")
  )

  saveRDS(obj, outfile)

  elapsed <- (proc.time() - t0)["elapsed"]

  cat("Saved", outfile, "\n")
  cat("Done in", round(elapsed, 1), "seconds\n")
  cat("Avg tries:", mean(obj$n_tries), "\n")
}

# Baseline simulations for Figures 4--7.
for (T_max in seq(100, 1000, by = 100)) {
  simulate_stats_set(
    theta = theta_baseline,
    T_max = T_max,
    M = 200,
    outfile = sprintf("data/simulated/stats_baseline_T%d_M200.rds", T_max),
    seed_fun = function(T, m) 1000 * T + m
  )
}

# Null simulations for Figure 8.
simulate_stats_set(
  theta = theta_null,
  T_max = 1000,
  M = 1000,
  outfile = "data/simulated/stats_null_T1000_M1000.rds",
  seed_fun = function(T, m) 1000 * T + m
)

cat("\nAll done.\n")