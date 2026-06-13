# scripts/compare_estimators_parallel.R
#
# Compute M = 200 estimates for T = 100, 200, ..., 1000:
#   - naive unconditional MLE
#   - conditional MLE
#   - QMLE
#
# Input for T = 100, 200, ..., 1000:
#   data/simulated/stats_baseline_T{T}_M200.rds
#
# Output for T = 100, 200, ..., 1000:
#   data/comparison_estimators_{T}.csv
#
# The file data/comparison_estimators_1000.csv is used for the unmarked panel
# of Figure 4. 
#
# The files for the whole T-grid are used for Figure 7.
#
# The columns b1_mle and b2_mle for T = 200, 500, 1000 are used for Figure 5.
#
# The precomputed statistics contain U_k, V_k, D_k, T_k and B_k = U_k + V_k,
# so the same saved data can be used for unconditional MLE, conditional MLE,
# QMLE, and marked-birth estimators.
#

library(parallel)

source("R/01_rates.R")
source("R/generator.R")
source("R/pf.R")
source("R/qmle.R")
source("R/mle.R")
source("R/unc_mle.R")


T_grid <- seq(100, 1000, by = 100)

init <- c(b1 = 1.0, b2 = 3.5, mu = 0.9)

ncores <- detectCores()
cat("Using", ncores, "cores\n")


run_one <- function(m, stats_list, N) {
  stats <- stats_list[[m]]

  fit_unc  <- fit_unc_mle_complete(stats, N, init = init)
  fit_qmle <- fit_qmle_complete(stats, N, init = init, verbose = FALSE)
  fit_mle  <- fit_mle_complete(stats, N, init = fit_qmle$par_hat)

  c(
    b1_unc  = unname(fit_unc$par_hat["b1"]),
    b2_unc  = unname(fit_unc$par_hat["b2"]),
    mu_unc  = unname(fit_unc$par_hat["mu"]),

    b1_mle  = unname(fit_mle$par_hat["b1"]),
    b2_mle  = unname(fit_mle$par_hat["b2"]),
    mu_mle  = unname(fit_mle$par_hat["mu"]),

    b1_qmle = unname(fit_qmle$par_hat["b1"]),
    b2_qmle = unname(fit_qmle$par_hat["b2"]),
    mu_qmle = unname(fit_qmle$par_hat["mu"])
  )
}


for (T_max in T_grid) {
  cat("\n=== T =", T_max, "===\n")
  t0 <- proc.time()

  infile <- sprintf("data/simulated/stats_baseline_T%d_M200.rds", T_max)
  obj <- readRDS(infile)

  N <- obj$N
  M <- obj$M
  stats_list <- obj$stats

  results_list <- mclapply(
    1:M,
    function(m) run_one(m, stats_list, N),
    mc.cores = ncores
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  cat("Done in", round(elapsed, 1), "seconds\n")

  results_mat <- do.call(rbind, results_list)

  results_df <- data.frame(
    replicate = 1:M,

    b1_unc  = results_mat[, "b1_unc"],
    b2_unc  = results_mat[, "b2_unc"],
    mu_unc  = results_mat[, "mu_unc"],

    b1_mle  = results_mat[, "b1_mle"],
    b2_mle  = results_mat[, "b2_mle"],
    mu_mle  = results_mat[, "mu_mle"],

    b1_qmle = results_mat[, "b1_qmle"],
    b2_qmle = results_mat[, "b2_qmle"],
    mu_qmle = results_mat[, "mu_qmle"],

    n_tries = obj$n_tries
  )

  outfile <- sprintf("data/comparison_estimators_%d.csv", T_max)
  write.csv(results_df, file = outfile, row.names = FALSE)

  cat("Saved", outfile, "\n")
  cat("True:", obj$theta, "\n")
  cat(
    "Uncond. MLE mean:",
    colMeans(results_df[, c("b1_unc", "b2_unc", "mu_unc")], na.rm = TRUE),
    "\n"
  )
  cat(
    "Cond. MLE mean:  ",
    colMeans(results_df[, c("b1_mle", "b2_mle", "mu_mle")], na.rm = TRUE),
    "\n"
  )
  cat(
    "QMLE mean:       ",
    colMeans(results_df[, c("b1_qmle", "b2_qmle", "mu_qmle")], na.rm = TRUE),
    "\n"
  )
  cat("Avg tries:", mean(results_df$n_tries, na.rm = TRUE), "\n")
}

cat("\nAll done.\n")
