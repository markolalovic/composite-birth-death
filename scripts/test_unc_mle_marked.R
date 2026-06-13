# scripts/test_unc_mle_marked.R
#
# UPDATED
# 
# Compute naive unconditional MLEs for marked births at T = 1000.
#
# Input:
#   data/simulated/stats_baseline_T1000_M200.rds
#
# Output:
#   data/comparison_estimators_marked_1000.csv
#
# The input file contains sufficient statistics from 200 trajectories of the
# baseline complete-graph SIS model conditioned to survive up to T = 1000.
# The same precomputed statistics can be used for marked and unmarked
# estimators, since compute_sufficient_stats() stores both U_k, V_k and
# B_k = U_k + V_k.

source("R/01_rates.R")
source("R/unc_mle_marked.R")

obj <- readRDS("data/simulated/stats_baseline_T1000_M200.rds")

N <- obj$N
M <- obj$M
stats_list <- obj$stats

results <- lapply(1:M, function(m) {
  fit <- fit_unc_mle_marked_complete(stats_list[[m]], N)

  c(
    replicate = m,
    b1_unc = unname(fit$par_hat["b1"]),
    b2_unc = unname(fit$par_hat["b2"]),
    mu_unc = unname(fit$par_hat["mu"]),
    n_tries = obj$n_tries[m]
  )
})

results_mat <- do.call(rbind, results)

results_df <- data.frame(
  replicate = results_mat[, "replicate"],
  b1_unc    = results_mat[, "b1_unc"],
  b2_unc    = results_mat[, "b2_unc"],
  mu_unc    = results_mat[, "mu_unc"],
  n_tries   = results_mat[, "n_tries"]
)

write.csv(
  results_df,
  file = "data/comparison_estimators_marked_1000.csv",
  row.names = FALSE
)

cat("Saved data/comparison_estimators_marked_1000.csv\n")
cat("Sample means:\n")
print(colMeans(results_df[, c("b1_unc", "b2_unc", "mu_unc")], na.rm = TRUE))
cat("Average tries for a survivor:", mean(results_df$n_tries, na.rm = TRUE), "\n")
