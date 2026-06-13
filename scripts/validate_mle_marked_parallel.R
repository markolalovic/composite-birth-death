# scripts/validate_mle_marked_parallel.R
#
# UPDATED
# 
# Compute conditional MLEs for marked births from precomputed sufficient
# statistics at T = 200, 500, 1000.
#
# Input for T = 200, 500, 1000:
#   data/simulated/stats_baseline_T{T}_M200.rds
#
# Output for T = 200, 500, 1000:
#   data/estimates_mle_marked_{T}.csv
#
# These files are used for Figure 6.
#
# The precomputed statistics contain U_k and V_k separately, so the same saved
# data can be used for marked conditional MLEs, unmarked conditional MLEs, and
# unconditional MLEs.

library(parallel)

source("R/01_rates.R")
source("R/generator.R")
source("R/pf.R")
source("R/mle.R")
source("R/mle_marked.R")


T_grid <- c(200, 500, 1000)

init <- c(b1 = 1.0, b2 = 3.5, mu = 0.9)

ncores <- detectCores()
cat("Using", ncores, "cores\n")


run_one <- function(m, stats_list, N) {
  fit <- fit_mle_marked_complete(stats_list[[m]], N, init = init)

  c(
    b1_mle = unname(fit$par_hat["b1"]),
    b2_mle = unname(fit$par_hat["b2"]),
    mu_mle = unname(fit$par_hat["mu"])
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
    b1_mle = results_mat[, "b1_mle"],
    b2_mle = results_mat[, "b2_mle"],
    mu_mle = results_mat[, "mu_mle"],
    n_tries = obj$n_tries
  )

  n_failed <- sum(is.na(results_df$b1_mle))

  if (n_failed > 0) {
    cat("WARNING:", n_failed, "replicates failed (NA)\n")
  }

  outfile <- sprintf("data/estimates_mle_marked_%d.csv", T_max)
  write.csv(results_df, file = outfile, row.names = FALSE)

  cat("Saved", outfile, "\n")
  cat("True:", obj$theta, "\n")
  cat(
    "Cond. marked MLE mean:",
    colMeans(results_df[, c("b1_mle", "b2_mle", "mu_mle")], na.rm = TRUE),
    "\n"
  )
  cat("Avg tries:", mean(results_df$n_tries, na.rm = TRUE), "\n")
}

cat("\nAll done.\n")
