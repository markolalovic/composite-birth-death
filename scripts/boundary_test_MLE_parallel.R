# scripts/boundary_test_MLE_parallel.R
#
# Compute the boundary-test statistics used in Figure 8.
#
# Test:
#   H0: b2 = 0
#   H1: b2 > 0
#
# Input:
#   data/simulated/stats_null_T1000_M1000.rds
#
# Output:
#   data/boundary_test_MLE.csv
#
# The input file contains sufficient statistics from 1000 trajectories of the
# null complete-graph SIS model conditioned to survive up to T = 1000
#
# For each replicate script computes the unconstrained conditional MLE
# and the Fisher-information-based standard error of b2
#
# The standardized statistic is
#   Z_b2 = b2_hat / se_b2
#
# The projected one-sided Wald statistic is also saved:
#   W_iT = max(0, Z_b2)^2
#

library(parallel)

source("R/01_rates.R")
source("R/generator.R")
source("R/pf.R")
source("R/mle.R")
source("R/fisher.R")


obj <- readRDS("data/simulated/stats_null_T1000_M1000.rds")

N <- obj$N
M <- obj$M
stats_list <- obj$stats

alpha <- 0.05
z_crit <- qnorm(1 - alpha)

init <- c(b1 = 2.5, b2 = 0.2, mu = 0.9)

ncores <- detectCores()
cat("Using", ncores, "cores for M =", M, "replicates\n")


run_one <- function(m) {
  stats <- stats_list[[m]]

  # unconstrained conditional MLE
  fit <- fit_mle_complete(stats, N, init = init)

  V_obj <- vcov_mle_complete(fit$par_hat, stats, N)

  if (is.null(V_obj)) {
    return(c(
      b1_hat = unname(fit$par_hat["b1"]),
      b2_hat = unname(fit$par_hat["b2"]),
      mu_hat = unname(fit$par_hat["mu"]),
      se_b2 = NA_real_,
      Z_b2 = NA_real_,
      W_iT = NA_real_,
      reject = NA_real_,
      convergence = fit$convergence
    ))
  }

  Sigma_hat <- V_obj$vcov

  se_b2 <- sqrt(Sigma_hat["b2", "b2"])
  Z_b2 <- fit$par_hat["b2"] / se_b2

  W_iT <- max(0, Z_b2)^2
  reject <- Z_b2 > z_crit

  c(
    b1_hat = unname(fit$par_hat["b1"]),
    b2_hat = unname(fit$par_hat["b2"]),
    mu_hat = unname(fit$par_hat["mu"]),
    se_b2 = unname(se_b2),
    Z_b2 = unname(Z_b2),
    W_iT = unname(W_iT),
    reject = as.numeric(reject),
    convergence = fit$convergence
  )
}


cat("Starting parallel computation...\n")
t0 <- proc.time()

results_list <- mclapply(1:M, run_one, mc.cores = ncores)

elapsed <- (proc.time() - t0)["elapsed"]
cat("Done in", round(elapsed, 1), "seconds\n")


results_mat <- do.call(rbind, results_list)

results <- data.frame(
  replicate = 1:M,
  b1_hat = results_mat[, "b1_hat"],
  b2_hat = results_mat[, "b2_hat"],
  mu_hat = results_mat[, "mu_hat"],
  se_b2 = results_mat[, "se_b2"],
  Z_b2 = results_mat[, "Z_b2"],
  W_iT = results_mat[, "W_iT"],
  reject = as.logical(results_mat[, "reject"]),
  convergence = as.integer(results_mat[, "convergence"]),
  n_tries = obj$n_tries
)

write.csv(
  results,
  file = "data/boundary_test_MLE.csv",
  row.names = FALSE
)

emp_reject_rate <- mean(results$reject, na.rm = TRUE)

cat("\nTrue parameters under H0:\n")
print(obj$theta)

cat("\nMean of Z_b2:", mean(results$Z_b2, na.rm = TRUE), "\n")
cat("SD of Z_b2:  ", sd(results$Z_b2, na.rm = TRUE), "\n")
cat("Proportion with Z_b2 < 0:", mean(results$Z_b2 < 0, na.rm = TRUE), "\n")
cat("Empirical rejection rate at alpha =", alpha, ":", emp_reject_rate, "\n")
cat("Optim convergence code 0:", sum(results$convergence == 0), "out of", M, "\n")
cat("Average tries for a survivor:", mean(results$n_tries), "\n")
