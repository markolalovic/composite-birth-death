# scripts/check_qem_estimators.R
#
# Diagnostic comparison of the stable Q-EM algorithms in R/qem.R with
# the other estimator implementations
#
# Comparisons:
#   Unmarked births:
#     fit_qmle_complete()
#     fit_qem_unmarked_complete()
#
#   Marked births:
#     fit_mle_marked_complete()
#     fit_qem_marked_complete()
#
# Input:
#   data/simulated/stats_baseline_T1000_M200.rds
#

source("R/01_rates.R")
source("R/generator.R")
source("R/pf.R")
source("R/qmle.R")
source("R/mle.R")
source("R/mle_marked.R")
source("R/qem.R")


obj <- readRDS("data/simulated/stats_baseline_T1000_M200.rds")

N <- obj$N
stats_list <- obj$stats

init <- c(b1 = 1.0, b2 = 3.5, mu = 0.9)

# small test
replicates <- 1:20

run_one <- function(m) {
  stats <- stats_list[[m]]

  fit_qmle <- fit_qmle_complete(
    stats,
    N,
    init = init,
    verbose = FALSE
  )

  fit_qem_unmarked <- fit_qem_unmarked_complete(
    stats,
    N,
    init = init
  )

  fit_mle_marked <- fit_mle_marked_complete(
    stats,
    N,
    init = init
  )

  fit_qem_marked <- fit_qem_marked_complete(
    stats,
    N,
    init = init
  )

  c(
    replicate = m,

    b1_qmle = unname(fit_qmle$par_hat["b1"]),
    b2_qmle = unname(fit_qmle$par_hat["b2"]),
    mu_qmle = unname(fit_qmle$par_hat["mu"]),

    b1_qem_unmarked = unname(fit_qem_unmarked$par_hat["b1"]),
    b2_qem_unmarked = unname(fit_qem_unmarked$par_hat["b2"]),
    mu_qem_unmarked = unname(fit_qem_unmarked$par_hat["mu"]),

    b1_mle_marked = unname(fit_mle_marked$par_hat["b1"]),
    b2_mle_marked = unname(fit_mle_marked$par_hat["b2"]),
    mu_mle_marked = unname(fit_mle_marked$par_hat["mu"]),

    b1_qem_marked = unname(fit_qem_marked$par_hat["b1"]),
    b2_qem_marked = unname(fit_qem_marked$par_hat["b2"]),
    mu_qem_marked = unname(fit_qem_marked$par_hat["mu"]),

    qem_unmarked_converged = as.numeric(fit_qem_unmarked$converged),
    qem_unmarked_iterations = fit_qem_unmarked$iterations,

    qem_marked_converged = as.numeric(fit_qem_marked$converged),
    qem_marked_iterations = fit_qem_marked$iterations
  )
}


results <- as.data.frame(do.call(rbind, lapply(replicates, run_one)))

results$diff_b1_qmle_qem_unmarked <-
  results$b1_qmle - results$b1_qem_unmarked

results$diff_b2_qmle_qem_unmarked <-
  results$b2_qmle - results$b2_qem_unmarked

results$diff_mu_qmle_qem_unmarked <-
  results$mu_qmle - results$mu_qem_unmarked

results$diff_b1_mle_qem_marked <-
  results$b1_mle_marked - results$b1_qem_marked

results$diff_b2_mle_qem_marked <-
  results$b2_mle_marked - results$b2_qem_marked

results$diff_mu_mle_qem_marked <-
  results$mu_mle_marked - results$mu_qem_marked


cat("Unmarked: QMLE minus Q-EM summary\n")
print(summary(results[, c(
  "diff_b1_qmle_qem_unmarked",
  "diff_b2_qmle_qem_unmarked",
  "diff_mu_qmle_qem_unmarked"
)]))

cat("\nMarked: full conditional MLE minus marked Q-EM summary\n")
print(summary(results[, c(
  "diff_b1_mle_qem_marked",
  "diff_b2_mle_qem_marked",
  "diff_mu_mle_qem_marked"
)]))

cat("\nQ-EM convergence counts\n")
cat(
  "Unmarked:",
  sum(results$qem_unmarked_converged == 1),
  "out of",
  nrow(results),
  "\n"
)
cat(
  "Marked:",
  sum(results$qem_marked_converged == 1),
  "out of",
  nrow(results),
  "\n"
)
