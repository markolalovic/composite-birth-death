# R/fisher.R
#
# Fisher information and covariance matrix for the full conditional MLE
#
# Used for Figure 8 to compute the Fisher-information-based standard error
# of the unconstrained conditional MLE of b2
#
# Requires:
#  * R/01_rates.R
#  * R/generator.R
#  * R/pf.R
#  * R/mle.R
#
# This file uses:
#  * compute_tilted_rates_complete()
#  * grad_log_R_num_complete()
#  * compute_pi()
#

compute_full_score_weights_complete <- function(par, N, eps = 1e-5) {
  par <- par[c("b1", "b2", "mu")]

  rates_obj <- compute_tilted_rates_complete(par, N)
  grad_obj <- grad_log_R_num_complete(par, N, eps = eps)

  if (is.null(rates_obj) || is.null(grad_obj)) {
    return(NULL)
  }

  lambda_base <- rates_obj$lambda_base
  grad_log_R_plus <- grad_obj$grad_log_R_plus
  grad_log_R_minus <- grad_obj$grad_log_R_minus

  G_plus <- matrix(NA_real_, nrow = N, ncol = 3)
  G_minus <- matrix(NA_real_, nrow = N, ncol = 3)

  rownames(G_plus) <- rownames(G_minus) <- 1:N
  colnames(G_plus) <- colnames(G_minus) <- c("b1", "b2", "mu")

  k_birth <- 1:(N - 1)

  phi1_k <- phi1_complete(k_birth, N)
  phi2_k <- phi2_complete(k_birth, N)

  G_plus[k_birth, "b1"] <-
    phi1_k / lambda_base[k_birth] +
    grad_log_R_plus[k_birth, "b1"]

  G_plus[k_birth, "b2"] <-
    phi2_k / lambda_base[k_birth] +
    grad_log_R_plus[k_birth, "b2"]

  G_plus[k_birth, "mu"] <-
    grad_log_R_plus[k_birth, "mu"]

  k_death <- 2:N

  G_minus[k_death, "b1"] <-
    grad_log_R_minus[k_death, "b1"]

  G_minus[k_death, "b2"] <-
    grad_log_R_minus[k_death, "b2"]

  G_minus[k_death, "mu"] <-
    1 / unname(par["mu"]) +
    grad_log_R_minus[k_death, "mu"]

  list(
    G_plus = G_plus,
    G_minus = G_minus,
    rates = rates_obj
  )
}

compute_fisher_complete <- function(par, N, eps = 1e-5) {
  par <- par[c("b1", "b2", "mu")]

  obj <- compute_full_score_weights_complete(par, N, eps = eps)

  if (is.null(obj)) {
    return(NULL)
  }

  G_plus <- obj$G_plus
  G_minus <- obj$G_minus

  lambda_tilde <- obj$rates$lambda_tilde
  mu_tilde <- obj$rates$mu_tilde

  pi_tilde <- compute_pi(obj$rates$v, obj$rates$h)

  I_hat <- matrix(0, nrow = 3, ncol = 3)
  rownames(I_hat) <- colnames(I_hat) <- c("b1", "b2", "mu")

  for (k in 1:(N - 1)) {
    gk <- matrix(G_plus[k, ], ncol = 1)
    I_hat <- I_hat + pi_tilde[k] * lambda_tilde[k] * (gk %*% t(gk))
  }

  for (k in 2:N) {
    gk <- matrix(G_minus[k, ], ncol = 1)
    I_hat <- I_hat + pi_tilde[k] * mu_tilde[k] * (gk %*% t(gk))
  }

  0.5 * (I_hat + t(I_hat))
}

vcov_mle_complete <- function(par_hat, stats, N, eps = 1e-5, ridge = 0) {
  T_obs <- sum(as.numeric(stats$T_k))

  I_hat <- compute_fisher_complete(par_hat, N, eps = eps)

  if (is.null(I_hat)) {
    return(NULL)
  }

  I_reg <- I_hat + ridge * diag(3)

  V_hat <- tryCatch(
    (1 / T_obs) * solve(I_reg),
    error = function(e) NULL
  )

  if (is.null(V_hat)) {
    return(NULL)
  }

  rownames(V_hat) <- colnames(V_hat) <- c("b1", "b2", "mu")

  list(
    vcov = V_hat,
    fisher = I_hat
  )
}

se_mle_complete <- function(par_hat, stats, N, eps = 1e-5, ridge = 0) {
  V_obj <- vcov_mle_complete(par_hat, stats, N, eps = eps, ridge = ridge)
  if (is.null(V_obj)) {
    return(NULL)
  }
  sqrt(diag(V_obj$vcov))
}
