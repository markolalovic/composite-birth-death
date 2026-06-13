# R/mle.R
#
# Full conditional MLE
#
# This file implements:
#   - tilted rates via the Doob h-transform
#   - numerical derivatives of log R^+ and log R^-
#   - the full score, including cross-sensitivity terms
#   - the full conditional log-likelihood
#   - an optimizer-based MLE fit
#
# Requires:
#   R/01_rates.R
#   R/generator.R
#   R/pf.R
#
# Input stats should come from compute_sufficient_stats() in
# R/03_path_statistics.R.
#

# small helper: extract states 1, ..., N from named vectors 0, ..., N
extract_1_to_N <- function(x, N) {
  if (!is.null(names(x))) {
    as.numeric(x[as.character(1:N)])
  } else {
    as.numeric(x[2:(N + 1)])
  }
}

# compute log R^+, log R^-
compute_log_R_pm_complete <- function(par, N) {
  par <- par[c("b1", "b2", "mu")]

  mu <- unname(par["mu"])
  lambda_base <- lambda_complete(1:(N - 1), par[c("b1", "b2")], N)

  if (!is.finite(mu) || mu <= 0) {
    return(NULL)
  }

  if (any(!is.finite(lambda_base)) || any(lambda_base <= 0)) {
    return(NULL)
  }

  Q_plus <- build_Q_plus_complete(par, N)
  pf <- pf_decomp_Qplus(Q_plus)
  R_obj <- compute_R_pm(pf$h)

  R_plus <- R_obj$R_plus
  R_minus <- R_obj$R_minus

  log_R_plus <- rep(NA_real_, N)
  log_R_minus <- rep(NA_real_, N)

  if (any(R_plus[1:(N - 1)] <= 0, na.rm = TRUE)) {
    return(NULL)
  }

  if (any(R_minus[2:N] <= 0, na.rm = TRUE)) {
    return(NULL)
  }

  log_R_plus[1:(N - 1)] <- log(R_plus[1:(N - 1)])
  log_R_minus[2:N] <- log(R_minus[2:N])

  names(log_R_plus) <- 1:N
  names(log_R_minus) <- 1:N

  list(
    gamma = pf$gamma,
    h = pf$h,
    v = pf$v,
    R_plus = R_plus,
    R_minus = R_minus,
    log_R_plus = log_R_plus,
    log_R_minus = log_R_minus
  )
}


# tilted rates
compute_tilted_rates_complete <- function(par, N) {
  par <- par[c("b1", "b2", "mu")]

  R_obj <- compute_log_R_pm_complete(par, N)
  if (is.null(R_obj)) return(NULL)

  lambda_base <- rep(NA_real_, N)
  mu_base <- rep(NA_real_, N)

  lambda_tilde <- rep(NA_real_, N)
  mu_tilde <- rep(NA_real_, N)

  lambda_base[1:(N - 1)] <- lambda_complete(
    1:(N - 1),
    par[c("b1", "b2")],
    N
  )

  lambda_tilde[1:(N - 1)] <- lambda_base[1:(N - 1)] *
    R_obj$R_plus[1:(N - 1)]

  mu_base[2:N] <- mu_complete(2:N, unname(par["mu"]))
  mu_tilde[2:N] <- mu_base[2:N] * R_obj$R_minus[2:N]

  if (any(lambda_tilde[1:(N - 1)] <= 0, na.rm = TRUE)) return(NULL)
  if (any(mu_tilde[2:N] <= 0, na.rm = TRUE)) return(NULL)

  names(lambda_base) <- 1:N
  names(mu_base) <- 1:N
  names(lambda_tilde) <- 1:N
  names(mu_tilde) <- 1:N

  list(
    gamma = R_obj$gamma,
    h = R_obj$h,
    v = R_obj$v,
    R_plus = R_obj$R_plus,
    R_minus = R_obj$R_minus,
    log_R_plus = R_obj$log_R_plus,
    log_R_minus = R_obj$log_R_minus,
    lambda_base = lambda_base,
    mu_base = mu_base,
    lambda_tilde = lambda_tilde,
    mu_tilde = mu_tilde
  )
}

# numerical gradients of log R^+, log R^-
grad_log_R_num_complete <- function(par, N, eps = 1e-5) {
  par <- par[c("b1", "b2", "mu")]
  param_names <- c("b1", "b2", "mu")

  grad_log_R_plus <- matrix(NA_real_, nrow = N, ncol = 3)
  grad_log_R_minus <- matrix(NA_real_, nrow = N, ncol = 3)

  rownames(grad_log_R_plus) <- rownames(grad_log_R_minus) <- 1:N
  colnames(grad_log_R_plus) <- colnames(grad_log_R_minus) <- param_names

  for (j in param_names) {
    step <- eps * max(1, abs(par[j]))

    par_plus <- par
    par_minus <- par

    par_plus[j] <- par_plus[j] + step
    par_minus[j] <- par_minus[j] - step

    obj_plus <- compute_log_R_pm_complete(par_plus, N)
    obj_minus <- compute_log_R_pm_complete(par_minus, N)

    if (is.null(obj_plus) || is.null(obj_minus)) {
      return(NULL)
    }

    grad_log_R_plus[, j] <-
      (obj_plus$log_R_plus - obj_minus$log_R_plus) / (2 * step)

    grad_log_R_minus[, j] <-
      (obj_plus$log_R_minus - obj_minus$log_R_minus) / (2 * step)
  }

  list(
    grad_log_R_plus = grad_log_R_plus,
    grad_log_R_minus = grad_log_R_minus
  )
}


# full conditional log-likelihood
loglik_mle_complete <- function(par, stats, N) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  rates_obj <- compute_tilted_rates_complete(par, N)
  if (is.null(rates_obj)) return(-Inf)

  lambda_tilde <- rates_obj$lambda_tilde
  mu_tilde <- rates_obj$mu_tilde

  ll_birth <- sum(
    B_k[1:(N - 1)] * log(lambda_tilde[1:(N - 1)]) -
      T_k[1:(N - 1)] * lambda_tilde[1:(N - 1)]
  )

  ll_death <- sum(
    D_k[2:N] * log(mu_tilde[2:N]) -
      T_k[2:N] * mu_tilde[2:N]
  )

  ll_birth + ll_death
}

# full score
score_mle_complete <- function(par, stats, N, eps = 1e-5) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  rates_obj <- compute_tilted_rates_complete(par, N)
  grad_obj <- grad_log_R_num_complete(par, N, eps = eps)

  if (is.null(rates_obj) || is.null(grad_obj)) {
    return(c(U_b1 = NA_real_, U_b2 = NA_real_, U_mu = NA_real_))
  }

  lambda_base <- rates_obj$lambda_base
  lambda_tilde <- rates_obj$lambda_tilde
  mu_tilde <- rates_obj$mu_tilde

  grad_log_R_plus <- grad_obj$grad_log_R_plus
  grad_log_R_minus <- grad_obj$grad_log_R_minus

  phi1_k <- rep(NA_real_, N)
  phi2_k <- rep(NA_real_, N)

  phi1_k[1:(N - 1)] <- phi1_complete(1:(N - 1), N)
  phi2_k[1:(N - 1)] <- phi2_complete(1:(N - 1), N)

  res_birth <- rep(NA_real_, N)
  res_death <- rep(NA_real_, N)

  res_birth[1:(N - 1)] <- B_k[1:(N - 1)] -
    T_k[1:(N - 1)] * lambda_tilde[1:(N - 1)]

  res_death[2:N] <- D_k[2:N] -
    T_k[2:N] * mu_tilde[2:N]

  g_plus_b1 <- rep(NA_real_, N)
  g_plus_b2 <- rep(NA_real_, N)
  g_plus_mu <- rep(NA_real_, N)

  g_minus_b1 <- rep(NA_real_, N)
  g_minus_b2 <- rep(NA_real_, N)
  g_minus_mu <- rep(NA_real_, N)

  g_plus_b1[1:(N - 1)] <- phi1_k[1:(N - 1)] /
    lambda_base[1:(N - 1)] +
    grad_log_R_plus[1:(N - 1), "b1"]

  g_plus_b2[1:(N - 1)] <- phi2_k[1:(N - 1)] /
    lambda_base[1:(N - 1)] +
    grad_log_R_plus[1:(N - 1), "b2"]

  g_plus_mu[1:(N - 1)] <- grad_log_R_plus[1:(N - 1), "mu"]

  g_minus_b1[2:N] <- grad_log_R_minus[2:N, "b1"]
  g_minus_b2[2:N] <- grad_log_R_minus[2:N, "b2"]
  g_minus_mu[2:N] <- 1 / unname(par["mu"]) +
    grad_log_R_minus[2:N, "mu"]

  U_b1 <- sum(res_birth[1:(N - 1)] * g_plus_b1[1:(N - 1)]) +
    sum(res_death[2:N] * g_minus_b1[2:N])

  U_b2 <- sum(res_birth[1:(N - 1)] * g_plus_b2[1:(N - 1)]) +
    sum(res_death[2:N] * g_minus_b2[2:N])

  U_mu <- sum(res_birth[1:(N - 1)] * g_plus_mu[1:(N - 1)]) +
    sum(res_death[2:N] * g_minus_mu[2:N])

  c(U_b1 = U_b1, U_b2 = U_b2, U_mu = U_mu)
}


# optimizer-based full MLE
fit_mle_complete <- function(stats, N, init,
                             method = "BFGS",
                             control = list(maxit = 200),
                             eps = 1e-5) {
  init <- init[c("b1", "b2", "mu")]

  obj_fn <- function(par_vec) {
    par <- c(b1 = par_vec[1], b2 = par_vec[2], mu = par_vec[3])
    ll <- loglik_mle_complete(par, stats, N)

    if (!is.finite(ll)) return(1e12)

    -ll
  }

  opt <- optim(
    par = unname(init),
    fn = obj_fn,
    method = method,
    control = control
  )

  par_hat <- c(b1 = opt$par[1], b2 = opt$par[2], mu = opt$par[3])
  list(
    par_hat = par_hat,
    loglik = loglik_mle_complete(par_hat, stats, N),
    score = score_mle_complete(par_hat, stats, N, eps = eps),
    convergence = opt$convergence,
    value = opt$value,
    counts = opt$counts,
    message = opt$message,
    optim = opt
  )
}
