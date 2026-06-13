# R/unc_mle.R
#
# Naive unconditional MLE for model with unmarked births
#
# The unconditional log-likelihood is
#   sum_k [ B_k log(lambda_k) - T_k lambda_k ] +
#   sum_k [ D_k log(mu r(k)) - T_k mu r(k) ]
#
# where
#   lambda_k = b1 phi1_complete(k, N) + b2 phi2_complete(k, N)
#   r(k) = k
#
# Requires:
#   R/01_rates.R
#
# Input stats are produced by compute_sufficient_stats() in R/03_path_statistics.R
#

# small helper: extract states 1, ..., N from named vectors 0, ..., N
extract_1_to_N <- function(x, N) {
  if (!is.null(names(x))) {
    as.numeric(x[as.character(1:N)])
  } else {
    as.numeric(x[2:(N + 1)])
  }
}

# closed-form unconditional MLE for mu
mu_hat_unc_complete <- function(stats, N) {
  T_k <- extract_1_to_N(stats$T_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_all <- 1:N

  denom <- sum(T_k[k_all] * r_complete(k_all))

  if (abs(denom) < .Machine$double.eps) {
    return(NA_real_)
  }

  sum(D_k[k_all]) / denom
}

# unconditional log-likelihood
loglik_unc_complete <- function(par, stats, N) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_birth <- 1:(N - 1)
  k_all <- 1:N

  lambda_k <- lambda_complete(k_birth, par[c("b1", "b2")], N)
  mu <- unname(par["mu"])

  if (any(lambda_k <= 0) || mu <= 0) {
    return(-Inf)
  }

  ll_birth <- sum(
    B_k[k_birth] * log(lambda_k) -
      T_k[k_birth] * lambda_k
  )

  ll_death <- sum(
    D_k[k_all] * log(mu * r_complete(k_all)) -
      T_k[k_all] * mu * r_complete(k_all)
  )

  ll_birth + ll_death
}

# unconditional score
score_unc_complete <- function(par, stats, N) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_birth <- 1:(N - 1)
  k_all <- 1:N

  phi1_k <- phi1_complete(k_birth, N)
  phi2_k <- phi2_complete(k_birth, N)
  lambda_k <- lambda_complete(k_birth, par[c("b1", "b2")], N)
  mu <- unname(par["mu"])

  if (any(lambda_k <= 0) || abs(mu) < .Machine$double.eps) {
    return(c(U_b1 = NA_real_, U_b2 = NA_real_, U_mu = NA_real_))
  }

  U_b1 <- sum(B_k[k_birth] * phi1_k / lambda_k) -
    sum(T_k[k_birth] * phi1_k)

  U_b2 <- sum(B_k[k_birth] * phi2_k / lambda_k) -
    sum(T_k[k_birth] * phi2_k)

  U_mu <- sum(D_k[k_all]) / mu -
    sum(T_k[k_all] * r_complete(k_all))

  c(U_b1 = U_b1, U_b2 = U_b2, U_mu = U_mu)
}


# log-likelihood in (b1, b2), plugging in mu_hat
loglik_unc_profile_complete <- function(b, stats, N) {
  par <- c(b1 = b[1], b2 = b[2], mu = mu_hat_unc_complete(stats, N))
  loglik_unc_complete(par, stats, N)
}

# fit unconditional MLE
fit_unc_mle_complete <- function(stats, N, init,
                                 method = "BFGS",
                                 control = list(maxit = 200)) {
  init <- init[c("b1", "b2", "mu")]

  mu_hat <- mu_hat_unc_complete(stats, N)

  if (!is.finite(mu_hat) || mu_hat <= 0) {
    stop("Closed-form unconditional mu estimate is not well-defined.")
  }

  obj_fn <- function(b_vec) {
    ll <- loglik_unc_profile_complete(b_vec, stats, N)
    if (!is.finite(ll)) return(1e12)
    -ll
  }

  opt <- optim(
    par = unname(init[c("b1", "b2")]),
    fn = obj_fn,
    method = method,
    control = control
  )

  par_hat <- c(
    b1 = opt$par[1],
    b2 = opt$par[2],
    mu = mu_hat
  )

  list(
    par_hat = par_hat,
    loglik = loglik_unc_complete(par_hat, stats, N),
    score = score_unc_complete(par_hat, stats, N),
    convergence = opt$convergence,
    value = opt$value,
    counts = opt$counts,
    message = opt$message,
    optim = opt
  )
}
