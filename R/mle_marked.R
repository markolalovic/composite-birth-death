# R/mle_marked.R
#
# Full conditional MLE with marked births
#
# The birth mechanism is observed for each birth event:
#   U_k = number of pairwise births from state k
#   V_k = number of triadic births from state k
#
# Conditional marked intensities under the Q-process:
#   pairwise birth: b1 * phi1(k) * R_plus(k)
#   triadic birth:  b2 * phi2(k) * R_plus(k)
#   death:          mu * r(k)    * R_minus(k)
#
# Requires:
#   R/01_rates.R
#   R/generator.R
#   R/pf.R
#   R/mle.R
#
# Input stats should come from compute_sufficient_stats() in
# R/03_path_statistics.R.
#

# some helpers
extract_1_to_N <- function(x, N) {
  if (!is.null(names(x))) {
    as.numeric(x[as.character(1:N)])
  } else {
    as.numeric(x[2:(N + 1)])
  }
}

safe_xlog <- function(x, y) {
  result <- rep(0, length(x))
  pos <- x > 0
  result[pos] <- x[pos] * log(y[pos])
  result
}


loglik_mle_marked_complete <- function(par, stats, N) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  U_k <- extract_1_to_N(stats$U_k, N)
  V_k <- extract_1_to_N(stats$V_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  rates_obj <- compute_tilted_rates_complete(par, N)
  if (is.null(rates_obj)) return(-Inf)

  R_plus <- rates_obj$R_plus
  mu_tilde <- rates_obj$mu_tilde

  k_birth <- 1:(N - 1)

  phi1 <- phi1_complete(k_birth, N)
  phi2 <- phi2_complete(k_birth, N)

  lam1_tilde <- unname(par["b1"]) * phi1 * R_plus[k_birth]
  lam2_tilde <- unname(par["b2"]) * phi2 * R_plus[k_birth]

  if (any(U_k[k_birth] > 0 & lam1_tilde <= 0)) return(-Inf)
  if (any(V_k[k_birth] > 0 & lam2_tilde <= 0)) return(-Inf)
  if (any(D_k[2:N] > 0 & mu_tilde[2:N] <= 0)) return(-Inf)

  ll_birth1 <- sum(
    safe_xlog(U_k[k_birth], lam1_tilde) -
      T_k[k_birth] * lam1_tilde
  )

  ll_birth2 <- sum(
    safe_xlog(V_k[k_birth], lam2_tilde) -
      T_k[k_birth] * lam2_tilde
  )

  ll_death <- sum(
    safe_xlog(D_k[2:N], mu_tilde[2:N]) -
      T_k[2:N] * mu_tilde[2:N]
  )

  ll_birth1 + ll_birth2 + ll_death
}

score_mle_marked_complete <- function(par, stats, N, eps = 1e-5) {
  par <- par[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  U_k <- extract_1_to_N(stats$U_k, N)
  V_k <- extract_1_to_N(stats$V_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  rates_obj <- compute_tilted_rates_complete(par, N)
  grad_obj <- grad_log_R_num_complete(par, N, eps = eps)

  if (is.null(rates_obj) || is.null(grad_obj)) {
    return(c(U_b1 = NA_real_, U_b2 = NA_real_, U_mu = NA_real_))
  }

  R_plus <- rates_obj$R_plus
  mu_tilde <- rates_obj$mu_tilde

  grad_log_R_plus <- grad_obj$grad_log_R_plus
  grad_log_R_minus <- grad_obj$grad_log_R_minus

  k_birth <- 1:(N - 1)

  phi1 <- phi1_complete(k_birth, N)
  phi2 <- phi2_complete(k_birth, N)

  lam1_tilde <- rep(0, N)
  lam2_tilde <- rep(0, N)

  lam1_tilde[k_birth] <- unname(par["b1"]) * phi1 * R_plus[k_birth]
  lam2_tilde[k_birth] <- unname(par["b2"]) * phi2 * R_plus[k_birth]

  res1 <- rep(0, N)
  res2 <- rep(0, N)

  res1[k_birth] <- U_k[k_birth] - T_k[k_birth] * lam1_tilde[k_birth]
  res2[k_birth] <- V_k[k_birth] - T_k[k_birth] * lam2_tilde[k_birth]

  res_birth <- res1 + res2

  res_death <- rep(0, N)
  res_death[2:N] <- D_k[2:N] - T_k[2:N] * mu_tilde[2:N]

  U_b1 <- sum(res1[k_birth]) / unname(par["b1"]) +
    sum(res_birth[k_birth] * grad_log_R_plus[k_birth, "b1"]) +
    sum(res_death[2:N] * grad_log_R_minus[2:N, "b1"])

  U_b2 <- sum(res2[k_birth]) / unname(par["b2"]) +
    sum(res_birth[k_birth] * grad_log_R_plus[k_birth, "b2"]) +
    sum(res_death[2:N] * grad_log_R_minus[2:N, "b2"])

  U_mu <- sum(res_birth[k_birth] * grad_log_R_plus[k_birth, "mu"]) +
    sum(res_death[2:N]) / unname(par["mu"]) +
    sum(res_death[2:N] * grad_log_R_minus[2:N, "mu"])

  c(U_b1 = U_b1, U_b2 = U_b2, U_mu = U_mu)
}

fit_mle_marked_complete <- function(stats, N, init,
                                    method = "Nelder-Mead",
                                    control = list(maxit = 1000),
                                    eps = 1e-5) {
  init <- init[c("b1", "b2", "mu")]

  obj_fn <- function(par_vec) {
    if (any(par_vec <= 0)) return(1e12)

    par <- c(b1 = par_vec[1], b2 = par_vec[2], mu = par_vec[3])
    ll <- loglik_mle_marked_complete(par, stats, N)

    if (!is.finite(ll)) return(1e12)

    -ll
  }

  opt <- optim(
    par = unname(init),
    fn = obj_fn,
    method = method,
    control = control
  )

  par_hat <- c(
    b1 = opt$par[1],
    b2 = opt$par[2],
    mu = opt$par[3]
  )

  list(
    par_hat = par_hat,
    loglik = loglik_mle_marked_complete(par_hat, stats, N),
    score = score_mle_marked_complete(par_hat, stats, N, eps = eps),
    convergence = opt$convergence,
    value = opt$value,
    counts = opt$counts,
    message = opt$message,
    optim = opt
  )
}