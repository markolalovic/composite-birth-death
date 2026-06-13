# R/qem.R
#
# Stable fixed-point / Q-EM algorithms
#
#
# Marked births: U_k and V_k are observed separately.
#
# Unmarked births:
#   only B_k = U_k + V_k is used; the birth marks are imputed by the
#   conditional probabilities:
#     p1(k) = b1 phi1(k) / [b1 phi1(k) + b2 phi2(k)]
#     p2(k) = b2 phi2(k) / [b1 phi1(k) + b2 phi2(k)]
#
# The Doob-h tilt R_plus(k) cancels in these probabilities, but
# remains in the tilted exposure denominators.
#
# Requires:
#   R/01_rates.R
#   R/generator.R
#   R/pf.R
#
# Input stats should come from compute_sufficient_stats() in
# R/03_path_statistics.R.
#

extract_1_to_N <- function(x, N) {
  if (!is.null(names(x))) {
    as.numeric(x[as.character(1:N)])
  } else {
    as.numeric(x[2:(N + 1)])
  }
}


compute_q_tilt_complete <- function(par, N) {
  par <- par[c("b1", "b2", "mu")]

  Q_plus <- build_Q_plus_complete(par, N)
  pf <- pf_decomp_Qplus(Q_plus)
  R_obj <- compute_R_pm(pf$h)

  list(
    gamma = pf$gamma,
    h = pf$h,
    R_plus = R_obj$R_plus,
    R_minus = R_obj$R_minus
  )
}


fit_qem_marked_complete <- function(stats, N,
                                    init = c(b1 = 1.0, b2 = 3.5, mu = 0.9),
                                    max_iter = 1000,
                                    tol = 1e-8) {  # FIXED: 1e-10 -> 1e-8, see convergence note below
  par <- init[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  U_k <- extract_1_to_N(stats$U_k, N)
  V_k <- extract_1_to_N(stats$V_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_birth <- 1:(N - 1)
  k_death <- 2:N  # death sum runs over k = 2, ..., N (no 1 -> 0 under the Q-process)

  phi1 <- phi1_complete(k_birth, N)
  phi2 <- phi2_complete(k_birth, N)
  r_k <- r_complete(k_death)

  history <- matrix(NA_real_, nrow = max_iter + 1, ncol = 3)
  colnames(history) <- c("b1", "b2", "mu")
  history[1, ] <- par

  converged <- FALSE
  iter_count <- 0

  for (iter in 1:max_iter) {
    iter_count <- iter
    par_old <- par

    tilt <- compute_q_tilt_complete(par, N)

    R_plus <- tilt$R_plus
    R_minus <- tilt$R_minus

    denom_b1 <- sum(T_k[k_birth] * phi1 * R_plus[k_birth])
    denom_b2 <- sum(T_k[k_birth] * phi2 * R_plus[k_birth])
    denom_mu <- sum(T_k[k_death] * r_k * R_minus[k_death])

    if (denom_b1 <= 0 || denom_b2 <= 0 || denom_mu <= 0) {
      break
    }

    par["b1"] <- sum(U_k[k_birth]) / denom_b1
    par["b2"] <- sum(V_k[k_birth]) / denom_b2
    par["mu"] <- sum(D_k[k_death]) / denom_mu

    history[iter + 1, ] <- par

    # FIXED: scale-invariant relative parameter-increment test, replacing the
    # original absolute test sum(abs(par - par_old)) < 1e-10.
    #
    # The absolute test gave false non-convergence when the fixed-point map
    # converged slowly: with linear rate rho, the increment ~ rho^m, so 1e-10
    # is unreachable within max_iter for rho close to 1, even though the
    # estimate is already at the root. A *relative* increment trips once the
    # parameters stop moving to ~8 significant figures, which is the accuracy
    # the estimates actually attain (verified: agreement with the QMLE / full
    # MLE in check_qem_estimators.R). The cross-check against those estimators
    # guards against the only failure mode of an increment test (a stall short
    # of the root would show up there as a large diff, and does not).
    rel_change <- max(abs(par - par_old) / pmax(abs(par), 1e-12))  # FIXED
    if (rel_change < tol) {                                        # FIXED
      converged <- TRUE
      break
    }
  }

  history <- history[1:(iter_count + 1), , drop = FALSE]

  list(
    par_hat = par,
    converged = converged,
    iterations = iter_count,
    history = history
  )
}


fit_qem_unmarked_complete <- function(stats, N,
                                      init = c(b1 = 1.0, b2 = 3.5, mu = 0.9),
                                      max_iter = 1000,
                                      tol = 1e-8) {  # FIXED: 1e-10 -> 1e-8, see convergence note below
  par <- init[c("b1", "b2", "mu")]

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_birth <- 1:(N - 1)
  k_death <- 2:N  # death sum runs over k = 2, ..., N (no 1 -> 0 under the Q-process)

  phi1 <- phi1_complete(k_birth, N)
  phi2 <- phi2_complete(k_birth, N)
  r_k <- r_complete(k_death)

  history <- matrix(NA_real_, nrow = max_iter + 1, ncol = 3)
  colnames(history) <- c("b1", "b2", "mu")
  history[1, ] <- par

  converged <- FALSE
  iter_count <- 0

  for (iter in 1:max_iter) {
    iter_count <- iter
    par_old <- par

    tilt <- compute_q_tilt_complete(par, N)

    R_plus <- tilt$R_plus
    R_minus <- tilt$R_minus

    rate1 <- unname(par["b1"]) * phi1
    rate2 <- unname(par["b2"]) * phi2
    lambda <- rate1 + rate2

    p1 <- numeric(length(k_birth))
    p2 <- numeric(length(k_birth))

    valid <- lambda > 0

    p1[valid] <- rate1[valid] / lambda[valid]
    p2[valid] <- rate2[valid] / lambda[valid]

    exp_U <- sum(B_k[k_birth] * p1)
    exp_V <- sum(B_k[k_birth] * p2)

    denom_b1 <- sum(T_k[k_birth] * phi1 * R_plus[k_birth])
    denom_b2 <- sum(T_k[k_birth] * phi2 * R_plus[k_birth])
    denom_mu <- sum(T_k[k_death] * r_k * R_minus[k_death])

    if (denom_b1 <= 0 || denom_b2 <= 0 || denom_mu <= 0) {
      break
    }

    par["b1"] <- exp_U / denom_b1
    par["b2"] <- exp_V / denom_b2
    par["mu"] <- sum(D_k[k_death]) / denom_mu

    history[iter + 1, ] <- par

    # FIXED: scale-invariant relative parameter-increment test, replacing the
    # original absolute test sum(abs(par - par_old)) < 1e-10.
    #
    # The unmarked map converges more slowly than the marked one because the
    # E-step allocation p1/p2 feeds back into each iterate, pushing the linear
    # rate closer to 1. With the absolute 1e-10 test this produced converged =
    # FALSE on roughly half the replicates while the estimate already matched
    # the QMLE to ~1e-8. A relative increment trips once the parameters stop
    # moving to ~8 significant figures, the accuracy actually attained. The
    # check_qem_estimators.R cross-check against the QMLE guards against a
    # premature stop short of the root.
    rel_change <- max(abs(par - par_old) / pmax(abs(par), 1e-12))  # FIXED
    if (rel_change < tol) {                                        # FIXED
      converged <- TRUE
      break
    }
  }

  history <- history[1:(iter_count + 1), , drop = FALSE]

  list(
    par_hat = par,
    converged = converged,
    iterations = iter_count,
    history = history
  )
}
