# R/qmle.R
#
# QMLE for the complete-graph SIS composite birth-death model.
#
# Requires:
#   R/01_rates.R
#   R/generator.R
#   R/pf.R
#
# Input stats should come from compute_sufficient_stats(traj, N).
# Only T_k, B_k, and D_k are used here.
# 

# small helper: extract states 1, ..., N from named vectors 0, ..., N
extract_1_to_N <- function(x, N) {
  if (!is.null(names(x))) {
    as.numeric(x[as.character(1:N)])
  } else {
    as.numeric(x[2:(N + 1)])
  }
}

# working score U_T(theta)
score_qmle_complete <- function(par, stats, N) {
  par <- par[c("b1", "b2", "mu")]

  if (any(!is.finite(par)) || any(par <= 0)) {
    return(c(U_b1 = NA_real_, U_b2 = NA_real_, U_mu = NA_real_))
  }

  T_k <- extract_1_to_N(stats$T_k, N)
  B_k <- extract_1_to_N(stats$B_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  Q_plus <- build_Q_plus_complete(par, N)
  pf <- pf_decomp_Qplus(Q_plus)
  R_obj <- compute_R_pm(pf$h)

  R_plus <- R_obj$R_plus
  R_minus <- R_obj$R_minus

  k_birth <- 1:(N - 1)
  k_death <- 2:N  # FIXED: death sum runs over k = 2, ..., N (no 1 -> 0 under the Q-process), consistent with qem.R and the docs

  phi1_k <- phi1_complete(k_birth, N)
  phi2_k <- phi2_complete(k_birth, N)
  lambda_k <- lambda_complete(k_birth, par[c("b1", "b2")], N)

  if (any(lambda_k <= 0)) {
    return(c(U_b1 = NA_real_, U_b2 = NA_real_, U_mu = NA_real_))
  }

  U_b1 <- sum(B_k[k_birth] * phi1_k / lambda_k) -
    sum(T_k[k_birth] * R_plus[k_birth] * phi1_k)

  U_b2 <- sum(B_k[k_birth] * phi2_k / lambda_k) -
    sum(T_k[k_birth] * R_plus[k_birth] * phi2_k)

  r_k <- r_complete(k_death)

  U_mu <- sum(D_k[k_death]) / par["mu"] -
    sum(T_k[k_death] * r_k * R_minus[k_death])  # FIXED: sum over k_death = 2:N instead of 1:N

  c(U_b1 = U_b1, U_b2 = U_b2, U_mu = U_mu)
}


# numerical Jacobian for the birth-score block
jacobian_beta_qmle_complete <- function(par, stats, N, eps = 1e-6) {
  par <- par[c("b1", "b2", "mu")]

  J <- matrix(NA_real_, nrow = 2, ncol = 2)
  rownames(J) <- c("U_b1", "U_b2")
  colnames(J) <- c("b1", "b2")

  for (j in c("b1", "b2")) {
    step <- eps * max(1, abs(par[j]))

    par_plus <- par
    par_minus <- par

    par_plus[j] <- par_plus[j] + step
    par_minus[j] <- par_minus[j] - step

    U_plus <- score_qmle_complete(par_plus, stats, N)[c("U_b1", "U_b2")]
    U_minus <- score_qmle_complete(par_minus, stats, N)[c("U_b1", "U_b2")]

    J[, j] <- (U_plus - U_minus) / (2 * step)
  }
  J
}


# score-norm objective, useful for diagnostics
qmle_score_norm_complete <- function(par, stats, N) {
  U <- score_qmle_complete(par, stats, N)
  if (any(!is.finite(U))) return(Inf)
  sum(U^2)
}


# block-Newton QMLE
fit_qmle_complete <- function(stats, N, init,
                              max_iter = 50,
                              tol = 1e-8,
                              eps = 1e-6,
                              ridge = 1e-8,
                              step_halving = TRUE,
                              verbose = FALSE) {
  par <- init[c("b1", "b2", "mu")]

  trace <- matrix(NA_real_, nrow = max_iter + 1, ncol = 6)
  colnames(trace) <- c("b1", "b2", "mu", "U_b1", "U_b2", "U_mu")

  U <- score_qmle_complete(par, stats, N)
  trace[1, ] <- c(par, U)

  if (verbose) {
    cat(
      "iter 0:",
      "b1 =", par["b1"],
      "b2 =", par["b2"],
      "mu =", par["mu"],
      "|score| =", sqrt(sum(U^2)),
      "\n"
    )
  }

  converged <- FALSE

  for (iter in 1:max_iter) {
    U <- score_qmle_complete(par, stats, N)

    if (any(!is.finite(U))) {
      warning("Current score is not finite; stopping.")
      break
    }

    if (max(abs(U)) < tol) {
      converged <- TRUE
      break
    }

    J_beta <- jacobian_beta_qmle_complete(par, stats, N, eps = eps)

    if (any(!is.finite(J_beta))) {
      warning("Numerical Jacobian is not finite; stopping.")
      break
    }

    J_beta_reg <- J_beta + ridge * diag(2)

    delta_beta <- tryCatch(
      solve(J_beta_reg, U[c("U_b1", "U_b2")]),
      error = function(e) rep(NA_real_, 2)
    )

    if (any(!is.finite(delta_beta))) {
      warning("Beta Newton step failed; stopping.")
      break
    }

    par_candidate <- par
    par_candidate[c("b1", "b2")] <- par[c("b1", "b2")] - delta_beta

    if (step_halving) {
      old_norm <- sum(U[c("U_b1", "U_b2")]^2)
      step_factor <- 1

      repeat {
        par_try <- par
        par_try[c("b1", "b2")] <- par[c("b1", "b2")] -
          step_factor * delta_beta

        U_try <- score_qmle_complete(par_try, stats, N)

        new_norm <- if (all(is.finite(U_try))) {
          sum(U_try[c("U_b1", "U_b2")]^2)
        } else {
          Inf
        }

        if (is.finite(new_norm) && new_norm <= old_norm) {
          par_candidate <- par_try
          break
        }

        step_factor <- step_factor / 2

        if (step_factor < 1e-6) {
          par_candidate <- par
          break
        }
      }
    }

    par[c("b1", "b2")] <- par_candidate[c("b1", "b2")]

    Q_plus <- build_Q_plus_complete(par, N)
    pf <- pf_decomp_Qplus(Q_plus)
    R_obj <- compute_R_pm(pf$h)
    R_minus <- R_obj$R_minus

    T_k <- extract_1_to_N(stats$T_k, N)
    D_k <- extract_1_to_N(stats$D_k, N)
    k_death <- 2:N  # mu update sums over k = 2, ..., N

    # sum over k_death = 2:N instead of 1:N
    denom_mu <- sum(T_k[k_death] * r_complete(k_death) * R_minus[k_death])

    if (abs(denom_mu) < .Machine$double.eps) {
      warning("mu update denominator is too small; stopping.")
      break
    }

    # numerator sums over k_death = 2:N
    par["mu"] <- sum(D_k[k_death]) / denom_mu 

    U_new <- score_qmle_complete(par, stats, N)
    trace[iter + 1, ] <- c(par, U_new)

    if (verbose) {
      cat(
        "iter", iter, ":",
        "b1 =", par["b1"],
        "b2 =", par["b2"],
        "mu =", par["mu"],
        "|score| =", sqrt(sum(U_new^2)),
        "\n"
      )
    }

    if (all(is.finite(U_new)) && max(abs(U_new)) < tol) {
      converged <- TRUE
      break
    }
  }

  last_iter <- max(which(!is.na(trace[, 1])))

  list(
    par_hat = par,
    score = score_qmle_complete(par, stats, N),
    converged = converged,
    iterations = last_iter - 1,
    trace = trace[1:last_iter, , drop = FALSE]
  )
}