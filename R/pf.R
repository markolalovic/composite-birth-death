# R/pf.R
#
# Perron-Frobenius objects for the killed generator Q^+
#
# Q_plus is the killed generator on states {1, ..., N}
# 
# `gamma` is the maximal real eigenvalue of Q_plus, so gamma < 0
#
# h is the right eigenvector:
#  Q_plus %*% h = gamma * h
# 
# v is the left eigenvector:
#       t(v) %*% Q_plus = gamma * t(v)
# The normalization is sum(v * h) = 1
#
# Vectors h, v, R_plus, R_minus, and pi are indexed by states 1, ..., N
#

pf_decomp_Qplus <- function(Q_plus) {
  N <- nrow(Q_plus)

  if (ncol(Q_plus) != N) {
    stop("Q_plus must be a square matrix.")
  }

  eig_r <- eigen(Q_plus)
  idx_r <- which.max(Re(eig_r$values))

  gamma <- Re(eig_r$values[idx_r])
  h <- Re(eig_r$vectors[, idx_r])

  eig_l <- eigen(t(Q_plus))
  idx_l <- which.max(Re(eig_l$values))

  v <- Re(eig_l$vectors[, idx_l])

  if (sum(h) < 0) h <- -h
  if (sum(v) < 0) v <- -v

  scale_vh <- sum(v * h)

  h <- h / sqrt(scale_vh)
  v <- v / sqrt(scale_vh)

  states <- 1:N
  names(h) <- states
  names(v) <- states

  list(
    gamma = gamma,
    h = h,
    v = v
  )
}

compute_R_pm <- function(h) {
  N <- length(h)

  h_ext <- c(0, h)

  R_plus <- rep(NA_real_, N)
  R_minus <- rep(NA_real_, N)

  if (N >= 2) {
    R_plus[1:(N - 1)] <- h[2:N] / h[1:(N - 1)]
  }

  R_minus[1] <- 0

  if (N >= 2) {
    R_minus[2:N] <- h_ext[2:N] / h[2:N]
  }

  states <- 1:N
  names(R_plus) <- states
  names(R_minus) <- states

  list(
    R_plus = R_plus,
    R_minus = R_minus
  )
}

compute_pi <- function(v, h) {
  pi <- v * h
  pi <- pi / sum(pi)

  names(pi) <- names(h)

  pi
}