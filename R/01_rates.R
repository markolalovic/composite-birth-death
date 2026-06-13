# src/01_rates.R
# 
# Complete-graph SIS model with two birth mechanisms:
#   - pairwise (PW)
#   - higher-order / triadic (HO)
#
# We work in scaled coordinates:
#   b1 = N * beta1
#   b2 = N^2 * beta2
#
# Hence the total birth rate is
#   lambda_k(b) = b1 * phi1_complete(k, N) + b2 * phi2_complete(k, N)
# where:
#   phi1_complete(k, N) = f1_complete(k, N) / N
#   phi2_complete(k, N) = f2_complete(k, N) / N^2
#
# Recovery rate:
#   mu_k(mu) = mu * r_complete(k),   with r_complete(k) = k
# 

# unscaled birth basis functions
f1_complete <- function(k, N) {
  ifelse(k > 0 & k < N, k * (N - k), 0)
}

f2_complete <- function(k, N) {
  ifelse(k > 1 & k < N, 0.5 * k * (k - 1) * (N - k), 0)
}


# scaled birth basis functions
phi1_complete <- function(k, N) {
  f1_complete(k, N) / N
}

phi2_complete <- function(k, N) {
  f2_complete(k, N) / (N^2)
}

# recovery basis function
r_complete <- function(k) {
  k
}


# rates
lambda_complete <- function(k, b, N) {
  # b must be named vector with entries b1, b2
  b[["b1"]] * phi1_complete(k, N) + b[["b2"]] * phi2_complete(k, N)
}

mu_complete <- function(k, mu) {
  mu * r_complete(k)
}

total_rate_complete <- function(k, b, mu, N) {
  lambda_complete(k, b, N) + mu_complete(k, mu)
}
