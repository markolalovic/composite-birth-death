# R/unc_mle_marked.R
#
# Naive unconditional MLE
#
# The birth mechanism is observed for every birth event:
#   U_k = number of pairwise births from state k,
#   V_k = number of triadic births from state k.
#
# Closed-form MLEs:
#   b1_hat = sum_k U_k / sum_k T_k phi1(k)
#   b2_hat = sum_k V_k / sum_k T_k phi2(k)
#   mu_hat = sum_k D_k / sum_k T_k r(k)
#
# Requires:
#   R/01_rates.R
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

fit_unc_mle_marked_complete <- function(stats, N) {
  T_k <- extract_1_to_N(stats$T_k, N)
  U_k <- extract_1_to_N(stats$U_k, N)
  V_k <- extract_1_to_N(stats$V_k, N)
  D_k <- extract_1_to_N(stats$D_k, N)

  k_birth <- 1:(N - 1)
  k_all <- 1:N

  b1_hat <- sum(U_k[k_birth]) /
    sum(T_k[k_birth] * phi1_complete(k_birth, N))

  b2_hat <- sum(V_k[k_birth]) /
    sum(T_k[k_birth] * phi2_complete(k_birth, N))

  mu_hat <- sum(D_k[k_all]) /
    sum(T_k[k_all] * r_complete(k_all))

  list(
    par_hat = c(b1 = b1_hat, b2 = b2_hat, mu = mu_hat)
  )
}