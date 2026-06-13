# R/generator.R
#
# Build the killed generator Q^+ on states {1, ..., N}
# for the complete-graph SIS composite birth-death model.
#
# Requires:
#   R/01_rates.R
#

build_Q_plus_complete <- function(par, N) {
  k_vector <- 1:N

  b <- par[c("b1", "b2")]
  mu <- par["mu"]

  b_rate <- lambda_complete(k_vector, b, N)
  d_rate <- mu_complete(k_vector, mu)

  Q_plus <- matrix(0, nrow = N, ncol = N)

  for (k in k_vector) {
    # diagonal is the total rate out of state k in the original process
    Q_plus[k, k] <- -(b_rate[k] + d_rate[k])

    # Birth: k -> k + 1.
    if (k < N) {
      Q_plus[k, k + 1] <- b_rate[k]
    }

    # death within {1, ..., N}: k -> k - 1
    # for k = 1 death goes to 0 and is not included off-diagonal
    if (k > 1) {
      Q_plus[k, k - 1] <- d_rate[k]
    }
  }
  Q_plus
}