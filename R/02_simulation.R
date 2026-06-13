# src/02_simulation.R
#
# Gillespie simulation for the complete-graph SIS model
#
# The model functions f1_complete(), f2_complete(), phi1_complete(),
# phi2_complete(), r_complete(), lambda_complete(), and mu_complete()
# are defined in src/01_rates.R
#
# Parameter: theta = c(b1 = ..., b2 = ..., mu = ...)
#
# where:
#   * b1 = N beta1,
#   * b2 = N^2 beta2.
#
# The event_type column records birth marks:
#   * "PW_birth" for pairwise births,
#   * "HO_birth" for triadic births.
#
# For unmarked-birth inference, ignore event_type 
# except for distinguishing births from recoveries
#

simulate_path <- function(theta, N, I0, time_max) {
  k <- I0
  t <- 0.0

  history_list <- vector("list", 1000)

  history_list[[1]] <- list(
    time = t,
    k = k,
    event_type = "initial",
    pre_jump_k = NA,
    jump_type = 0
  )

  i <- 2
  while (t < time_max && k > 0) {

    rate_pw_birth <- theta[["b1"]] * phi1_complete(k, N)
    rate_ho_birth <- theta[["b2"]] * phi2_complete(k, N)
    rate_recovery <- theta[["mu"]] * r_complete(k)

    total_rate <- rate_pw_birth + rate_ho_birth + rate_recovery

    if (total_rate < .Machine$double.eps) {
      break
    }

    t <- t + rexp(1, rate = total_rate)

    if (t >= time_max) {
      break
    }

    event_draw <- sample(
      c("PW_birth", "HO_birth", "Recovery"),
      size = 1,
      prob = c(rate_pw_birth, rate_ho_birth, rate_recovery)
    )

    pre_jump_k <- k

    if (event_draw == "PW_birth" || event_draw == "HO_birth") {
      k <- k + 1
      jump_type <- 1
    } else {
      k <- k - 1
      jump_type <- -1
    }

    if (i > length(history_list)) {
      length(history_list) <- 2 * length(history_list)
    }

    history_list[[i]] <- list(
      time = t,
      k = k,
      event_type = event_draw,
      pre_jump_k = pre_jump_k,
      jump_type = jump_type
    )

    i <- i + 1
  }

  history_list[[i]] <- list(
    time = time_max,
    k = k,
    event_type = "end",
    pre_jump_k = NA,
    jump_type = 0
  )

  dplyr::bind_rows(history_list)
}