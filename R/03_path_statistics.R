# src/03_path_statistics.R

compute_sufficient_stats <- function(traj, N) {
  # initialize
  states <- 0:N

  T_k <- setNames(numeric(N + 1), states)
  U_k <- setNames(numeric(N + 1), states)
  V_k <- setNames(numeric(N + 1), states)
  D_k <- setNames(numeric(N + 1), states)

  # durations T_k
  if (nrow(traj) > 1) {
    durations <- diff(traj$time)

    # traj$k[i] holds for interval (time[i], time[i+1])
    current_states <- traj$k[-nrow(traj)]

    # sum durations by state
    dur_sums <- tapply(
      durations,
      factor(current_states, levels = states),
      sum)

    # tapply returns NA for unobserved levels
    dur_sums[is.na(dur_sums)] <- 0
    T_k[] <- dur_sums
  }

  # event counts U_k, V_k, D_k by pre-jump state
  events <- traj[!is.na(traj$pre_jump_k), ]

  # the column selection needs to be explicit
  # if some birth type never happens
  if (nrow(events) > 0) {
    counts <- table(
      factor(events$pre_jump_k, levels = states),
      events$event_type
    )

    if ("PW_birth" %in% colnames(counts)) U_k[] <- counts[, "PW_birth"]
    if ("HO_birth" %in% colnames(counts)) V_k[] <- counts[, "HO_birth"]
    if ("Recovery" %in% colnames(counts)) D_k[] <- counts[, "Recovery"]
  }

  list(
    T_k = T_k,      # waiting times 
    U_k = U_k,      # PW births
    V_k = V_k,      # HO births
    D_k = D_k,      # deaths
    B_k = U_k + V_k # for unmarked process: total births
  )
}