# scripts/fig03_single_trajectory.R
#
# For Figure 3:
# simulated surviving trajectory for the complete-graph simplicial SIS model.
#

source("R/01_rates.R")
source("R/02_simulation.R")

set.seed(2227)

# parameters
N <- 100
I0 <- 10
T_max <- 100

theta0 <- c(b1 = 1.01, b2 = 3.70, mu = 1)

# Simulate until survival up to T_max.
repeat {
  traj <- simulate_path(theta0, N = N, I0 = I0, time_max = T_max)
  if (tail(traj$k, 1) > 0) break
}

# Plotting style
col_main <- "#4994df"

pdf_family    <- "Times"
pdf_pointsize <- 12

axis_cex <- 1.05
lab_cex  <- 1.10
axis_lwd <- 0.8
line_lwd <- 1.4

pdf(
  "figures/single_trajectory_base.pdf",
  width = 4.8,
  height = 3.1,
  family = pdf_family,
  pointsize = pdf_pointsize,
  useDingbats = FALSE
)

par(
  mar = c(3.6, 4.2, 0.5, 1.0),
  mgp = c(2.2, 0.7, 0),
  tcl = -0.25,
  las = 1
)

plot(
  NA,
  xlim = c(0, T_max + 2),
  ylim = c(0, N),
  xlab = "Time",
  ylab = "Infected count",
  xaxs = "i",
  yaxs = "i",
  axes = FALSE,
  cex.lab = lab_cex
)

axis(
  1,
  at = seq(0, T_max, by = 10),
  cex.axis = axis_cex,
  lwd = axis_lwd,
  lwd.ticks = axis_lwd
)

axis(
  2,
  at = seq(0, N, by = 10),
  cex.axis = axis_cex,
  lwd = axis_lwd,
  lwd.ticks = axis_lwd,
  las = 1
)

box(lwd = axis_lwd)

lines(
  traj$time,
  traj$k,
  type = "s",
  lwd = line_lwd,
  col = col_main
)

dev.off()