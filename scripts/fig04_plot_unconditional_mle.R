# scripts/fig04_plot_unconditional_mle.R
#
# Plot Figure 4.
#
# This script generates the two scatter-plot panels showing naive
# unconditional MLEs computed from 200 trajectories of the complete-graph
# simplicial SIS model process conditioned to survive up to T = 1000.
#
# Output:
#   figures/unconditional_mle_1000.pdf
#   figures/unconditional_mle_marked_1000.pdf
#
# Required input:
#   data/comparison_estimators_1000.csv
#   data/comparison_estimators_marked_1000.csv
#
# The first input file contains unconditional MLEs. 
# The second contains unconditional MLEs computed with observed birth marks.
#
# Both files must contain columns b1_unc and b2_unc,
# in the scaled coordinates b1 = N beta1 and b2 = N^2 beta2.
#
# TODO: add the script that generates these two input files.
# 

par_true <- c(b1 = 1.01, b2 = 3.70, mu = 1.0)

df_unmarked <- read.csv("data/comparison_estimators_1000.csv")
df_marked   <- read.csv("data/comparison_estimators_marked_1000.csv")

col_est      <- "#4994df"
col_est_fill <- adjustcolor(col_est, alpha.f = 0.70)
col_mean     <- "red"

pdf_family    <- "Times"
pdf_pointsize <- 12

axis_cex <- 1.00
lab_cex  <- 1.10
axis_lwd <- 0.8

pt_cex    <- 0.95
truth_cex <- 2.0
mean_cex  <- 2.4


compute_limits <- function(df, truth, pad_frac = 0.06) {
  x_rng <- range(c(truth["b1"], df$b1_unc), finite = TRUE)
  y_rng <- range(c(truth["b2"], df$b2_unc), finite = TRUE)

  x_pad <- pad_frac * diff(x_rng)
  y_pad <- pad_frac * diff(y_rng)

  if (x_pad == 0) x_pad <- 0.05
  if (y_pad == 0) y_pad <- 0.05

  list(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  )
}


plot_unconditional_cloud <- function(df, outfile, xlim, ylim) {
  mean_est <- c(
    b1 = mean(df$b1_unc, na.rm = TRUE),
    b2 = mean(df$b2_unc, na.rm = TRUE)
  )

  pdf(
    outfile,
    width = 3.8,
    height = 3.3,
    family = pdf_family,
    pointsize = pdf_pointsize,
    useDingbats = FALSE
  )

  par(
    mar = c(3.6, 4.2, 0.5, 0.5),
    mgp = c(2.2, 0.7, 0),
    tcl = -0.25,
    las = 1
  )

  plot(
    df$b1_unc,
    df$b2_unc,
    xlim = xlim,
    ylim = ylim,
    xlab = expression(b[1]),
    ylab = expression(b[2]),
    pch = 16,
    cex = pt_cex,
    col = col_est_fill,
    axes = FALSE,
    cex.lab = lab_cex,
    cex.axis = axis_cex
  )

  axis(1, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd)
  axis(2, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd, las = 1)
  box(lwd = axis_lwd)

  abline(v = par_true["b1"], lty = 2, lwd = 0.9, col = "grey35")
  abline(h = par_true["b2"], lty = 2, lwd = 0.9, col = "grey35")

  points(
    par_true["b1"],
    par_true["b2"],
    pch = 4,
    cex = truth_cex,
    lwd = 1.8,
    col = "black"
  )

  points(
    mean_est["b1"],
    mean_est["b2"],
    pch = 1,
    cex = mean_cex,
    lwd = 1.6,
    col = col_mean
  )

  dev.off()
}


lims_unmarked <- compute_limits(df_unmarked, par_true)
lims_marked   <- compute_limits(df_marked, par_true)

plot_unconditional_cloud(
  df_unmarked,
  outfile = "figures/unconditional_mle_1000.pdf",
  xlim = lims_unmarked$xlim,
  ylim = lims_unmarked$ylim
)

plot_unconditional_cloud(
  df_marked,
  outfile = "figures/unconditional_mle_marked_1000.pdf",
  xlim = lims_marked$xlim,
  ylim = lims_marked$ylim
)
