# scripts/fig05-06_plot_mle_consistency_panels.R
#
# Generates the six standalone scatter plots used in Figures 5 and 6.
#
# Outputs for T = 200, 500, 1000:
#   figures/mle_consistency_unmarked_{T}.pdf
#   figures/mle_consistency_marked_{T}.pdf
#
# Required input for T = 200, 500, 1000:
#   data/comparison_estimators_{T}.csv
#   data/estimates_mle_marked_{T}.csv
#
# The comparison_estimators files contain conditional MLEs for unmarked births
# in columns b1_mle and b2_mle.
#
# The estimates_mle_marked files contain conditional MLEs for marked births
# in columns b1_mle and b2_mle.


par_true <- c(b1 = 1.01, b2 = 3.70)
Ts <- c(200, 500, 1000)

col_est      <- "#4994df"
col_est_fill <- adjustcolor(col_est, alpha.f = 0.70)
col_mean     <- "red"

pdf_family    <- "Times"
pdf_pointsize <- 12

axis_cex <- 1.00
lab_cex  <- 1.08
axis_lwd <- 0.8

pt_cex    <- 0.80
truth_cex <- 1.9
mean_cex  <- 2.2


dfs_unmarked <- lapply(Ts, function(TT) {
  read.csv(sprintf("data/comparison_estimators_%d.csv", TT))
})

dfs_marked <- lapply(Ts, function(TT) {
  read.csv(sprintf("data/estimates_mle_marked_%d.csv", TT))
})


compute_lims <- function(dfs, truth, pad_frac = 0.07) {
  all_b1 <- c(unlist(lapply(dfs, function(d) d$b1_mle)), truth["b1"])
  all_b2 <- c(unlist(lapply(dfs, function(d) d$b2_mle)), truth["b2"])

  x_rng <- range(all_b1, na.rm = TRUE)
  y_rng <- range(all_b2, na.rm = TRUE)

  x_pad <- diff(x_rng) * pad_frac
  y_pad <- diff(y_rng) * pad_frac

  if (x_pad == 0) x_pad <- 0.05
  if (y_pad == 0) y_pad <- 0.05

  list(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  )
}


plot_panel <- function(df, lims, outfile) {
  m1 <- mean(df$b1_mle, na.rm = TRUE)
  m2 <- mean(df$b2_mle, na.rm = TRUE)

  pdf(
    outfile,
    width = 2.55,
    height = 2.35,
    family = pdf_family,
    pointsize = pdf_pointsize,
    useDingbats = FALSE
  )

  par(
    mar = c(3.25, 3.9, 0.35, 0.35),
    mgp = c(2.0, 0.55, 0),
    tcl = -0.23,
    las = 1
  )

  plot(
    df$b1_mle, df$b2_mle,
    xlim = lims$xlim,
    ylim = lims$ylim,
    xlab = expression(b[1]),
    ylab = expression(b[2]),
    cex.lab = lab_cex,
    cex.axis = axis_cex,
    pch = 16,
    cex = pt_cex,
    col = col_est_fill,
    axes = FALSE
  )

  axis(1, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd)
  axis(2, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd, las = 1)
  box(lwd = axis_lwd)

  abline(v = par_true["b1"], lty = 2, lwd = 0.9, col = "grey35")
  abline(h = par_true["b2"], lty = 2, lwd = 0.9, col = "grey35")

  points(
    par_true["b1"], par_true["b2"],
    pch = 4,
    cex = truth_cex,
    lwd = 1.8,
    col = "black"
  )

  points(
    m1, m2,
    pch = 1,
    cex = mean_cex,
    lwd = 1.6,
    col = col_mean
  )

  dev.off()
  cat("Saved", outfile, "\n")
}


lims_unmarked <- compute_lims(dfs_unmarked, par_true, pad_frac = 0.07)
lims_marked   <- compute_lims(dfs_marked,   par_true, pad_frac = 0.07)

for (i in seq_along(Ts)) {
  plot_panel(
    dfs_unmarked[[i]],
    lims_unmarked,
    sprintf("figures/mle_consistency_unmarked_%d.pdf", Ts[i])
  )

  plot_panel(
    dfs_marked[[i]],
    lims_marked,
    sprintf("figures/mle_consistency_marked_%d.pdf", Ts[i])
  )
}