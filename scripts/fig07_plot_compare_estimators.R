# scripts/fig07_plot_compare_estimators.R
#
# Plot Figure 7.
#
# Generates two standalone panels for the sample means of three estimators
# as functions of T:
#
#   figures/comparison_estimator_means_a.pdf   # panel (a): b1
#   figures/comparison_estimator_means_b.pdf   # panel (b): b2
#
# Required input for T = 100, 200, ..., 1000:
#   data/comparison_estimators_{T}.csv
#
# Each input file contains estimates from M = 200 surviving trajectories:
#   b1_unc,  b2_unc,  mu_unc
#   b1_mle,  b2_mle,  mu_mle
#   b1_qmle, b2_qmle, mu_qmle


par_true <- c(b1 = 1.01, b2 = 3.70, mu = 1.0)

T_grid   <- seq(100, 1000, by = 100)
T_labels <- c(100, 300, 500, 700, 1000)


col_unc  <- "#cb5b6b"
col_mle  <- "#4994df"
col_qmle <- "black"

pch_unc  <- 15
pch_mle  <- 16
pch_qmle <- 17

pdf_family    <- "Times"
pdf_pointsize <- 12

axis_cex <- 1.00
lab_cex  <- 1.10
axis_lwd <- 0.8

line_lwd <- 1.3
pt_cex   <- 0.95


summary_df <- data.frame(
  T_max = T_grid,
  mean_b1_unc  = NA_real_,
  mean_b1_mle  = NA_real_,
  mean_b1_qmle = NA_real_,
  mean_b2_unc  = NA_real_,
  mean_b2_mle  = NA_real_,
  mean_b2_qmle = NA_real_
)

for (i in seq_along(T_grid)) {
  df <- read.csv(sprintf("data/comparison_estimators_%d.csv", T_grid[i]))

  summary_df$mean_b1_unc[i]  <- mean(df$b1_unc,  na.rm = TRUE)
  summary_df$mean_b1_mle[i]  <- mean(df$b1_mle,  na.rm = TRUE)
  summary_df$mean_b1_qmle[i] <- mean(df$b1_qmle, na.rm = TRUE)

  summary_df$mean_b2_unc[i]  <- mean(df$b2_unc,  na.rm = TRUE)
  summary_df$mean_b2_mle[i]  <- mean(df$b2_mle,  na.rm = TRUE)
  summary_df$mean_b2_qmle[i] <- mean(df$b2_qmle, na.rm = TRUE)
}


compute_ylim <- function(values, truth, pad_frac = 0.08) {
  rng <- range(c(values, truth), na.rm = TRUE)
  pad <- diff(rng) * pad_frac

  if (pad == 0) {
    pad <- 0.05
  }

  c(rng[1] - pad, rng[2] + pad)
}


plot_means_panel <- function(
  outfile,
  y_unc,
  y_mle,
  y_qmle,
  y_true,
  ylab_expr,
  legend_on = FALSE,
  legend_pos = "topright"
) {
  ylim <- compute_ylim(c(y_unc, y_mle, y_qmle), y_true, pad_frac = 0.08)

  pdf(
    outfile,
    width = 3.8,
    height = 3.2,
    family = pdf_family,
    pointsize = pdf_pointsize,
    useDingbats = FALSE
  )

  par(
    mar = c(3.6, 4.1, 0.45, 0.45),
    mgp = c(2.15, 0.65, 0),
    tcl = -0.25,
    las = 1
  )

  plot(
    summary_df$T_max,
    y_unc,
    type = "b",
    pch = pch_unc,
    col = col_unc,
    lwd = line_lwd,
    cex = pt_cex,
    xlab = expression(T),
    ylab = ylab_expr,
    xlim = c(80, 1020),
    ylim = ylim,
    xaxt = "n",
    yaxt = "n",
    cex.lab = lab_cex,
    cex.axis = axis_cex,
    axes = FALSE
  )

  axis(
    1,
    at = T_labels,
    labels = T_labels,
    cex.axis = axis_cex,
    lwd = axis_lwd,
    lwd.ticks = axis_lwd
  )

  axis(
    2,
    cex.axis = axis_cex,
    lwd = axis_lwd,
    lwd.ticks = axis_lwd,
    las = 1
  )

  box(lwd = axis_lwd)

  lines(
    summary_df$T_max,
    y_mle,
    type = "b",
    pch = pch_mle,
    col = col_mle,
    lwd = line_lwd,
    cex = pt_cex
  )

  lines(
    summary_df$T_max,
    y_qmle,
    type = "b",
    pch = pch_qmle,
    col = col_qmle,
    lwd = line_lwd,
    cex = pt_cex
  )

  abline(h = y_true, lty = 2, lwd = 0.9, col = "grey35")

  if (legend_on) {
    legend(
      legend_pos,
      legend = c("uncond. MLE", "cond. MLE", "QMLE"),
      col = c(col_unc, col_mle, col_qmle),
      pch = c(pch_unc, pch_mle, pch_qmle),
      lty = 1,
      lwd = line_lwd,
      pt.cex = pt_cex,
      bty = "n",
      cex = 0.88
    )
  }

  dev.off()
  cat("Saved", outfile, "\n")
}


plot_means_panel(
  outfile    = "figures/comparison_estimator_means_a.pdf",
  y_unc      = summary_df$mean_b1_unc,
  y_mle      = summary_df$mean_b1_mle,
  y_qmle     = summary_df$mean_b1_qmle,
  y_true     = par_true["b1"],
  ylab_expr  = expression(b[1]),
  legend_on  = TRUE,
  legend_pos = "topright"
)

plot_means_panel(
  outfile    = "figures/comparison_estimator_means_b.pdf",
  y_unc      = summary_df$mean_b2_unc,
  y_mle      = summary_df$mean_b2_mle,
  y_qmle     = summary_df$mean_b2_qmle,
  y_true     = par_true["b2"],
  ylab_expr  = expression(b[2]),
  legend_on  = FALSE
)