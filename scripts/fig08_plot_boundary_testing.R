# scripts/fig08_plot_boundary_testing.R
#
# Plot Figure 8.
#
# Input:
#   data/boundary_test_MLE.csv
#
# Output:
#   figures/boundary_testing.pdf
#
# The input file is produced by scripts/boundary_test_MLE_parallel.R and
# contains the standardized statistics Z_b2 from 1000 survival-conditioned
# trajectories generated under H0: b2 = 0.


h0_results <- read.csv("data/boundary_test_MLE.csv")

cat("Mean of empirical se(b2_hat):", mean(h0_results$se_b2, na.rm = TRUE), "\n")
cat("Empirical SD of b2_hat:      ", sd(h0_results$b2_hat, na.rm = TRUE), "\n")
cat(
  "Ratio:",
  mean(h0_results$se_b2, na.rm = TRUE) /
    sd(h0_results$b2_hat, na.rm = TRUE),
  "\n"
)

col_hist <- "grey85"
col_ref  <- "#4994df"

pdf_family    <- "Times"
pdf_pointsize <- 10

axis_cex <- 0.88
lab_cex  <- 0.95
axis_lwd <- 0.75


z_vals <- h0_results$Z_b2
z_vals <- z_vals[is.finite(z_vals)]

z_max <- max(abs(z_vals), na.rm = TRUE)

bar_width <- 2 * z_max / 15
n_half <- ceiling(z_max / bar_width + 0.5)
breaks <- seq(-(n_half - 0.5), (n_half - 0.5)) * bar_width


dir.create("figures", showWarnings = FALSE, recursive = TRUE)

pdf(
  "figures/boundary_testing.pdf",
  width = 3.15,
  height = 2.45,
  family = pdf_family,
  pointsize = pdf_pointsize,
  useDingbats = FALSE
)

par(
  mar = c(3.0, 3.5, 0.25, 0.25),
  mgp = c(1.85, 0.55, 0),
  tcl = -0.22,
  las = 1
)

hist(
  z_vals,
  breaks = breaks,
  freq = FALSE,
  xlab = expression(Z[list(2, T)]),
  ylab = "Density",
  main = "",
  col = col_hist,
  border = "white",
  axes = FALSE,
  cex.lab = lab_cex,
  cex.axis = axis_cex
)

axis(1, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd)
axis(2, cex.axis = axis_cex, lwd = axis_lwd, lwd.ticks = axis_lwd, las = 1)
box(lwd = axis_lwd)

curve(
  dnorm(x, mean = 0, sd = 1),
  add = TRUE,
  col = col_ref,
  lwd = 1.05
)

abline(v = 0, col = "grey35", lty = 2, lwd = 0.8)

dev.off()