#!/usr/bin/env bash
set -euo pipefail

# simulate data
Rscript scripts/precompute_surviving_stats.R

# calculate estimates
Rscript scripts/compare_estimators_parallel.R      # for Figures 4, 5, 7
Rscript scripts/test_unc_mle_marked.R              # for Figure 4 (marked)
Rscript scripts/validate_mle_marked_parallel.R     # for Figure 6
Rscript scripts/boundary_test_MLE_parallel.R       # for Figure 8

# generate figures
Rscript scripts/fig03_single_trajectory.R
Rscript scripts/fig04_plot_unconditional_mle.R
Rscript scripts/fig05-06_plot_mle_consistency_panels.R
Rscript scripts/fig07_plot_compare_estimators.R
Rscript scripts/fig08_plot_boundary_testing.R

echo "All done"
