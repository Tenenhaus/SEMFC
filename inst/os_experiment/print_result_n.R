library(dplyr)
library(ggplot2)

# ============================================================
# Load results
# ============================================================

x <- readRDS("inst/os_experiment/monte_carlo_scaling_n_full.rds")

results_mc <- x$results_mc


# ============================================================
# Summary by n
# ============================================================

summary_n <- results_mc %>%
  group_by(n) %>%
  summarise(
    M = n(),

    svd_rmse = sqrt(mean(svd_mse, na.rm = TRUE)),
    os_rmse  = sqrt(mean(os_mse,  na.rm = TRUE)),
    rml_rmse = sqrt(mean(rml_mse, na.rm = TRUE)),

    D_os_mean = mean(scaled_os_rml_distance, na.rm = TRUE),
    D_os_sd   = sd(scaled_os_rml_distance, na.rm = TRUE),

    D_svd_mean = mean(scaled_svd_rml_distance, na.rm = TRUE),
    D_svd_sd   = sd(scaled_svd_rml_distance, na.rm = TRUE),

    os_likelihood_gap =
      mean(likelihood_gap, na.rm = TRUE),

    os_constraint =
      max(os_constraint, na.rm = TRUE),

    rml_constraint =
      max(rml_constraint, na.rm = TRUE),

    # variance_error_os_mean = mean(variance_error_os, na.rm = TRUE),
    # variance_error_os_sd = sd(variance_error_os, na.rm = TRUE),
    # variance_error_rml_mean = mean(variance_error_rml, na.rm = TRUE),
    # variance_error_rml_sd = sd(variance_error_rml, na.rm = TRUE),

    svd_failures = sum(!svd_success),
    os_failures  = sum(!os_success),
    rml_failures = sum(!rml_success),

    .groups = "drop"
  ) %>%

  mutate(
    # variance_error_os_se = variance_error_os_sd / sqrt(M),
    # variance_error_os_lower = variance_error_os_mean - 1.96 * variance_error_os_se,
    # variance_error_os_upper = variance_error_os_mean + 1.96 * variance_error_os_se,
    #
    # variance_error_rml_se = variance_error_rml_sd / sqrt(M),
    # variance_error_rml_lower = variance_error_rml_mean - 1.96 * variance_error_rml_se,
    # variance_error_rml_upper = variance_error_rml_mean + 1.96 * variance_error_rml_se,

    D_os_se = D_os_sd / sqrt(M),
    D_os_lower = D_os_mean - 1.96 * D_os_se,
    D_os_upper = D_os_mean + 1.96 * D_os_se,

    D_svd_se = D_svd_sd / sqrt(M),
    D_svd_lower = D_svd_mean - 1.96 * D_svd_se,
    D_svd_upper = D_svd_mean + 1.96 * D_svd_se
  )


covariance_mc_summary <- data.frame(
  n = x$n_grid,
  variance_mc_error_os = x$variance_mc_error_os,
  variance_mc_error_rml = x$variance_mc_error_rml
)

summary_n <- summary_n %>%
  left_join(covariance_mc_summary, by = "n")
print(summary_n)


# ============================================================
# Main-paper table
# ============================================================

table_n <- summary_n %>%
  select(
    n,
    svd_rmse,
    os_rmse,
    rml_rmse,
    D_os_mean,
    D_svd_mean,
    variance_error_os_mean,
    variance_error_rml_mean,
    variance_mc_error_os,
    variance_mc_error_rml
  )

print(table_n)


# LaTeX rows
for (i in seq_len(nrow(table_n))) {
  cat(sprintf(
    "%d & %.4f & %.4f & %.4f & %.3f & %.3f \\\\\n",
    table_n$n[i],
    table_n$svd_rmse[i],
    table_n$os_rmse[i],
    table_n$rml_rmse[i],
    table_n$D_os_mean[i],
    table_n$D_svd_mean[i]
  ))
}


# ============================================================
# Figure data
# ============================================================

plot_n <- bind_rows(
  summary_n %>%
    transmute(
      n,
      estimator = "One-Step",
      mean = D_os_mean,
      lower = D_os_lower,
      upper = D_os_upper
    ),

  summary_n %>%
    transmute(
      n,
      estimator = "SVD",
      mean = D_svd_mean,
      lower = D_svd_lower,
      upper = D_svd_upper
    )
)

plot_n$estimator <- factor(
  plot_n$estimator,
  levels = c("One-Step", "SVD")
)



# ============================================================
# Useful diagnostics for the text
# ============================================================

summary_n %>%
  select(
    n,
    os_likelihood_gap,
    os_constraint,
    rml_constraint,
    svd_failures,
    os_failures,
    rml_failures
  ) %>%
  print()



# ============================================================
# 1. RMSE : estimation, MCSE et intervalle Monte Carlo
# ============================================================

rmse_svd <- table_rmse(estimates_svd, theta_true, n_grid)
rmse_os  <- table_rmse(estimates_os,  theta_true, n_grid)
rmse_rml <- table_rmse(estimates_rml, theta_true, n_grid)


# ============================================================
# 2. Distances normalisées : sqrt(n) * ||theta_a - theta_b||
# Moyenne, MCSE et intervalle Monte Carlo
# ============================================================

distance_os_rml <- table_distance(
  estimates_os, estimates_rml, n_grid,
  scaled = TRUE
)

distance_svd_rml <- table_distance(
  estimates_svd, estimates_rml, n_grid,
  scaled = TRUE
)


# ============================================================
# 3. Erreurs relatives des covariances plug-in
# Moyenne, MCSE et intervalle Monte Carlo
# ============================================================

covariance_plugin_os <- table_covariance_plugin(
  Vhat_os, V_0, n_grid
)

covariance_plugin_rml <- table_covariance_plugin(
  Vhat_rml, V_0, n_grid
)


# ============================================================
# 4. Covariances Monte Carlo : sans intervalle
# ============================================================

cov_mc_os <- table_covariance_mc(
  estimates_os, V_0, n_grid
)

cov_mc_rml <- table_covariance_mc(
  estimates_rml, V_0, n_grid
)

V_mc_os  <- cov_mc_os$V_mc
V_mc_rml <- cov_mc_rml$V_mc

table_variance_mc <- data.frame(
  n = n_grid,
  OS = cov_mc_os$table$variance_mc_error,
  RML = cov_mc_rml$table$variance_mc_error,
  n_valid_os = cov_mc_os$table$n_valid,
  n_valid_rml = cov_mc_rml$table$n_valid
)


# ============================================================
# 5. Coverage par paramètre
# Coverage, MCSE et intervalle Monte Carlo de Wilson
# ============================================================

coverage_os <- table_coverage(
  estimates_os, Vhat_os, theta_true, n_grid
)

coverage_rml <- table_coverage(
  estimates_rml, Vhat_rml, theta_true, n_grid
)


# ============================================================
# 6. Temps : médiane et quartiles Q25–Q75
# ============================================================

time_svd <- table_time(results_mc, "svd_time")
time_os  <- table_time(results_mc, "os_time")
time_rml <- table_time(results_mc, "rml_time")


# ============================================================
# 7. Écarts de critère
# Moyenne, MCSE et intervalle Monte Carlo
# ============================================================

gap_os <- table_scalar(
  results_mc, "likelihood_gap"
)

gap_svd <- table_scalar(
  results_mc, "likelihood_gap_svd"
)
