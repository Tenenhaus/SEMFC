library(dplyr)
library(ggplot2)

# ============================================================
# Load results
# ============================================================

x <- readRDS("inst/os_experiment/monte_carlo_scaling_J_full.rds")
results_mc <- if (is.data.frame(x)) x else x$results_mc


# ============================================================
# Summary by J
# ============================================================

summary_J <- results_mc %>%
  group_by(J) %>%
  summarise(
    d = unique(d),

    svd_rmse = sqrt(mean(svd_mse, na.rm = TRUE)),
    os_rmse  = sqrt(mean(os_mse,  na.rm = TRUE)),
    rml_rmse = sqrt(mean(rml_mse, na.rm = TRUE)),

    os_time = median(os_time, na.rm = TRUE),
    os_q25  = quantile(os_time, 0.25, na.rm = TRUE),
    os_q75  = quantile(os_time, 0.75, na.rm = TRUE),

    rml_time = median(rml_time, na.rm = TRUE),
    rml_q25  = quantile(rml_time, 0.25, na.rm = TRUE),
    rml_q75  = quantile(rml_time, 0.75, na.rm = TRUE),

    speedup = median(
      ifelse(os_time > 0, rml_time / os_time, NA_real_),
      na.rm = TRUE
    ),

    svd_failures = sum(!svd_success),
    os_failures  = sum(!os_success),
    rml_failures = sum(!rml_success),

    os_constraint  = max(os_constraint,  na.rm = TRUE),
    rml_constraint = max(rml_constraint, na.rm = TRUE),

    .groups = "drop"
  )


print(summary_J)


# ============================================================
# Main-paper table
# ============================================================

table_J <- summary_J %>%
  select(
    J,
    d,
    svd_rmse,
    os_rmse,
    rml_rmse,
    os_time,
    rml_time,
    speedup
  )

print(table_J)


# LaTeX rows
for (i in seq_len(nrow(table_J))) {
  cat(sprintf(
    "%d & %d & %.4f & %.4f & %.4f & %.2f & %.2f & %.1f \\\\\n",
    table_J$J[i],
    table_J$d[i],
    table_J$svd_rmse[i],
    table_J$os_rmse[i],
    table_J$rml_rmse[i],
    table_J$os_time[i],
    table_J$rml_time[i],
    table_J$speedup[i]
  ))
}


# ============================================================
# Runtime figure data
# ============================================================

plot_J <- bind_rows(
  summary_J %>%
    transmute(
      J,
      estimator = "One-Step",
      median = os_time,
      lower = os_q25,
      upper = os_q75
    ),

  summary_J %>%
    transmute(
      J,
      estimator = "RML",
      median = rml_time,
      lower = rml_q25,
      upper = rml_q75
    )
)

plot_J$estimator <- factor(
  plot_J$estimator,
  levels = c("One-Step", "RML")
)


# ============================================================
# Camera-ready runtime figure
# ============================================================

p <- ggplot(
  plot_J,
  aes(
    x = J,
    y = median,
    linetype = estimator,
    group = estimator
  )
) +
  geom_line(linewidth = 0.75) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.5,
    linewidth = 0.45
  ) +
  scale_linetype_manual(
    values = c(
      "One-Step" = "solid",
      "RML" = "dashed"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_J$J))
  ) +
  scale_y_log10() +
  labs(
    x = "Number of blocks J",
    y = "Runtime (seconds)",
    linetype = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  "inst/os_experiment/figure_scaling_J_runtime.pdf",
  p,
  device = cairo_pdf,
  width = 4.6,
  height = 3.0,
  units = "in"
)


# ============================================================
# Diagnostics
# ============================================================

summary_J %>%
  select(
    J,
    svd_failures,
    os_failures,
    rml_failures,
    os_constraint,
    rml_constraint
  ) %>%
  print()