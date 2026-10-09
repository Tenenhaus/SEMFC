library(ggplot2)

# ============================================================
# 1. Load saved experiment
# ============================================================

save_path <- "inst/os_experiment/monte_carlo_scaling_J_full.rds"

mc_J <- readRDS(save_path)

# Compatible with either:
# saveRDS(results_mc, ...)
# or
# saveRDS(list(results_mc = results_mc, ...), ...)
if (is.data.frame(mc_J)) {
  
  results_mc_J <- mc_J
  
} else if (!is.null(mc_J$results_mc)) {
  
  results_mc_J <- mc_J$results_mc
  
} else {
  
  stop("Could not find results_mc in the saved RDS object.")
}

cat("Number of rows:", nrow(results_mc_J), "\n")
cat("J values:", sort(unique(results_mc_J$J)), "\n")

print(names(results_mc_J))


# ============================================================
# 2. Helper functions
# ============================================================

finite_values <- function(x) {
  x[is.finite(x)]
}

rmse_from_mse <- function(x) {
  x <- finite_values(x)
  sqrt(mean(x))
}

median_na <- function(x) {
  median(finite_values(x))
}

q25_na <- function(x) {
  quantile(
    finite_values(x),
    probs = 0.25,
    names = FALSE
  )
}

q75_na <- function(x) {
  quantile(
    finite_values(x),
    probs = 0.75,
    names = FALSE
  )
}


# ============================================================
# 3. Summary for one value of J
# ============================================================

summarize_one_J_final <- function(data_J) {
  
  # Paired runtime ratio:
  # computed replication by replication
  speedup <- data_J$rml_time / data_J$os_time
  
  data.frame(
    
    J = unique(data_J$J),
    q = unique(data_J$q),
    p = unique(data_J$p),
    d = unique(data_J$d),
    n = unique(data_J$n),
    
    replications = nrow(data_J),
    
    # --------------------------------------------------------
    # Failure rates
    # --------------------------------------------------------
    
    svd_failure_rate =
      mean(!data_J$svd_success, na.rm = TRUE),
    
    os_failure_rate =
      mean(!data_J$os_success, na.rm = TRUE),
    
    rml_failure_rate =
      mean(!data_J$rml_success, na.rm = TRUE),
    
    # --------------------------------------------------------
    # Statistical accuracy
    # --------------------------------------------------------
    
    svd_rmse =
      rmse_from_mse(data_J$svd_mse),
    
    os_rmse =
      rmse_from_mse(data_J$os_mse),
    
    rml_rmse =
      rmse_from_mse(data_J$rml_mse),
    
    # --------------------------------------------------------
    # OS-RML discrepancy
    # --------------------------------------------------------
    
    os_rml_rmse =
      rmse_from_mse(data_J$os_rml_mse),
    
    likelihood_gap =
      mean(
        finite_values(data_J$likelihood_gap)
      ),
    
    # --------------------------------------------------------
    # SVD runtime
    # --------------------------------------------------------
    
    svd_time_median =
      median_na(data_J$svd_time),
    
    svd_time_q25 =
      q25_na(data_J$svd_time),
    
    svd_time_q75 =
      q75_na(data_J$svd_time),
    
    # --------------------------------------------------------
    # One-Step runtime
    # --------------------------------------------------------
    
    os_time_median =
      median_na(data_J$os_time),
    
    os_time_q25 =
      q25_na(data_J$os_time),
    
    os_time_q75 =
      q75_na(data_J$os_time),
    
    # --------------------------------------------------------
    # RML runtime
    # --------------------------------------------------------
    
    rml_time_median =
      median_na(data_J$rml_time),
    
    rml_time_q25 =
      q25_na(data_J$rml_time),
    
    rml_time_q75 =
      q75_na(data_J$rml_time),
    
    # --------------------------------------------------------
    # Paired speedup
    # --------------------------------------------------------
    
    speedup_median =
      median_na(speedup),
    
    speedup_q25 =
      q25_na(speedup),
    
    speedup_q75 =
      q75_na(speedup),
    
    # --------------------------------------------------------
    # Constraints
    # --------------------------------------------------------
    
    os_constraint =
      mean(
        finite_values(data_J$os_constraint)
      ),
    
    rml_constraint =
      mean(
        finite_values(data_J$rml_constraint)
      )
  )
}

ggsave(
  filename =
    "inst/os_experiment/figure_scaling_J_runtime.pdf",
  plot = fig_runtime_J,
  width = 4.6,
  height = 3.0,
  device = cairo_pdf
)
