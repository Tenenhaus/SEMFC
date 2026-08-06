


# ============================================================
# 2. Metric functions
# ============================================================

# Euclidean estimation error
l2_error <- function(theta_hat, theta_true) {
  sqrt(sum((theta_hat - theta_true)^2))
}


# Root mean squared error across parameters
parameter_rmse <- function(theta_hat, theta_true) {
  sqrt(mean((theta_hat - theta_true)^2))
}


# Mean absolute error across parameters
parameter_mae <- function(theta_hat, theta_true) {
  mean(abs(theta_hat - theta_true))
}


# Average signed bias across parameters
mean_bias <- function(theta_hat, theta_true) {
  mean(theta_hat - theta_true)
}


# Distance between One-Step and RML
os_rml_distance <- function(theta_os, theta_rml) {
  sqrt(sum((theta_os - theta_rml)^2))
}


# Root-n scaled distance between One-Step and RML
scaled_os_rml_distance <- function(theta_os, theta_rml, n) {
  sqrt(n) * os_rml_distance(theta_os, theta_rml)
}


# Objective-function gap between One-Step and RML
likelihood_gap <- function(theta_os, theta_rml, S, objective_fun) {
  objective_fun(theta_os, S) -
    objective_fun(theta_rml, S)
}


# Maximum absolute constraint residual
constraint_residual <- function(theta_hat, constraint_fun) {

  residuals <- constraint_fun(theta_hat)

  if (length(residuals) == 0L) {
    return(0)
  }

  max(abs(residuals))
}


# Relative distance to RML
relative_os_rml_distance <- function(theta_os, theta_rml) {

  denominator <- sqrt(sum(theta_rml^2))

  if (denominator < .Machine$double.eps) {
    return(NA_real_)
  }

  os_rml_distance(theta_os, theta_rml) / denominator
}





# ============================================================
# 3. Utility functions
# ============================================================



safe_estimation <- function(model, infer = FALSE, ...) {

  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    {

      model$fit(
        infer = infer,
        ...
      )

      theta_hat <- model$get_estimate("theta")

      list(
        model = model,
        theta = theta_hat,
        success = TRUE,
        elapsed_time =
          proc.time()[["elapsed"]] - start_time,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        model = model,
        theta = NULL,
        success = FALSE,
        elapsed_time =
          proc.time()[["elapsed"]] - start_time,
        error_message = conditionMessage(e)
      )
    }
  )

  result
}





# ============================================================
# 6. Summary by sample size
# ============================================================



summarize_one_n <- function(data_n) {

  data.frame(
    n = unique(data_n$n),
    replications = nrow(data_n),

    svd_failure_rate = mean(!data_n$svd_success),
    os_failure_rate = mean(!data_n$os_success),
    rml_failure_rate = mean(!data_n$rml_success),

    svd_rmse_mean = mean_na(data_n$svd_rmse),
    svd_rmse_sd = sd_na(data_n$svd_rmse),

    os_rmse_mean = mean_na(data_n$os_rmse),
    os_rmse_sd = sd_na(data_n$os_rmse),

    rml_rmse_mean = mean_na(data_n$rml_rmse),
    rml_rmse_sd = sd_na(data_n$rml_rmse),

    scaled_distance_mean =
      mean_na(data_n$scaled_os_rml_distance),

    scaled_distance_sd =
      sd_na(data_n$scaled_os_rml_distance),

    likelihood_gap_mean =
      mean_na(data_n$likelihood_gap),

    likelihood_gap_sd =
      sd_na(data_n$likelihood_gap),

    scaled_distance_svd_rml_mean =
      mean_na(data_n$scaled_svd_rml_distance),
    scaled_distance_svd_rml_sd =
      sd_na(data_n$scaled_svd_rml_distance),
    likelihood_gap_svd_mean =
      mean_na(data_n$likelihood_gap_svd),
    likelihood_gap_svd_sd =
      sd_na(data_n$likelihood_gap_svd),

    svd_time_mean = mean_na(data_n$svd_time),
    os_time_mean = mean_na(data_n$os_time),
    rml_time_mean = mean_na(data_n$rml_time),

    os_constraint_mean =
      mean_na(data_n$os_constraint),

    rml_constraint_mean =
      mean_na(data_n$rml_constraint)
  )
}




mean_na <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}


sd_na <- function(x) {
  if (sum(!is.na(x)) <= 1L) {
    return(NA_real_)
  }

  sd(x, na.rm = TRUE)
}


# ============================================================
# 7. Componentwise Monte-Carlo bias and RMSE
# ============================================================



componentwise_metrics <- function(estimates, theta_true, n_grid) {

  output <- vector("list", length(n_grid))

  for (i in seq_along(n_grid)) {

    theta_matrix <- estimates[i, , , drop = FALSE]
    theta_matrix <- matrix(
      theta_matrix,
      ncol = length(theta_true)
    )

    valid_rows <- apply(
      theta_matrix,
      MARGIN = 1L,
      FUN = function(x) all(is.finite(x))
    )

    theta_matrix <- theta_matrix[valid_rows, , drop = FALSE]

    if (nrow(theta_matrix) == 0L) {
      next
    }

    bias <- colMeans(
      sweep(theta_matrix, 2L, theta_true, FUN = "-")
    )

    rmse <- sqrt(
      colMeans(
        sweep(theta_matrix, 2L, theta_true, FUN = "-")^2
      )
    )

    output[[i]] <- data.frame(
      n = n_grid[i],
      parameter = if (is.null(names(theta_true))) {
        paste0("theta_", seq_along(theta_true))
      } else {
        names(theta_true)
      },
      true_value = theta_true,
      bias = bias,
      rmse = rmse,
      successful_replications = nrow(theta_matrix)
    )
  }

  do.call(rbind, output)
}