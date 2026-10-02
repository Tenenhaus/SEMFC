


# ============================================================
# 2. Metric functions
# ============================================================


# Root mean squared error across parameters
parameter_mse <- function(theta_hat, theta_true) {
  mean((theta_hat - theta_true)^2)
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



os_rml_mse <- function(theta_os, theta_rml) {
  mean((theta_os - theta_rml)^2)
}


variance_plug_in_error <- function(Vhat, V_0){
  norm(Vhat - V_0, type = "F") / norm(V_0, type = "F")
}



# ============================================================
# 3. Utility functions
# ============================================================



safe_estimation <- function(model, compute_cov = T, infer = FALSE, ...) {

  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    {

      model$fit(
        infer = infer,
        ...
      )
      elapsed_time <- proc.time()[["elapsed"]] - start_time
      theta_hat <- model$get_estimate("theta")
      if (compute_cov) {
        Vhat <- get_covariance_matrix_plug_in(model$get_estimate('all'),
                                            model$get_model(), model$get_data())
      } else {
        Vhat <- NULL
      }

      list(
        model = model,
        theta = theta_hat,
        V = Vhat,
        success = TRUE,
        elapsed_time = elapsed_time,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        model = model,
        theta = NULL,
        V = NULL,
        success = FALSE,
        elapsed_time =
          proc.time()[["elapsed"]] - start_time,
        error_message = conditionMessage(e)
      )
    }
  )

  result
}

safe_estimation_lavaan <- function(
    model,
    data,
    composites_cov = "fixed"
) {

  start_time <- proc.time()[["elapsed"]]

  fit <- tryCatch(
    lavaan::sem(
      model = model,
      data = data,
      estimator = "ML",
      likelihood = "wishart",
      composites.cov = composites_cov,
      se = "none",
      test = "none"
    ),
    error = function(e) e
  )

  elapsed_time <- proc.time()[["elapsed"]] - start_time

  if (inherits(fit, "error")) {
    return(
      list(
        success = fit$optim$converged,
        elapsed_time = elapsed_time,
        fit = NULL,
        error_message = conditionMessage(fit)
      )
    )
  }

  list(
    success = TRUE,
    elapsed_time = elapsed_time,
    fit = fit,
    error_message = NA_character_
  )
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

    svd_rmse = sqrt(mean_na(data_n$svd_mse)),

    os_rmse = sqrt(mean_na(data_n$os_mse)),

    rml_rmse = sqrt(mean_na(data_n$rml_mse)),

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


summarize_one_q <- function(data_q) {

  data.frame(
    q = unique(data_q$q),
    p = unique(data_q$p),
    d = unique(data_q$d),
    n = unique(data_q$n),
    replications = nrow(data_q),

    svd_failure_rate = mean(!data_q$svd_success),
    os_failure_rate = mean(!data_q$os_success),
    rml_failure_rate = mean(!data_q$rml_success),

    svd_rmse = sqrt(mean_na(data_q$svd_mse)),

    os_rmse = sqrt(mean_na(data_q$os_mse)),

    rml_rmse = sqrt(mean_na(data_q$rml_mse)),

    rmse_os_rml =
      sqrt(mean_na(data_q$os_rml_mse)),


    os_rml_distance_mean =
      mean_na(data_q$os_rml_distance),

    relative_os_rml_distance_mean =
      mean_na(data_q$relative_os_rml_distance),

    scaled_distance_mean =
      mean_na(data_q$scaled_os_rml_distance),

    likelihood_gap_mean =
      mean_na(data_q$likelihood_gap),

    svd_time_mean =
      mean_na(data_q$svd_time),

    os_time_mean =
      mean_na(data_q$os_time),

    rml_time_mean =
      mean_na(data_q$rml_time),

    lavaan_time_mean =
      mean_na(data_q$lavaan_time),

    speedup =
      mean_na(data_q$rml_time) /
      mean_na(data_q$os_time),

    speedup_lavaan =
      mean_na(data_q$lavaan_time) /
      mean_na(data_q$os_time),

    svd_constraint_mean =
      mean_na(data_q$svd_constraint),

    os_constraint_mean =
      mean_na(data_q$os_constraint),

    rml_constraint_mean =
      mean_na(data_q$rml_constraint)
  )
}

componentwise_metrics_indicators <- function(
    estimates,
    theta_true_list,
    q_grid
) {

  output <- vector("list", length(q_grid))

  for (i in seq_along(q_grid)) {

    theta_matrix <- estimates[[i]]
    theta_true <- theta_true_list[[i]]

    if (is.null(theta_matrix) || is.null(theta_true)) {
      next
    }

    valid_rows <- apply(
      theta_matrix,
      MARGIN = 1L,
      FUN = function(x) all(is.finite(x))
    )

    theta_matrix <-
      theta_matrix[valid_rows, , drop = FALSE]

    if (nrow(theta_matrix) == 0L) {
      next
    }

    errors <- sweep(
      theta_matrix,
      MARGIN = 2L,
      STATS = theta_true,
      FUN = "-"
    )

    bias <- colMeans(errors)

    rmse <- sqrt(
      colMeans(errors^2)
    )

    output[[i]] <- data.frame(
      q = q_grid[i],
      p = 6 * q_grid[i],
      d = length(theta_true),

      parameter =
        if (is.null(names(theta_true))) {
          paste0(
            "theta_",
            seq_along(theta_true)
          )
        } else {
          names(theta_true)
        },

      true_value = theta_true,
      bias = bias,
      rmse = rmse,

      successful_replications =
        nrow(theta_matrix)
    )
  }

  do.call(rbind, output)
}


summarize_one_J <- function(data_J) {

  data.frame(
    J = unique(data_J$J),
    q = unique(data_J$q),
    p = unique(data_J$p),
    d = unique(data_J$d),
    n = unique(data_J$n),
    replications = nrow(data_J),

    svd_failure_rate = mean(!data_J$svd_success),
    os_failure_rate  = mean(!data_J$os_success),
    rml_failure_rate = mean(!data_J$rml_success),

    svd_rmse = sqrt(mean_na(data_J$svd_mse)),

    os_rmse = sqrt(mean_na(data_J$os_mse)),

    rml_rmse = sqrt(mean_na(data_J$rml_mse)),

    os_rml_rmse =
      sqrt(mean_na(data_J$os_rml_mse)),


    likelihood_gap_mean =
      mean_na(data_J$likelihood_gap),

    likelihood_gap_sd =
      sd_na(data_J$likelihood_gap),

    svd_time_mean =
      mean_na(data_J$svd_time),

    os_time_mean =
      mean_na(data_J$os_time),

    rml_time_mean =
      mean_na(data_J$rml_time),

    speedup =
      mean_na(data_J$rml_time) /
      mean_na(data_J$os_time),

    svd_constraint_mean =
      mean_na(data_J$svd_constraint),

    os_constraint_mean =
      mean_na(data_J$os_constraint),

    rml_constraint_mean =
      mean_na(data_J$rml_constraint)
  )
}

componentwise_metrics_J <- function(
    estimates,
    theta_true_list,
    J_grid,
    q = 3
) {

  output <- vector("list", length(J_grid))

  for (i in seq_along(J_grid)) {

    theta_matrix <- estimates[[i]]
    theta_true <- theta_true_list[[i]]

    if (is.null(theta_matrix) || is.null(theta_true)) {
      next
    }

    valid_rows <- apply(
      theta_matrix,
      1L,
      function(x) all(is.finite(x))
    )

    theta_matrix <-
      theta_matrix[valid_rows, , drop = FALSE]

    if (nrow(theta_matrix) == 0L) {
      next
    }

    errors <- sweep(
      theta_matrix,
      2L,
      theta_true,
      FUN = "-"
    )

    bias <- colMeans(errors)

    rmse <- sqrt(
      colMeans(errors^2)
    )

    output[[i]] <- data.frame(
      J = J_grid[i],
      q = q,
      p = J_grid[i] * q,
      d = length(theta_true),

      parameter =
        if (is.null(names(theta_true))) {
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


get_covariance_matrix_plug_in <- function(init_estimate, model, data) {

  lambda <- init_estimate$lambda
  S_diag_composites <- data$S_diag_composites
  block_sizes <- model$block_sizes
  lengths_theta <- model$lengths_theta
  mode <- model$mode
  BETA <- init_estimate$beta
  GAMMA <- init_estimate$gamma
  sigma_implied <- init_estimate$sigma_implied
  r <- sum(mode == "formative")

  jac_Sigma_theta <- compute_J_Sigma_theta(block_sizes,
                                         lengths_theta,
                                         mode, model$dag,
                                         init_estimate$p_implied,
                                         model$relation_matrix,
                                         lambda, BETA, GAMMA)
  H <- compute_Hessian(sigma_implied, jac_Sigma_theta)
  J_h_eq <- compute_gradient_constraint(lambda,
                                      S_diag_composites,
                                      block_sizes,
                                      lengths_theta,
                                      mode)
  ktt_mat <- rbind(cbind(H, t(J_h_eq)),
                   cbind(J_h_eq, matrix(0, r, r)))
  return (2*(solve(ktt_mat)[1:sum(lengths_theta), 1:sum(lengths_theta)]))
}



coverage <- function(estimate, theta_true, se, alpha = 0.05) {
  lower <- estimate$theta - qnorm(1 - alpha / 2) * sqrt(diag(estimate$V))
  upper <- estimate$theta + qnorm(1 - alpha / 2) * sqrt(diag(estimate$V))
  coverage <- (theta_true >= lower) & (theta_true <= upper)
  return(coverage)
}


coverage_table <- function(estimates, Vhat, theta_true, n_grid) {
  dims <- dim(estimates)
  d <- dims[3]

  # Extraire les diagonales des matrices de covariance :
  # résultat de dimensions (n, réplication, paramètre)
  var_diag <- array(
    matrix(Vhat, nrow = dims[1] * dims[2])[
      , seq.int(1L, d * d, by = d + 1L),
      drop = FALSE
    ],
    dim = dims,
    dimnames = dimnames(estimates)
  )
  # Variances invalides ou nulles : pas d'IC exploitable
  var_diag[!is.finite(var_diag) | var_diag <= 0] <- NA_real_

  # SE = sqrt(V_jj / n)
  se <- sqrt(sweep(var_diag, MARGIN = 1, STATS = n_grid, FUN = "/"))

  # Erreurs d'estimation, par rapport au vrai paramètre
  error <- sweep(
    estimates,
    MARGIN = 3,
    STATS = as.numeric(theta_true),
    FUN = "-"
  )
  error[!is.finite(error)] <- NA_real_

  # TRUE si le vrai paramètre appartient à l'IC à 95 %
  covered <- abs(error) <= qnorm(0.975) * se

  # Agrégation sur les réplications (dimension 2)
  n_valid <- apply(!is.na(covered), c(1, 3), sum)

  coverage <- apply(covered, c(1, 3), sum, na.rm = TRUE) /
    n_valid

  coverage[n_valid == 0] <- NA_real_

  # Erreur-type Monte Carlo du coverage
  mcse <- sqrt(coverage * (1 - coverage) / n_valid)

  parameter_names <- dimnames(estimates)[[3]]

  if (is.null(parameter_names)) {
    parameter_names <- paste0("theta_", seq_len(d))
  }

  coverage_table <- data.frame(
    n         = rep(n_grid, times = d),
    parameter = rep(parameter_names, each = length(n_grid)),
    coverage  = as.vector(coverage),
    n_valid   = as.vector(n_valid),
    mcse      = as.vector(mcse)
  )
  return(coverage_table)
}
