# ============================================================
# Monte-Carlo experiment:
# variation of the number of indicators per block, with a fixed number of blocks.
# ============================================================
devtools::load_all()
rm(list = ls())
source('inst/os_experiment/utils_os.R')
source('inst/os_experiment/population.R')

set.seed(123)


# ============================================================
# 1. Simulation parameters
# ============================================================

q_grid <- c(3,5,10,15,20)
N <- 2000
n_rep  <- 2




# ============================================================
# 4. Storage objects
# ============================================================

results <- vector(
  mode = "list",
  length = length(q_grid) * n_rep
)

# Optional storage of all parameter estimates.
# Useful for componentwise bias and coverage analyses.
estimates_svd <- vector("list", length(q_grid))
estimates_os  <- vector("list", length(q_grid))
estimates_rml <- vector("list", length(q_grid))
theta_true_list <- vector("list", length(q_grid))

names(estimates_svd) <-
  names(estimates_os) <-
  names(estimates_rml) <-
  names(theta_true_list) <-
  as.character(q_grid)

# ============================================================
# 5. Monte-Carlo loop
# ============================================================

row_id <- 0L

for (q_index in seq_along(q_grid)) {

  q <- q_grid[q_index]
  source('inst/simulations/data_simulation_mixed.R')
  lavaan_model <- generate_sem_model_J(6, q)

  # True parameter vector.
  # It must follow exactly the same ordering as the estimators.
  theta_true <- true_param_with_S

  d <- length(theta_true)

  theta_true_list[[q_index]] <- theta_true

  make_storage <- function() {
    matrix(
      NA_real_,
      nrow = n_rep,
      ncol = d,
      dimnames = list(
        replication = seq_len(n_rep),
        parameter = names(theta_true)
      )
    )
  }

  estimates_svd[[q_index]] <- make_storage()
  estimates_os[[q_index]]  <- make_storage()
  estimates_rml[[q_index]] <- make_storage()

  message(
    "Number of indicators per block q = ", q,
    " (", q_index, "/", length(q_grid), ")"
  )

  for (rep_id in seq_len(n_rep)) {

    row_id <- row_id + 1L

    # --------------------------------------------------------
    # Simulate one sample
    # --------------------------------------------------------

    # Expected output:
    #   either the empirical covariance matrix S directly,
    #   or adapt this block if simulate_sample_fun returns raw data.
    source('inst/model/model_mixed.R')


    # Lavaan model estimation for comparison
    fit_lavaan <- safe_estimation_lavaan(
      model = lavaan_model,
      data = X_2,
      composites_cov = "fixed"
    )



    # --------------------------------------------------------
    # SVD-SEM
    # --------------------------------------------------------
    modelsvd <- SemFC$new(Y_2, relation_matrix = C, mode=mode, estimator = "svd")
    fit_svd <- safe_estimation(
      modelsvd
    )



    # ------------------------------------------------------
    # One-Step estimator
    # ------------------------------------------------------
    modelos <- SemFC$new(Y_2, relation_matrix = C, mode=mode, estimator = "one_step")
    fit_os <- safe_estimation(
      modelos
    )

    # ------------------------------------------------------
    # Restricted maximum likelihood
    # ------------------------------------------------------
    modelml <- SemFC$new(Y_2, relation_matrix = C, mode=mode, estimator = "ml")
    fit_rml <- safe_estimation(
      modelml
    )


    objective_fun <- function(x, S) {
      F1(x, S, modelsvd$get_model())
    }

    constraint_fun <- function(x) {
      heq(x, cov(X_2), modelsvd$get_model())
    }






    # --------------------------------------------------------
    # Store complete parameter vectors
    # --------------------------------------------------------

    if (fit_svd$success) {
      estimates_svd[[q_index]][rep_id, ] <- fit_svd$theta
    }

    if (fit_os$success) {
      estimates_os[[q_index]][rep_id, ] <- fit_os$theta
    }

    if (fit_rml$success) {
      estimates_rml[[q_index]][rep_id, ] <- fit_rml$theta
    }


    # --------------------------------------------------------
    # Individual estimator metrics
    # --------------------------------------------------------

    svd_mse <- if (fit_svd$success) {
      parameter_mse(fit_svd$theta, theta_true)
    } else {
      NA_real_
    }

    os_mse <- if (fit_os$success) {
      parameter_mse(fit_os$theta, theta_true)
    } else {
      NA_real_
    }

    rml_mse <- if (fit_rml$success) {
      parameter_mse(fit_rml$theta, theta_true)
    } else {
      NA_real_
    }

    svd_constraint <- if (fit_svd$success) {
      constraint_residual(fit_svd$theta, constraint_fun)
    } else {
      NA_real_
    }

    os_constraint <- if (fit_os$success) {
      constraint_residual(fit_os$theta, constraint_fun)
    } else {
      NA_real_
    }

    rml_constraint <- if (fit_rml$success) {
      constraint_residual(fit_rml$theta, constraint_fun)
    } else {
      NA_real_
    }


    # --------------------------------------------------------
    # Direct comparison between One-Step and RML
    # --------------------------------------------------------

    both_success <- fit_os$success && fit_rml$success

    distance_os_rml <- if (both_success) {
      os_rml_distance(fit_os$theta, fit_rml$theta)
    } else {
      NA_real_
    }

    scaled_distance_os_rml <- if (both_success) {
      scaled_os_rml_distance(
        theta_os = fit_os$theta,
        theta_rml = fit_rml$theta,
        n = N
      )
    } else {
      NA_real_
    }

    relative_distance <- if (both_success) {
      relative_os_rml_distance(
        theta_os = fit_os$theta,
        theta_rml = fit_rml$theta
      )
    } else {
      NA_real_
    }

    objective_gap <- if (both_success) {
      likelihood_gap(
        theta_os = fit_os$theta,
        theta_rml = fit_rml$theta,
        S = cov(X_2),
        objective_fun = objective_fun
      )
    } else {
      NA_real_
    }

    mse_os_rml <- if (fit_os$success && fit_rml$success) {
      os_rml_mse(
        fit_os$theta,
        fit_rml$theta
      )
    } else {
      NA_real_
    }


    # --------------------------------------------------------
    # Direct comparison between svd and RML
    # --------------------------------------------------------

    both_success <- fit_svd$success && fit_rml$success

    distance_svd_rml <- if (both_success) {
      os_rml_distance(fit_svd$theta, fit_rml$theta)
    } else {
      NA_real_
    }

    scaled_distance_svd_rml <- if (both_success) {
      scaled_os_rml_distance(
        theta_os = fit_svd$theta,
        theta_rml = fit_rml$theta,
        n = N
      )
    } else {
      NA_real_
    }

    relative_distance_svd <- if (both_success) {
      relative_os_rml_distance(
        theta_os = fit_svd$theta,
        theta_rml = fit_rml$theta
      )
    } else {
      NA_real_
    }

    objective_gap_svd <- if (both_success) {
      likelihood_gap(
        theta_os = fit_svd$theta,
        theta_rml = fit_rml$theta,
        S = cov(X_2),
        objective_fun = objective_fun
      )
    } else {
      NA_real_
    }





    # --------------------------------------------------------
    # Store one row per Monte-Carlo replication
    # --------------------------------------------------------

    results[[row_id]] <- data.frame(
      q = q,
      p = 6 * q,
      d = d,
      n = N,
      replication = rep_id,

      svd_success = fit_svd$success,
      os_success = fit_os$success,
      rml_success = fit_rml$success,
      lavaan_success = fit_lavaan$success,

      svd_mse = svd_mse,
      os_mse = os_mse,
      rml_mse = rml_mse,

      svd_constraint = svd_constraint,
      os_constraint = os_constraint,
      rml_constraint = rml_constraint,

      os_rml_distance = distance_os_rml,
      scaled_os_rml_distance = scaled_distance_os_rml,
      relative_os_rml_distance = relative_distance,
      likelihood_gap = objective_gap,

      svd_rml_distance = distance_svd_rml,
      scaled_svd_rml_distance = scaled_distance_svd_rml,
      relative_svd_rml_distance = relative_distance_svd,
      likelihood_gap_svd = objective_gap_svd,
      os_rml_mse = mse_os_rml,

      svd_time = fit_svd$elapsed_time,
      os_time = fit_os$elapsed_time,
      rml_time = fit_rml$elapsed_time,
      lavaan_time = fit_lavaan$elapsed_time,

      svd_error = fit_svd$error_message,
      os_error = fit_os$error_message,
      rml_error = fit_rml$error_message,
      lavaan_error = fit_lavaan$error_message,

      stringsAsFactors = FALSE
    )
  }
}


# Combine all replications
results_mc <- do.call(rbind, results)

rownames(results_mc) <- NULL





# ============================================================
# 6. Summary by sample size
# ============================================================






summary_by_q <- do.call(
  rbind,
  lapply(
    split(results_mc, results_mc$q),
    summarize_one_q
  )
)

rownames(summary_by_q) <- NULL

print(summary_by_q)





# ============================================================
# 8. Export
# ============================================================

write.csv(
  results_mc,
  file = "inst/os_experiment/scaling_q_replications.csv",
  row.names = FALSE
)

write.csv(
  summary_by_q,
  file = "inst/os_experiment/scaling_q_summary.csv",
  row.names = FALSE
)




saveRDS(
  list(
    results_mc = results_mc,
    summary_by_q = summary_by_q,

    estimates_svd = estimates_svd,
    estimates_os = estimates_os,
    estimates_rml = estimates_rml,

    theta_true_list = theta_true_list,
    q_grid = q_grid,
    N = N,
    n_rep = n_rep
  ),
  file = "inst/os_experiment/monte_carlo_scaling_q_full.rds"
)