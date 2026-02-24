




#' Extract Non-Zero Matrix Elements with Row and Column Names
#'
#' Creates a grid of row-column pairs for non-zero elements in a matrix.
#'
#' @param matrix A matrix with named rows and columns.
#'
#' @return Data frame with columns: lhs (row names), rhs (column names).
#'
#' @keywords internal
extract_matrix_pairs <- function(matrix) {
  indices <- which(matrix != 0, arr.ind = TRUE)
  data.frame(
    lhs = rownames(matrix)[indices[, 1]],
    rhs = colnames(matrix)[indices[, 2]]
  )
}



#' Create Grid Based on Type
#'
#' Creates a grid of lhs-rhs pairs based on the estimation type.
#'
#' @param type Character string indicating the type: "lambda", "gamma", "beta", "residual", "total", "indirect", or "omega".
#' @param fit_component The relevant component from the fit object.
#'
#' @return Data frame with columns: lhs, rhs, op.
#'
#' @keywords internal
create_grid <- function(type, fit_component) {
  switch(type,
    "lambda" = {
      grid <- do.call(rbind, lapply(names(fit_component), function(lv_name) {
        expand.grid(lhs = lv_name, rhs = names(fit_component[[lv_name]]))
      }))
      grid$op <- "=~"
      grid
    },
    "regression" = {
      grid <- extract_matrix_pairs(fit_component)
      grid$op <- "~"
      grid
    },
    "residual" = {
      grid_temp <- do.call(rbind, lapply(names(fit_component), function(lv_name) {
        expand.grid(lhs = lv_name, rhs = gsub("^\\.", "", names(fit_component[[lv_name]])))
      }))
      grid <- data.frame(lhs = grid_temp$rhs, rhs = grid_temp$rhs, op = "~~")
      grid
    },
    "effect" = {
      grid <- expand.grid(lhs = rownames(fit_component), rhs = colnames(fit_component))
      grid$op <- "~"
      grid
    },
    "omega" = {
      grid <- do.call(rbind, lapply(names(fit_component), function(lv_name) {
        expand.grid(lhs = lv_name, rhs = rownames(fit_component[[lv_name]]))
      }))
      grid$op <- "<~"
      grid
    }
  )
}




#' Format Parameter Estimates Table
#'
#' Creates a standardized data frame with parameter estimates, standard errors,
#' z-scores, confidence intervals and p-values based on the estimation type.
#'
#' @param type Character string indicating the type: "lambda", "regression", "residual", "effect", or "omega".
#' @param fit_component The relevant component from the fit object (e.g., fit$lambda, fit$gamma, fit$beta).
#' @param fit_std_component Standardized estimates component (optional, defaults to NULL).
#' @param se Numeric vector of standard errors (optional, defaults to NA).
#' @param alpha Significance level for confidence intervals (default: 0.05).
#'
#' @return Data frame with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all.
#'   Returns an empty data frame with the same structure if length(est) = 0.
#' @importFrom stats qnorm
#' @keywords internal
format_estimates_table <- function(type, fit_component, fit_std_component = NULL, se = NA, alpha = 0.05){

  est <- switch(type,
    "lambda" = unlist(fit_component),
    "regression" = fit_component[fit_component != 0],
    "residual" = unlist(unname(fit_component)),
    "effect" = as.vector(fit_component),
    "omega" = unlist(lapply(names(fit_component), function(lv) {
      as.vector(fit_component[[lv]])
    }))
  )


  # Return empty data frame if no estimates
  if (length(est) == 0) {
    return(data.frame())
  }

  grid <- create_grid(type, fit_component)


  std.all <- if (!is.null(fit_std_component)) {
    switch(type,
      "lambda" = unlist(fit_std_component),
      "regression" = fit_std_component,
      "omega" = unlist(fit_std_component),
      NA
    )
  } else {
    NA
  }

  z_critical <- qnorm(1 - alpha / 2)
  if (all(is.na(se))) {
    se <- rep(NA, length(est))
    z_scores <- rep(NA, length(est))
    ci_lower <- rep(NA, length(est))
    ci_upper <- rep(NA, length(est))
    pvalue <- rep(NA, length(est))
  } else {
    z_scores <- est / se
    ci_lower <- est - z_critical * se
    ci_upper <- est + z_critical * se
    pvalue <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
  }

  data.frame(
    lhs = grid$lhs,
    op = grid$op,
    rhs = grid$rhs,
    est = est,
    se = se,
    z = z_scores,
    ci.lower = ci_lower,
    ci.upper = ci_upper,
    pvalue = pvalue,
    std.all = if (all(is.na(std.all))) NA else std.all
  )
}

#' Format Model Estimates
#'
#' Processes a fitted model object to extract and format parameter estimates
#' into standardized data frames. Creates tables for loadings, path coefficients,
#' regression coefficients, residual variances, total effects, indirect effects,
#' and composite weights.
#'
#' @param fit A fitted model object containing the following components:
#'   - `lambda`: A named list of loadings for each latent variable.
#'   - `std_lambda`: Standardized loadings (optional).
#'   - `gamma`: A matrix of path coefficients between latent variables.
#'   - `beta`: A matrix of regression coefficients.
#'   - `residual_variance`: A named list of residual variances.
#'   - `effect$total_effect`: A matrix of total effects.
#'   - `effect$indirect_effect`: A matrix of indirect effects.
#'   - `omega`: A list of matrices containing weights of indicators for composite blocks.
#' @param se_list A list containing standard errors for each component (optional, defaults to empty list).
#'   Expected elements: `sd_lambda`, `sd_gamma`, `sd_beta`, `sd_residual_variance`,
#'   `sd_total_effects`, `sd_indirect_effects`, `sd_omega`.
#'
#' @return A list of data frames with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all:
#'   - `lambda`: Loadings (op = "=~").
#'   - `gamma`: Path coefficients between latent variables (op = "~").
#'   - `beta`: Regression coefficients (op = "~").
#'   - `residual_variance`: Residual variances (op = "~~").
#'   - `total_effects`: Total effects (op = "~").
#'   - `indirect_effects`: Indirect effects (op = "~").
#'   - `omega`: Composite weights (op = "<~").
#'
#' @keywords internal

formatting_estimate <- function(fit, se_list = list()){
  # Lambda
  table_lambda <- format_estimates_table(
    type = "lambda",
    fit_component = fit$lambda,
    fit_std_component = fit$std_lambda,
    se = se_list$sd_lambda
  )

  # Gamma
  table_gamma <- format_estimates_table(
    type = "regression",
    fit_component = fit$gamma,
    fit_std_component = fit$gamma[fit$gamma != 0],
    se = se_list$sd_gamma
  )

  # Beta
  table_beta <- format_estimates_table(
    type = "regression",
    fit_component = fit$beta,
    fit_std_component = fit$beta[fit$beta != 0],
    se = se_list$sd_beta
  )

  # Residual variance
  table_residual_variance <- format_estimates_table(
    type = "residual",
    fit_component = fit$residual_variance,
    se = se_list$sd_residual_variance
  )
  table_residual_variance$std.all <- 1 - (table_lambda[table_lambda$rhs %in% table_residual_variance$rhs, "std.all"])^2

  # Total effects
  table_total_effects <- format_estimates_table(
    type = "effect",
    fit_component = fit$effect$total_effect,
    se = se_list$sd_total_effects
  )
  rownames(table_total_effects) <- paste(table_total_effects$lhs, table_total_effects$rhs, sep = " ~ ")

  # Indirect effects
  table_indirect_effects <- format_estimates_table(
    type = "effect",
    fit_component = fit$effect$indirect_effect,
    se = se_list$sd_indirect_effects
  )
  rownames(table_indirect_effects) <- paste(table_indirect_effects$lhs, table_indirect_effects$rhs, sep = " ~ ")

  # Omega
  table_omega <- format_estimates_table(
    type = "omega",
    fit_component = fit$omega,
    fit_std_component = fit$std_omega,
    se = se_list$sd_omega
  )

  return(list(lambda = table_lambda,
              gamma = table_gamma,
              beta = table_beta,
              residual_variance = table_residual_variance,
              total_effects = table_total_effects,
              indirect_effects = table_indirect_effects,
              omega = table_omega))

}