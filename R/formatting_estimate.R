




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
#' z-scores, confidence intervals and p-values.
#'
#' @param lhs Character vector of left-hand side variable names.
#' @param op Character vector of operators (e.g., "~", "=~", "~~").
#' @param rhs Character vector of right-hand side variable names.
#' @param est Numeric vector of parameter estimates.
#' @param se Numeric vector of standard errors (optional, defaults to NA).
#' @param std.all Numeric vector of standardized estimates (optional, defaults to NA).
#' @param alpha Significance level for confidence intervals (default: 0.05).
#'
#' @return Data frame with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all.
#'   Returns an empty data frame with the same structure if length(est) = 0.
#'
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

#' Formatting Estimates for Model Fit
#'
#' This function processes a fitted model object to extract and format its estimates
#' (lambda, gamma, beta, residual variance, total effects, and indirect effects)
#' into data frames. Each data frame includes columns for estimates, standard errors,
#' z-scores, and p-values, which are initialized as `NA`.
#'
#' @param fit A fitted model object containing the following components:
#'   - `lambda`: A named vector of loadings.
#'   - `gamma`: A matrix of path coefficients.
#'   - `beta`: A matrix of regression coefficients.
#'   - `residual_variance`: A vector of residual variances.
#'   - `effect$total_effect`: A matrix of total effects (non-zero values are used).
#'   - `effect$indirect_effect`: A matrix of indirect effects (non-zero values are used).
#'   - `omega`: A list of matrices containing weigths of indicators for composites blocs.
#'
#'
#' @return A list of data frames:
#'   - `lambda`: Data frame of loadings with columns for estimates, standard errors, z-scores, and p-values.
#'   - `gamma`: Data frame of path coefficients with the same columns as `lambda`.
#'   - `beta`: Data frame of regression coefficients with the same columns as `lambda`.
#'   - `residual_variance`: Data frame of residual variances with the same columns as `lambda`.
#'   - `total_effects`: Data frame of total effects with the same columns as `lambda`.
#'   - `indirect_effects`: Data frame of indirect effects with the same columns as `lambda`.
#'   - `omega`: Data frame of omega values with the same columns as `lambda`.
#'
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