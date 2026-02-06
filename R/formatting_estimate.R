
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

formatting_estimate <- function(fit){

  lambda <- unlist(fit$lambda)
  std_lambda <- unlist(fit$std_lambda)
  gamma <- fit$gamma[fit$gamma!=0]
  beta <- fit$beta[fit$beta!=0]
  residual_variance <- unlist(unname(fit$residual_variance))
  total_effects <- as.vector(fit$effect$total_effect)
  indirect_effects <- as.vector(fit$effect$indirect_effect)
  omega <- unlist(lapply(names(fit$omega), function(lv) {
    setNames(as.vector(fit$omega[[lv]]), paste(lv, rownames(fit$omega[[lv]]), sep = "."))
  }))


  parts <- strsplit(names(lambda), "\\.")
  table_lambda <- data.frame(
    lhs = sapply(parts, `[`, 1),
    op = "=~",
    rhs = sapply(parts, `[`, 2),
    est = lambda,
    se = NA,
    z = NA,
    ci.lower = NA,
    ci.upper = NA,
    pvalue = NA,
    std.all = std_lambda)

  # rownames(table_lambda) <- gsub("\\.", "~", rownames(table_lambda))

  table_std_lambda <- data.frame(
    lhs = sapply(parts, `[`, 1),
    op = "=~",
    rhs = sapply(parts, `[`, 2),
    est = std_lambda,
    se = NA,
    z = NA,
    ci.lower = NA,
    ci.upper = NA,
    pvalue = NA)
  # rownames(table_std_lambda) <- gsub("\\.", "~", rownames(table_std_lambda))

  table_gamma <- data.frame(est = gamma,
                            se = NA,
                            z = NA,
                            ci.lower = NA,
                            ci.upper = NA,
                            pvalue = NA,
                            std.all = gamma)
  rownames(table_gamma) <- sapply(seq_len(NROW(table_gamma)),
                                  function(b)
                          paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]],
                                colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                sep = ".")
  )
  parts <- strsplit(rownames(table_gamma), "\\.")
  table_gamma <- cbind(
    data.frame(
      lhs = sapply(parts, `[`, 1),
      op = "~",
      rhs = sapply(parts, `[`, 2)
    ),
    table_gamma
  )

  table_beta <- data.frame()

  if (length(beta) != 0){

    table_beta <- data.frame(est = beta,
                             se = NA,
                             z = NA,
                             ci.lower = NA,
                             ci.upper = NA,
                             pvalue = NA,
                             std.all = beta
    )
    rownames(table_beta) <- sapply(seq_len(NROW(table_beta)),
                                   function(b)
           paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
                 collapse = ".")
         )
    parts <- strsplit(rownames(table_beta), "\\.")
    table_beta <- cbind(
      data.frame(
        lhs = sapply(parts, `[`, 1),
        op = "~",
        rhs = sapply(parts, `[`, 2)
      ),    table_beta
    )
  }

  if (length(residual_variance) != 0){
    parts <- strsplit(names(residual_variance), "\\.")
    table_residual_variance <- data.frame(
      lhs = sapply(parts, `[`, 2),
      op = "~~",
      rhs = sapply(parts, `[`, 2),
      est = residual_variance,
      se = NA,
      z = NA,
      ci.lower = NA,
      ci.upper = NA,
      pvalue = NA)
    std_loadings_in_residuals <- table_std_lambda[table_std_lambda$rhs %in% table_residual_variance$rhs, "est"]
    table_std_residual_variance <- table_residual_variance
    table_std_residual_variance$est <- 1 - std_loadings_in_residuals^2
    table_residual_variance$std.all <- table_std_residual_variance$est

  } else {
    table_residual_variance <- data.frame()
    table_std_residual_variance <- data.frame()
  }



  grid_total_effects <- expand.grid(
    LHS = rownames(fit$effect$total_effect), RHS = colnames(fit$effect$total_effect)
  )
  table_total_effects <- data.frame(
    lhs = grid_total_effects$LHS,
    op = "~",
    rhs = grid_total_effects$RHS,
    est = total_effects,
    se = NA,
    z = NA,
    ci.lower = NA,
    ci.upper = NA,
    pvalue = NA,
    std.all = NA)



  rownames(table_total_effects) <- paste(grid_total_effects$LHS, grid_total_effects$RHS, sep = " ~ ")


  grid_indirect_effects <- expand.grid(
    LHS = rownames(fit$effect$indirect_effect), RHS = colnames(fit$effect$indirect_effect)
  )
  table_indirect_effects <- data.frame(
    lhs = grid_indirect_effects$LHS,
    op = "~",
    rhs = grid_indirect_effects$RHS,
    est = indirect_effects,
    se = NA,
    z = NA,
    ci.lower = NA,
    ci.upper = NA,
    pvalue = NA,
    std.all = NA)
  rownames(table_indirect_effects) <- paste(grid_indirect_effects$LHS, grid_indirect_effects$RHS, sep = " ~ ")

  table_omega <- data.frame()

  if (length(omega) != 0){
    parts <- strsplit(names(omega), "\\.")
    table_omega <- data.frame(
      lhs = sapply(parts, `[`, 1),
      op = "<~",
      rhs = sapply(parts, `[`, 2),
      est = omega,
      se = NA,
      z = NA,
      ci.lower = NA,
      ci.upper = NA,
      pvalue = NA,
      std.all = omega)
  }

  return(list(lambda = table_lambda,
              std_lambda = table_std_lambda,
              gamma = table_gamma,
              beta = table_beta,
              residual_variance = table_residual_variance,
              std_residual_variance = table_std_residual_variance,
              total_effects = table_total_effects,
              indirect_effects = table_indirect_effects,
              omega = table_omega))

}