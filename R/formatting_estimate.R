
#' Formatting Estimates for Model Fit
#'
#' This function processes a fitted model object to extract and format its estimates
#' (lambda, gamma, beta, residual variance, total effects, and indirect effects)
#' into data frames. Each data frame includes columns for estimates, standard errors,
#' z-scores, and p-values, which are initialized as `NA`.
#'
#' @param fit A fitted model object containing the following components:
#'   - `lambda`: A named vector of loadings.
#'   - `gamma`: A matrix of path coefficients (non-zero values are used).
#'   - `beta`: A matrix of regression coefficients (non-zero values are used).
#'   - `residual_variance`: A vector of residual variances.
#'   - `effect$total_effect`: A matrix of total effects (non-zero values are used).
#'   - `effect$indirect_effect`: A matrix of indirect effects (non-zero values are used).
#'
#' @return A list of data frames:
#'   - `lambda`: Data frame of loadings with columns for estimates, standard errors, z-scores, and p-values.
#'   - `gamma`: Data frame of path coefficients with the same columns as `lambda`.
#'   - `beta`: Data frame of regression coefficients with the same columns as `lambda`.
#'   - `residual_variance`: Data frame of residual variances with the same columns as `lambda`.
#'   - `total_effects`: Data frame of total effects with the same columns as `lambda`.
#'   - `indirect_effects`: Data frame of indirect effects with the same columns as `lambda`.
#'
#' @examples
#' \dontrun{
#' # Assuming `fit` is a fitted model object:
#' formatted_estimates <- formatting_estimate(fit)
#' print(formatted_estimates$lambda)
#' }

formatting_estimate <- function(fit){

  lambda <- unlist(fit$lambda)
  gamma <- fit$gamma[fit$gamma!=0]
  beta <- fit$beta[fit$beta!=0]
  residual_variance <- unlist(unname(fit$residual_variance))
  total_effects <- fit$effect$total_effect[fit$effect$total_effect!=0]
  indirect_effects <- fit$effect$indirect_effect[fit$effect$total_effect!=0]


  table_lambda <- data.frame(Estimate = lambda,
                           std = NA,
                           z_score = NA,
                           pval = NA)

  rownames(table_lambda) <- gsub("\\.", "~", rownames(table_lambda))

  table_gamma <- data.frame(Estimate = gamma,
                            std = NA,
                            z_score = NA,
                            pval = NA)
  rownames(table_gamma) <- sapply(seq_len(NROW(table_gamma)),
                                  function(b)
                          paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]],
                                colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                sep = "~")
  )

  table_beta <- data.frame()

  if (length(beta) != 0){

    table_beta <- data.frame(Estimate = beta,
                             std = NA,
                             z_score = NA,
                             pval = NA
    )
    rownames(table_beta) <- sapply(seq_len(NROW(table_beta)),
                                   function(b)
           paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
                 collapse = "~")
         )
  }

  table_residual_variance <- data.frame(Estimate = residual_variance,
                           std = NA,
                           z_score = NA,
                           pval = NA)

  table_total_effects <- data.frame(Estimate = total_effects,
                          std = NA,
                          z_score = NA,
                          pval = NA)
  rownames(table_total_effects) <- sapply(seq_len(NROW(table_total_effects)),
                                  function(b)
                          paste(rownames(fit$effect$total_effect)[which(fit$effect$total_effect!=0, arr.ind = TRUE)[b, 1]],
                                colnames(fit$effect$total_effect)[which(fit$effect$total_effect!=0, arr.ind = TRUE)[b, 2]],
                                sep = "~")
  )

  table_indirect_effects <- data.frame(Estimate = indirect_effects,
                        std = NA,
                        z_score = NA,
                        pval = NA)
  rownames(table_indirect_effects) <- sapply(seq_len(NROW(table_indirect_effects)),
                                  function(b)
                          paste(rownames(fit$effect$indirect_effect)[which(fit$effect$indirect_effect!=0, arr.ind = TRUE)[b, 1]],
                                colnames(fit$effect$indirect_effect)[which(fit$effect$indirect_effect!=0, arr.ind = TRUE)[b, 2]],
                                sep = "~")
  )



  return(list(lambda = table_lambda,
              gamma = table_gamma,
              beta = table_beta,
              residual_variance = table_residual_variance,
              total_effects = table_total_effects,
              indirect_effects = table_indirect_effects))





}