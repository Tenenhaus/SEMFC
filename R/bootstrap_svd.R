

#' Bootstrap Inference and Goodness-of-Fit Test for SVD-Estimated SEM Models
#'
#' Performs bootstrap-based statistical inference for models estimated using the SVD method.
#' Computes bootstrap distributions for all model parameters and conducts a Bollen-Stine
#' bootstrap test for model fit assessment using the Yuan & Hayashi (2003) data transformation approach.
#'
#' @param fit A svdSEM fit object containing the fitted model information, including
#'   parameter estimates, model-implied covariance matrix, and data blocks.
#' @param B Integer specifying the number of bootstrap replications (default: 100).
#' @param verbose Logical indicating whether to display a progress bar during
#'   bootstrap resampling (default: `TRUE`).
#'
#' @return A list containing bootstrap replications for all parameters:
#'   - `boot_lambda`: Matrix of bootstrap replications for factor loadings.
#'   - `boot_std_loadings`: Matrix of bootstrap replications for standardized loadings.
#'   - `boot_beta`: Matrix of bootstrap replications for structural coefficients (endogenous).
#'   - `boot_gamma`: Matrix of bootstrap replications for regression coefficients (exogenous).
#'   - `boot_residual_variance`: Matrix of bootstrap replications for residual variances.
#'   - `boot_total_effects`: Matrix of bootstrap replications for total effects.
#'   - `boot_indirect_effects`: Matrix of bootstrap replications for indirect effects.
#'   - `boot_omega`: Matrix of bootstrap replications for formative block weights.
#'   - `boot_Tb_LS`: Vector of bootstrap test statistics for goodness-of-fit.
#'   - `improper`: Integer count of improper bootstrap solutions (negative eigenvalues).
#'
#' @details The function performs two parallel bootstrap procedures:
#'   1. Standard bootstrap on original data for parameter inference
#'   2. Bollen-Stine bootstrap on transformed data for goodness-of-fit testing
#'
#'   The data transformation follows Yuan & Hayashi (2003). Bootstrap samples with
#'   negative eigenvalues are flagged as improper and excluded from standard error calculations.
#'
#' @examples
#' \dontrun{
#' boot_results <- bootstrap_svd(fit, B = 500, verbose = TRUE)
#' # Access bootstrap distributions
#' head(boot_results$boot_lambda)
#' # Check improper solutions
#' boot_results$improper
#' }
#'
#' @importFrom stats sd pnorm setNames
#' @keywords internal

bootstrap_svd <- function(fit, B = 100, verbose = TRUE, seed = NULL){

  df <- data.frame(Reduce("cbind", fit$blocks))


  Z0 <- lapply(split(data.frame(t(df)),
                    as.factor(rep(seq_along(fit$blocks),
                                  sapply(fit$blocks, NCOL)))),
              t)

  Z0_bs <- scaleDataSet(df, fit$sigma_implied)
  Z0_bs <- lapply(split(data.frame(t(Z0_bs)),
                  as.factor(rep(seq_along(fit$blocks),
                                sapply(fit$blocks, NCOL)))),
            t)



  if(verbose == FALSE){
    pbapply::pboptions(type = "none")
  }else{
    pbapply::pboptions(type = "txt")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  L <- pbapply::pbsapply(1:B,
                        function(b){

                          ind <- sample(NROW(df), replace = TRUE)

                          Zb <- lapply(Z0, function(x) x[ind, ])
                          Zb_bs <- lapply(Z0_bs, function(x) x[ind, ])
                          names(Zb) <- names(Zb_bs) <- names(fit$blocks)

                          SD <- apply(Reduce("cbind", Zb), 2, sd)
                          S_bs <- cov2(Reduce("cbind", Zb_bs), bias = fit$bias)


                          fit_b <- svdSEM(Zb,
                                         C = fit$C,
                                         scale = fit$scale,
                                         mode = fit$mode,
                                         bias = fit$bias)



                          fit_bs_b <- svdSEM(Zb_bs,
                                         C = fit$C,
                                         scale = fit$scale,
                                         mode = fit$mode,
                                         bias = fit$bias)


                          if (!any(eigen(fit_b$p_tilde)$values<=0)){
                            effect_b <- compute_effect(fit_b$beta, fit_b$gamma)
                            total_effects_b <- effect_b$total_effect
                            indirect_effects_b <- effect_b$indirect_effect

                            res <- list(
                              boot_lambda  = as.vector(Reduce("c", fit_b$lambda)),
                              boot_std_loadings = as.vector(Reduce("c", fit_b$lambda)/SD),
                              boot_beta = fit_b$beta[fit_b$beta!=0],
                              boot_gamma = fit_b$gamma[fit_b$gamma!=0],
                              boot_residual_variance = as.vector(Reduce("c", fit_b$residual_variance)),
                              boot_total_effects = as.vector(total_effects_b),
                              boot_indirect_effects = as.vector(indirect_effects_b),
                              boot_omega = as.vector(Reduce("c", fit_b$omega))
                            )

                          }
                          else{
                            res <- list(
                              boot_lambda = NA,
                              boot_std_loadings = NA,
                              boot_beta = NA,
                              boot_gamma = NA,
                              boot_residual_variance = NA,
                              boot_total_effects = NA,
                              boot_indirect_effects = NA,
                              boot_omega = NA
                            )
                          }
                          if (!any(eigen(fit_bs_b$p_tilde)$values<=0)){
                            res$Tb_LS <- d_LS(S_bs, fit_bs_b$sigma_implied)
                          }else{
                            res$Tb_LS <- NA
                          }

                          return(res)
                        })

  boot_lambda <- Reduce("rbind", L[1, ][!is.na(L[1, ])])
  boot_std_loadings <- Reduce("rbind", L[2, ][!is.na(L[2, ])])
  boot_beta <- Reduce("rbind", L[3, ][!is.na(L[3, ])])
  boot_gamma <- Reduce("rbind", L[4, ][!is.na(L[4, ])])
  boot_residual_variance <- Reduce("rbind", L[5, ][!is.na(L[5, ])])
  boot_total_effects <- Reduce("rbind", L[6, ][!is.na(L[6, ])])
  boot_indirect_effects <- Reduce("rbind", L[7, ][!is.na(L[7, ])])
  boot_omega <- Reduce("rbind", L[8, ][!is.na(L[8, ])])

  boot_Tb_LS <- unlist(L[9, ])
  improper <- sum(is.na(L[1, ]))
  improper_gof <- sum(is.na(L[9, ]))



  boot <- list(
    boot_lambda = boot_lambda,
    boot_std_loadings = boot_std_loadings,
    boot_beta = boot_beta,
    boot_gamma = boot_gamma,
    boot_residual_variance = boot_residual_variance,
    boot_total_effects = boot_total_effects,
    boot_indirect_effects = boot_indirect_effects,
    boot_omega = boot_omega,
    boot_Tb_LS = boot_Tb_LS,
    improper = improper,
    improper_gof = improper_gof
  )


  return(boot)
}


#' Compute Bootstrap Standard Errors for Model Parameters
#'
#' Calculates standard errors for all model parameters based on bootstrap replications.
#' The standard errors are computed as the standard deviation of the bootstrap distribution
#' for each parameter.
#'
#' @param boot List containing bootstrap replications. Expected elements include:
#'   - `boot_lambda`: Matrix of bootstrap replications for factor loadings.
#'   - `boot_std_loadings`: Matrix of bootstrap replications for standardized loadings.
#'   - `boot_beta`: Matrix of bootstrap replications for structural coefficients (endogenous).
#'   - `boot_gamma`: Matrix of bootstrap replications for regression coefficients (exogenous).
#'   - `boot_residual_variance`: Matrix of bootstrap replications for residual variances.
#'   - `boot_total_effects`: Matrix of bootstrap replications for total effects.
#'   - `boot_indirect_effects`: Matrix of bootstrap replications for indirect effects.
#'   - `boot_omega`: Matrix of bootstrap replications for formative block weights.
#'
#' @return A list containing standard errors for all parameters:
#'   - `sd_lambda`: Vector of standard errors for factor loadings.
#'   - `sd_std_loadings`: Vector of standard errors for standardized loadings.
#'   - `sd_beta`: Vector of standard errors for structural coefficients.
#'   - `sd_gamma`: Vector of standard errors for regression coefficients.
#'   - `sd_residual_variance`: Vector of standard errors for residual variances.
#'   - `sd_total_effects`: Vector of standard errors for total effects.
#'   - `sd_indirect_effects`: Vector of standard errors for indirect effects.
#'   - `sd_omega`: Vector of standard errors for formative block weights.
#'
#' @keywords internal

get_se_boot <- function(boot){

  safe_sd <- function(x) if (is.null(x)) NA else apply(x, 2, sd)

  sd_lambda <- safe_sd(boot$boot_lambda)
  sd_std_loadings <- safe_sd(boot$boot_std_loadings)
  sd_beta <- safe_sd(boot$boot_beta)
  sd_gamma <- safe_sd(boot$boot_gamma)
  sd_residual_variance <- safe_sd(boot$boot_residual_variance)
  sd_total_effects <- safe_sd(boot$boot_total_effects)
  sd_indirect_effects <- safe_sd(boot$boot_indirect_effects)
  sd_omega <- safe_sd(boot$boot_omega)

  return(list(
    sd_lambda = sd_lambda,
    sd_std_loadings = sd_std_loadings,
    sd_beta = sd_beta,
    sd_gamma = sd_gamma,
    sd_residual_variance = sd_residual_variance,
    sd_total_effects = sd_total_effects,
    sd_indirect_effects = sd_indirect_effects,
    sd_omega = sd_omega
  ))
}


#' Format Statistical Inference Results for SVD Estimation
#'
#' Formats parameter estimates with bootstrap-based standard errors, z-statistics,
#' and p-values into structured tables for SVD-estimated models.
#'
#' @param fit A svdSEM fit object containing the fitted model information, including
#'   parameter estimates for loadings, structural coefficients, and other model parameters.
#' @param boot List containing bootstrap replications from `bootstrap_svd()`. Expected elements include:
#'   - `boot_lambda`: Matrix of bootstrap replications for factor loadings.
#'   - `boot_std_loadings`: Matrix of bootstrap replications for standardized loadings.
#'   - `boot_beta`: Matrix of bootstrap replications for structural coefficients.
#'   - `boot_gamma`: Matrix of bootstrap replications for regression coefficients.
#'   - `boot_residual_variance`: Matrix of bootstrap replications for residual variances.
#'   - `boot_total_effects`: Matrix of bootstrap replications for total effects.
#'   - `boot_indirect_effects`: Matrix of bootstrap replications for indirect effects.
#'   - `boot_omega`: Matrix of bootstrap replications for formative block weights.
#'
#' @return A structured table containing formatted parameter estimates with bootstrap
#'   standard errors, z-statistics, and p-values.
#'
#' @keywords internal

formatting_svd_infer <- function(fit, boot){

  se <- get_se_boot(boot)
  table <- formatting_estimate(fit, se)

  return(table)

}




#' Statistical Inference and Goodness-of-Fit for SVD-Estimated Models
#'
#' Performs bootstrap-based statistical inference and computes goodness-of-fit statistics
#' for models estimated using the SVD method. Combines bootstrap replications with
#' Bollen-Stine bootstrap test results.
#'
#' @param fit A svdSEM fit object containing the fitted model information, including
#'   parameter estimates and model-implied covariance matrix.
#' @param B Integer specifying the number of bootstrap replications to perform.
#' @param verbose Logical indicating whether to display a progress bar during
#'   bootstrap resampling (default: `TRUE`).
#'
#' @return A list containing two main components:
#'   - `result`: List with bootstrap and inference results:
#'     - `boot`: List of bootstrap replications for all parameters (see `bootstrap_svd()`).
#'     - `infer`: Formatted tables with parameter estimates, standard errors, z-statistics,
#'       and p-values (see `formatting_svd_infer()`).
#'   - `gof`: List with goodness-of-fit test results:
#'     - `T_LS`: Observed test statistic from the fitted model.
#'     - `Tb_LS`: Vector of bootstrap test statistics.
#'     - `pval`: Bootstrap p-value for the goodness-of-fit test.
#'     - `improper`: Number of improper bootstrap solutions (negative eigenvalues).
#'
#' @keywords internal
svdsem_infer <- function(fit, B, verbose = TRUE, seed = seed){

  boot <- bootstrap_svd(fit, B, verbose, seed = seed)
  infer <- formatting_svd_infer(fit, boot)
  gof <- list(T_LS = fit$T_LS,
              Tb_LS = boot$boot_Tb_LS,
              pval = mean(fit$T_LS <= boot$boot_Tb_LS, na.rm = TRUE),
              improper = boot$improper
  )



  return(list(
    result = list(
      boot = boot,
      infer = infer
    ),
    gof = gof
    )
  )
}
