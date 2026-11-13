#' Bootstrap Inference and Goodness-of-Fit Test for SVD-Estimated SEM Models
#'
#' Performs bootstrap-based statistical inference for models estimated using the SVD method.
#' Computes bootstrap distributions for all model parameters and conducts a Bollen-Stine
#' bootstrap test for model fit assessment using the Yuan & Hayashi (2003) data transformation approach.
#'
#' @param fit A svdSEM fit object containing the fitted model information, including
#'   parameter estimates, model-implied covariance matrix, and data blocks.
#' @param B Integer specifying the number of bootstrap replications (default: 500).
#' @param verbose Logical indicating whether to display a progress bar during
#'   bootstrap resampling (default: `TRUE`).
#' @param seed Optional integer for reproducibility of bootstrap results. If `NULL`, no seed is set (default: `NULL`).
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
#'   negative eigenvalues are flagged as improper and excluded from standard
#'   error calculations.
#'
#' \strong{References}
#'
#' \enumerate{
#' \item Yuan, K. H., & Hayashi, K. (2003). Bootstrap approach to inference and
#'   power analysis based on three test statistics for covariance structure
#'   models. British Journal of Mathematical and Statistical Psychology,
#'   56(1), 93-110.
#' }
#'
#'
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

bootstrap_svd <- function(fit, B = 500, verbose = TRUE, seed = NULL) {
  df <- data.frame(Reduce("cbind", fit$blocks)) # le dataframe des data


  Z0 <- lapply(
    split(
      data.frame(t(df)),
      as.factor(rep(
        seq_along(fit$blocks),
        sapply(fit$blocks, NCOL)
      ))
    ),
    t
  ) # maniere tres compliquee de juste separer df en plusieurs blocs... pourquoi ne pas juste utiliser fit$blocks.


  Z0_bs <- scaleDataSet(df, fit$SIGMA_IMPLIED, bias = fit$bias) # rescaling bollen stine bootstrap
  Z0_bs <- lapply(
    split(
      data.frame(t(Z0_bs)),
      as.factor(rep(
        seq_along(fit$blocks),
        sapply(fit$blocks, NCOL)
      ))
    ),
    t
  )

  # Maintenant on a le second jeu de data pour bollen stine

  if (verbose == FALSE) {
    pbapply::pboptions(type = "none")
  } else {
    pbapply::pboptions(type = "txt")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }


  L <- pbapply::pblapply(
    1:B,
    function(b) {
      ind <- sample(NROW(df), replace = TRUE)

      Zb <- lapply(Z0, function(x) x[ind, ])
      Zb_bs <- lapply(Z0_bs, function(x) x[ind, ])
      names(Zb) <- names(Zb_bs) <- names(fit$blocks)

      S_bs <- cov2(
        scale2(
          Reduce("cbind", Zb_bs),
          bias = fit$bias, scale = fit$scale
        ),
        bias = fit$bias
      )

      fit_b <- svdSEM_multi(Zb,
        fit$li_C,
        scale = fit$scale,
        mode = fit$mode,
        bias = fit$bias
      )

      fit_bs_b <- svdSEM_multi(Zb_bs,
        fit$li_C,
        scale = fit$scale,
        mode = fit$mode,
        bias = fit$bias
      )
      res <- list()
      proper_boot <- !any(eigen(fit_b$P_IMPLIED, symmetric = TRUE, only.values = TRUE)$values <= 0)
      if (proper_boot) {
        res$boot <- list(
          boot_Lambda <- fit_b$li_Lambda,
          boot_gamma <- fit_b$li_gamma,
          boot_beta <- fit_b$li_beta,
          boot_phi_endo <- fit_b$li_phi_endo,
          boot_phi_exo <- fit_b$li_phi_exo,
          boot_psi <- fit_b$li_PSI
        )
      } else {
        res$boot <- NULL
      }
      proper_bs <- !any(eigen(fit_bs_b$P_IMPLIED, symmetric = TRUE, only.values = TRUE)$values <= 0)
      if (proper_bs) {
        res$BS <- d_LS(S_bs, fit_bs_b$SIGMA_IMPLIED)
      } else {
        res$BS <- NULL
      }
      return(res)
    }
  )


  boot_Lambda <- lapply(L, function(x) if (!(is.null(x$boot))) x$boot[[1]] else NULL)
  boot_gamma <- lapply(L, function(x) if (!is.null(x$boot)) x$boot[[2]] else NULL)
  boot_beta <- lapply(L, function(x) if (!is.null(x$boot)) x$boot[[3]] else NULL)
  boot_phi_endo <- lapply(L, function(x) if (!is.null(x$boot)) x$boot[[4]] else NULL)
  boot_phi_exo <- lapply(L, function(x) if (!is.null(x$boot)) x$boot[[5]] else NULL)
  boot_psi <- lapply(L, function(x) if (!is.null(x$boot)) x$boot[[6]] else NULL)
  boot_BS <- unlist(lapply(L, function(x) x$BS))

  improper <- sum(sapply(L, function(x) is.null(x$boot)))
  improper_gof <- sum(sapply(L, function(x) is.null(x$BS)))


  # print(L)
  # print(improper)
  # print(improper_gof)
  J <- length(fit$li_Lambda)
  R <- ncol(fit$li_Lambda[[1]])
  n_endogene <- nrow(fit$li_gamma[[1]])
  n_exogene <- ncol(fit$li_gamma[[1]])

  boot <- list(
    boot_Lambda = boot_Lambda,
    boot_gamma = boot_gamma,
    boot_beta = boot_beta,
    boot_phi_endo = boot_phi_endo,
    boot_phi_exo = boot_phi_exo,
    boot_psi = boot_psi,
    boot_BS = boot_BS,
    improper = improper,
    improper_gof = improper_gof,
    R = R,
    J = J,
    n_endogene = n_endogene,
    n_exogene = n_exogene
  )


  return(boot)
}


#' Compute Bootstrap Standard Errors for Model Parameters


get_se_boot <- function(boot) {
  safe_sd <- function(x) if (is.null(x)) NA else apply(x, 2, sd)
  J <- boot$J
  R <- boot$R
  n_endogene <- boot$n_endogene
  n_exogene <- boot$n_exogene


  li_sd_Lambda <- lapply(1:J, function(j) {
    li_Lambda_j <- lapply(boot$boot_Lambda, function(x) x[[j]])
    lambda_j_sd <- matrix(NA, nrow = nrow(li_Lambda_j[[1]]), ncol = ncol(li_Lambda_j[[1]]))
    for (i in 1:vec_indicators_per_bloc[j]) {
      for (k in 1:R) {
        vec_values <- sapply(li_Lambda_j, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        lambda_j_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(lambda_j_sd)
  })

  li_sd_gamma <- lapply(1:R, function(r) {
    li_gamma_r <- lapply(boot$boot_gamma, function(x) x[[r]])
    gamma_r_sd <- matrix(NA, nrow = n_endogene, ncol = n_exogene)
    for (i in 1:n_endogene) {
      for (k in 1:n_exogene) {
        vec_values <- sapply(li_gamma_r, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        gamma_r_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(gamma_r_sd)
  })

  li_sd_beta <- lapply(1:R, function(r) {
    li_beta_r <- lapply(boot$boot_beta, function(x) x[[r]])
    beta_r_sd <- matrix(NA, nrow = n_endogene, ncol = n_endogene)
    for (i in 1:n_endogene) {
      for (k in 1:n_endogene) {
        vec_values <- sapply(li_beta_r, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        beta_r_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(beta_r_sd)
  })

  li_sd_phi_endo <- lapply(1:R, function(r) {
    li_phi_endo_r <- lapply(boot$boot_phi_endo, function(x) x[[r]])
    phi_endo_r_sd <- matrix(NA, nrow = n_endogene, ncol = n_endogene)
    for (i in 1:n_endogene) {
      for (k in 1:n_endogene) {
        vec_values <- sapply(li_phi_endo_r, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        phi_endo_r_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(phi_endo_r_sd)
  })

  li_sd_phi_exo <- lapply(1:R, function(r) {
    li_phi_exo_r <- lapply(boot$boot_phi_exo, function(x) x[[r]])
    phi_exo_r_sd <- matrix(NA, nrow = n_exogene, ncol = n_exogene)
    for (i in 1:n_exogene) {
      for (k in 1:n_exogene) {
        vec_values <- sapply(li_phi_exo_r, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        phi_exo_r_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(phi_exo_r_sd)
  })

  li_sd_psi <- lapply(1:R, function(r) {
    li_psi_r <- lapply(boot$boot_psi, function(x) x[[r]])
    psi_r_sd <- matrix(NA, nrow = n_endogene, ncol = n_endogene)
    for (i in 1:n_endogene) {
      for (k in 1:n_endogene) {
        vec_values <- sapply(li_psi_r, function(x) {
          if (is.null(x)) {
            return(NA)
          } else {
            return(x[i, k])
          }
        })
        psi_r_sd[i, k] <- sd(vec_values, na.rm = TRUE)
      }
    }
    return(psi_r_sd)
  })


  n_improper_sd <- boot$improper

  return(list(
    li_sd_Lambda = li_sd_Lambda,
    li_sd_gamma = li_sd_gamma,
    li_sd_beta = li_sd_beta,
    li_sd_phi_endo = li_sd_phi_endo,
    li_sd_phi_exo = li_sd_phi_exo,
    li_sd_psi = li_sd_psi,
    n_improper_sd = n_improper_sd
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

formatting_svd_infer <- function(fit, boot) {
  se <- get_se_boot(boot)
  table <- formatting_estimate(fit, se)

  return(table)
}



svdsem_infer <- function(fit, B, verbose = TRUE, seed = NULL) {
  boot <- bootstrap_svd(fit, B, verbose, seed = seed)
  infer <- formatting_svd_infer(fit, boot)
  gof <- list(
    BS_stat = fit$BS_stat,
    boot_BS = boot$boot_BS,
    pval = mean(fit$BS_stat <= boot$boot_BS, na.rm = TRUE),
    improper_gof = boot$improper_gof
  )


  return(list(
    result = list(
      boot = boot,
      infer = infer
    ),
    gof = gof
  ))
}
