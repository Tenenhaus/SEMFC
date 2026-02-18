


#' Statistical Inference  and Global Fit Assessment for svdSEM Models
#'
#' This function performs bootstrap-based statistical inference for svdSEM models,
#' computing standard errors, z-statistics, and p-values for all model parameters
#' including factor loadings, structural coefficients, and residual variances.
#' And in the same way, performs a bootstrap test to assess the significance of the fit.
#' The method is based on the Yuan & Hayashi (2003) approach for data transformation.
#'
#' @param fit A svdSEM fit object containing the fitted model information
#' @param B Number of bootstrap replications (default: 100)
#' @param verbose Logical indicating whether to display a progress bar (default: TRUE)
#'
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{infer}{A list containing statistical inference results with the following components:
#'     \itemize{
#'       \item \code{out}: List of bootstrap replications for lambda, standardized loadings,
#'         beta, gamma, residual variance, total effects, indirect effects, and omega
#'       \item \code{lambda}: Data frame with loadings, standard errors, z-statistics, and p-values
#'       \item \code{std_lambda}: Data frame with standardized loadings, standard errors,
#'         z-statistics, and p-values
#'       \item \code{beta}: Data frame with structural coefficients between latent variables,
#'         standard errors, z-statistics, and p-values
#'       \item \code{gamma}: Data frame with regression coefficients from observed to latent
#'         variables, standard errors, z-statistics, and p-values
#'       \item \code{residual_variance}: Data frame with residual variances, standard errors,
#'         z-statistics, and p-values
#'       \item \code{total_effects}: Data frame with total effects, standard errors,
#'         z-statistics, and p-values
#'       \item \code{indirect_effects}: Data frame with indirect effects, standard errors,
#'         z-statistics, and p-values
#'       \item \code{omega}: Data frame with error variances for latent variables,
#'         standard errors, z-statistics, and p-values
#'       \item \code{improper}: Number of improper bootstrap solutions (negative eigenvalues)
#'     }
#'   }
#'   \item{gof}{A list containing goodness-of-fit test results with the following components:
#'     \itemize{
#'       \item \code{T_LS}: The observed test statistic
#'       \item \code{Tb_LS}: Vector of bootstrap test statistics
#'       \item \code{pval}: The p-value of the fit test
#'       \item \code{improper}: Number of improper solutions (negative eigenvalues)
#'     }
#'   }
#' }
#'
#'
#' @details The function uses bootstrap resampling to estimate the sampling distribution
#' of all model parameters. Z-statistics are computed as parameter estimates divided by
#' bootstrap standard errors, and p-values are calculated using a two-tailed normal test.
#' Bootstrap samples with negative eigenvalues are considered improper and excluded from
#' the analysis. In the same time, the function  transforms the data according to Yuan & Hayashi (2003),
#' then performs bootstrap to estimate the empirical distribution of the test statistic.
#' Solutions with negative eigenvalues are considered improper and excluded.
#'
#' @examples
#' \dontrun{
#' bootstrap_results <- bootstrap_svd(fit, B = 100, bias = FALSE)
#' print(bootstrap_results$infer$out$boot_lambda)
#' print(bootstrap_results$gof$pval)
#' }
#'
#' @importFrom stats sd pnorm setNames
#' @export










bootstrap_svd <- function(fit, B = 100, verbose = TRUE){

  df <- data.frame(Reduce("cbind", fit$blocks))


  Z0 <- lapply(split(data.frame(t(df)),
                    as.factor(rep(seq_along(fit$blocks),
                                  sapply(fit$blocks, NCOL)))),
              t)

  Z0_bs <- scaleDataSet(df, fit$SIGMA_IMPLIED)
  Z0_bs <- lapply(split(data.frame(t(Z0_bs)),
                  as.factor(rep(seq_along(fit$blocks),
                                sapply(fit$blocks, NCOL)))),
            t)



  if(verbose == FALSE){
    pbapply::pboptions(type = "none")
  }else{
    pbapply::pboptions(type = "txt")
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


                          if (!any(eigen(fit_b$Ptilde)$values<=0)){
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
                          if (!any(eigen(fit_bs_b$Ptilde)$values<=0)){
                            res$Tb_LS <- NA
                          }else{
                            res$Tb_LS <- d_LS(S_bs, fit_bs_b$SIGMA_IMPLIED)
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
    improper = improper
  )


  return(boot)
}




get_se_boot <- function(boot){
  sd_lambda <- apply(boot$boot_lambda, 2, sd)
  sd_std_loadings <- apply(boot$boot_std_loadings, 2, sd)
  sd_beta <- apply(boot$boot_beta, 2, sd)
  sd_gamma <- apply(boot$boot_gamma, 2, sd)
  sd_residual_variance <- apply(boot$boot_residual_variance, 2, sd)
  sd_total_effects <- apply(boot$boot_total_effects, 2, sd)
  sd_indirect_effects <- apply(boot$boot_indirect_effects, 2, sd)
  sd_omega <- apply(boot$boot_omega, 2, sd)

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




formatting_svd_infer <- function(fit, boot){

  se <- get_se_boot(boot)
  table <- formatting_estimate(fit, se)

  return(table)

}





svdsem_infer <- function(fit, B, verbose = TRUE){

  boot <- bootstrap_svd(fit, B, verbose)
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
