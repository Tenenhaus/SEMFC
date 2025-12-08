


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
#'   \item{infer}{A list containing statistical inference results:
#'     \itemize{
#'       \item \code{out}: List of bootstrap replications for lambda, standardized loadings, beta, gamma, and residual variance
#'       \item \code{lambda}: Data frame with factor loadings, standard errors, z-statistics, and p-values
#'       \item \code{std_lambda}: Data frame with standardized loadings, standard errors, z-statistics, and p-values
#'       \item \code{beta}: Data frame with structural coefficients between latent variables, standard errors, z-statistics, and p-values
#'       \item \code{gamma}: Data frame with regression coefficients from observed to latent variables, standard errors, z-statistics, and p-values
#'       \item \code{residual_variance}: Data frame with residual variances, standard errors, z-statistics, and p-values
#'       \item \code{improper}: Number of improper bootstrap solutions (negative eigenvalues) encountered
#'     }
#'   }
#'   \item{gof}{A list containing goodness-of-fit test results:
#'     \itemize{
#'       \item \code{T_LS}: The observed test statistic
#'       \item \code{Tb_LS}: Vector of bootstrap test statistics
#'       \item \code{pval}: The p-value of the fit test
#'       \item \code{improper}: Number of improper solutions (negative eigenvalues) encountered
#'     }
#'   }
#'
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
#'
#' @export










bootstrap_svd <- function(fit, B = 100, verbose = TRUE){

  df <- data.frame(Reduce("cbind", fit$blocks))
  sd_init <- apply(df, 2, sd)

  beta <- fit$beta[fit$beta!=0]
  gamma <- fit$gamma[fit$gamma!=0]
  lambda <- unlist(unname(fit$lambda))
  # names(lambda) <- unlist(lapply(fit$lambda, names))
  std_loadings <- lambda/sd_init
  residual_variance <- unlist(unname(fit$residual_variance))

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
                            res <- list(
                              boot_lambda  = as.vector(Reduce("c", fit_b$lambda)),
                              boot_std_loadings = as.vector(Reduce("c", fit_b$lambda)/SD),
                              boot_beta = fit_b$beta[fit_b$beta!=0],
                              boot_gamma = fit_b$gamma[fit_b$gamma!=0],
                              boot_residual_variance = as.vector(Reduce("c", fit_b$residual_variance)))

                          }
                          else{
                            res <- list(
                              boot_lambda = NA,
                              boot_std_loadings = NA,
                              boot_beta = NA,
                              boot_gamma = NA,
                              boot_residual_variance = NA)
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
  boot_Tb_LS <- unlist(L[6, ])


  std_residual_variance <- apply(boot_residual_variance, 2, sd)
  t_ratio <- residual_variance/std_residual_variance
  pval_residual_variance <- sapply(seq_along(t_ratio),
                                   function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )
  residual_variance <- data.frame(lambda = residual_variance,
                      std = std_residual_variance,
                      z = t_ratio,
                      pval = pval_residual_variance)



  std_lambda <- apply(boot_lambda, 2, sd)

  t_ratio <- lambda/std_lambda
  pval_lambda <- sapply(seq_along(t_ratio),
                        function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )


  lambda <- data.frame(lambda = lambda,
                      std = std_lambda,
                      z = t_ratio,
                      pval = pval_lambda)

  std_std_loadings <- apply(boot_std_loadings, 2, sd)

  t_ratio <- std_loadings/std_std_loadings
  pval_std_loadings <- sapply(seq_along(t_ratio),
                              function(x)
                                2*pnorm(t_ratio[x],
                                        lower.tail = FALSE)
  )


  std_lambda <- data.frame(std_loadings = std_loadings,
                          std = std_std_loadings,
                          z = t_ratio,
                          pval = pval_std_loadings)

  std_beta <- apply(boot_beta, 2, sd)
  t_ratio <- beta[beta!=0]/apply(boot_beta, 2, sd)
  pval_beta <- sapply(seq_along(t_ratio),
                      function(x)
                           2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )

  # beta <- data.frame()
  if (length(beta) != 0){

    beta <- data.frame(beta = beta,
                    std = std_beta,
                    z = t_ratio,
                    pval = pval_beta)

    rownames(beta) <- sapply(seq_len(NROW(beta)),
                             function(b)
             paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
                   collapse = "~")
           )

  }
  beta <- data.frame(beta)


  std_gamma <- apply(boot_gamma, 2, sd)

  t_ratio <- gamma[gamma!=0]/std_gamma
  pval_gamma <- sapply(seq_along(t_ratio),
                       function(x)
                            2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )

  gamma <- data.frame(gamma = gamma,
                     std = std_gamma,
                     z = t_ratio,
                     pval = pval_gamma)

  rownames(gamma) <- sapply(seq_len(NROW(gamma)),
                            function(b)
                            paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]],
                                  colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                  sep = "~")
  )


  return(list(infer = list(out = list(boot_lambda,
                                      boot_std_loadings,
                                      boot_beta,
                                      boot_gamma,
                                      boot_residual_variance),
                           lambda = lambda,
                           std_lambda = std_lambda,
                           beta = beta,
                           gamma = gamma,
                           residual_variance = residual_variance,
                           improper = sum(is.na(L[1, ]))),
              gof = list(T_LS = fit$T_LS,
                         Tb_LS = boot_Tb_LS,
                         pval = mean(fit$T_LS <= boot_Tb_LS, na.rm = TRUE),
                         improper = sum(is.na(boot_Tb_LS))
              )
  ))
}
