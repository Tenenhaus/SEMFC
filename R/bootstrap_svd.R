


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
  sd_init <- apply(df, 2, sd)

  beta <- fit$beta[fit$beta!=0]
  gamma <- fit$gamma[fit$gamma!=0]
  lambda <- unlist(fit$lambda)
  std_loadings <- lambda/sd_init
  residual_variance <- unlist(unname(fit$residual_variance))
  total_effects <- as.vector(fit$effect$total_effect)
  indirect_effects <- as.vector(fit$effect$indirect_effect)
  omega <- unlist(lapply(names(fit$omega), function(lv) {
    setNames(as.vector(fit$omega[[lv]]), paste(lv, rownames(fit$omega[[lv]]), sep = "."))
  }))

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


  std_residual_variance <- apply(boot_residual_variance, 2, sd)
  t_ratio <- residual_variance/std_residual_variance
  pval_residual_variance <- sapply(seq_along(t_ratio),
                                   function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )

  parts <- strsplit(names(residual_variance), "\\.")
  residual_variance <- data.frame(
    lhs = sapply(parts, `[`, 2),
    op = "~~",
    rhs = sapply(parts, `[`, 2),
    est = residual_variance,
    se = std_residual_variance,
    z = t_ratio,
    pvalue = pval_residual_variance)



  std_lambda <- apply(boot_lambda, 2, sd)

  t_ratio <- lambda/std_lambda
  pval_lambda <- sapply(seq_along(t_ratio),
                        function(x)
                               2*pnorm(abs(t_ratio[x]),
                                       lower.tail = FALSE)
  )

  parts <- strsplit(names(lambda), "\\.")
  lambda <- data.frame(
    lhs = sapply(parts, `[`, 1),
    op = "=~",
    rhs = sapply(parts, `[`, 2),
    est = lambda,
    se = std_lambda,
    z = t_ratio,
    pvalue = pval_lambda)
  # rownames(lambda) <- gsub("\\.", "~", rownames(lambda))



  if (length(omega) != 0){
    std_omega <- apply(boot_omega, 2, sd)
    t_ratio <- omega/std_omega
    pval_omega <- sapply(seq_along(t_ratio),
                         function(x)
                           2*pnorm(abs(t_ratio[x]),
                                   lower.tail = FALSE)
    )
    parts <- strsplit(names(omega), "\\.")
    omega <- data.frame(
      lhs = sapply(parts, `[`, 1),
      op = "<~",
      rhs = sapply(parts, `[`, 2),
      est = omega,
      se = std_omega,
      z = t_ratio,
      pvalue = pval_omega)

  }
  omega <- data.frame(omega)





  std_std_loadings <- apply(boot_std_loadings, 2, sd)

  t_ratio <- std_loadings/std_std_loadings
  pval_std_loadings <- sapply(seq_along(t_ratio),
                              function(x)
                                2*pnorm(t_ratio[x],
                                        lower.tail = FALSE)
  )

  parts <- strsplit(names(std_loadings), "\\.")
  std_lambda <- data.frame(
    lhs = sapply(parts, `[`, 1),
    op = "=~",
    rhs = sapply(parts, `[`, 2),
    est = std_loadings,
    se = std_std_loadings,
    z = t_ratio,
    pvalue = pval_std_loadings)

  rownames(std_lambda) <- gsub("\\.", "~", rownames(lambda))

  std_beta <- apply(boot_beta, 2, sd)
  t_ratio <- beta[beta!=0]/apply(boot_beta, 2, sd)
  pval_beta <- sapply(seq_along(t_ratio),
                      function(x)
                           2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )

  # beta <- data.frame()
  if (length(beta) != 0){

    beta <- data.frame(est = beta,
                    se = std_beta,
                    z = t_ratio,
                    pvalue = pval_beta)

    rownames(beta) <- sapply(seq_len(NROW(beta)),
                             function(b)
             paste(colnames(fit$beta)[which(fit$beta!=0, arr.ind = TRUE)[b, ]],
                   collapse = ".")
    )
    parts <- strsplit(rownames(beta), "\\.")
    beta <- cbind(
      data.frame(
        lhs = sapply(parts, `[`, 1),
        op = "~",
        rhs = sapply(parts, `[`, 2)
      ),    beta
    )

  }
  beta <- data.frame(beta)


  std_gamma <- apply(boot_gamma, 2, sd)

  t_ratio <- gamma[gamma!=0]/std_gamma
  pval_gamma <- sapply(seq_along(t_ratio),
                       function(x)
                            2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )

  gamma <- data.frame(est = gamma,
                     se = std_gamma,
                     z = t_ratio,
                     pvalue = pval_gamma)

  rownames(gamma) <- sapply(seq_len(NROW(gamma)),
                            function(b)
                            paste(rownames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 1]],
                                  colnames(fit$gamma)[which(fit$gamma!=0, arr.ind = TRUE)[b, 2]],
                                  sep = ".")
  )

  parts <- strsplit(rownames(gamma), "\\.")
  gamma <- cbind(
    data.frame(
      lhs = sapply(parts, `[`, 1),
      op = "~",
      rhs = sapply(parts, `[`, 2)
    ),
    gamma
  )



  std_total_effects <- apply(boot_total_effects, 2, sd)
  t_ratio <- total_effects/std_total_effects
  pval_total_effects <- sapply(seq_along(t_ratio),
                               function(x)
                                   2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )



  grid_total_effects <- expand.grid(
    LHS = rownames(fit$effect$total_effect), RHS = colnames(fit$effect$total_effect)
  )

  total_effects <- data.frame(
    lhs = grid_total_effects$LHS,
    op = "~",
    rhs = grid_total_effects$RHS,
    est = total_effects,
    se = std_total_effects,
    z = t_ratio,
    pvalue = pval_total_effects)

  rownames(total_effects) <- paste(grid_total_effects$LHS, grid_total_effects$RHS, sep = " ~ ")

  # rownames(total_effects) <- sapply(seq_len(NROW(total_effects)),
  #                         function(b)
  #                         paste(
  #                           rownames(fit$effect$total_effect)[which(fit$effect$total_effect!=0, arr.ind = TRUE)[b, 1]],
  #                           colnames(fit$effect$total_effect)[which(fit$effect$total_effect!=0, arr.ind = TRUE)[b, 2]],
  #                           sep = "~"
  #                         )
  # )



  std_indirect_effects <- apply(boot_indirect_effects, 2, sd)
  t_ratio <- indirect_effects/std_indirect_effects
  pval_indirect_effects <- sapply(seq_along(t_ratio),
                                  function(x)
                                      2*pnorm(abs(t_ratio[x]), lower.tail = FALSE)
  )



  grid_indirect_effects <- expand.grid(
    LHS = rownames(fit$effect$indirect_effect), RHS = colnames(fit$effect$indirect_effect)
  )


  indirect_effects <- data.frame(
    lhs = grid_indirect_effects$LHS,
    op = "~",
    rhs = grid_indirect_effects$RHS,
    est = indirect_effects,
    se = std_indirect_effects,
    z = t_ratio,
    pvalue = pval_indirect_effects)

  rownames(indirect_effects) <- paste(grid_indirect_effects$LHS, grid_indirect_effects$RHS, sep = " ~ ")


  # rownames(indirect_effects) <- sapply(
  #   seq_len(NROW(indirect_effects)),
  #   function(b)
  #     paste(
  #       rownames(fit$effect$indirect_effect)[which(fit$effect$indirect_effect!=0, arr.ind = TRUE)[b, 1]],
  #       colnames(fit$effect$indirect_effect)[which(fit$effect$indirect_effect!=0, arr.ind = TRUE)[b, 2]],
  #       sep = "~"
  #     )
  # )



  return(list(
    infer = list(
      out = list(
        boot_lambda,
        boot_std_loadings,
        boot_beta,
        boot_gamma,
        boot_residual_variance,
        boot_total_effects,
        boot_indirect_effects,
        boot_omega
      ),
      lambda = lambda,
      std_lambda = std_lambda,
      beta = beta,
      gamma = gamma,
      residual_variance = residual_variance,
      total_effects = total_effects,
      indirect_effects = indirect_effects,
      omega = omega,
      improper = sum(is.na(L[1, ]))
    ),
    gof = list(T_LS = fit$T_LS,
               Tb_LS = boot_Tb_LS,
               pval = mean(fit$T_LS <= boot_Tb_LS, na.rm = TRUE),
               improper = sum(is.na(boot_Tb_LS))
    )
  ))
}
