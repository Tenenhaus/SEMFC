

#' Get Significance Stars
#'
#' Assigns significance stars based on the p-value of a statistical test.
#'
#' @param p Numeric value representing the p-value.
#'   - If `p` is `NA`, returns an empty string.
#'   - If `p < 0.001`, returns `"***"`.
#'   - If `p < 0.01`, returns `"**"`.
#'   - If `p < 0.05`, returns `"*"`.
#'   - If `p < 0.1`, returns `"."`.
#'   - Otherwise, returns an empty string.
#'
#' @return A character string representing the significance level.
#'
#' @examples
#' \dontrun{
#' get_stars(0.0005) # Returns "***"
#' }
#' @keywords internal
get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  return("")
}




#' Print Parameter Estimates by Latent Variable
#'
#' Formats and prints parameter estimates for a given data frame, grouped by latent variables.
#'
#' @param df Data frame containing parameter estimates. Expected columns include:
#'   - `lhs`: Left-hand side (latent variable names).
#'   - `rhs`: Right-hand side (observed variable names or predictors).
#'   - `est`: Estimated parameter values.
#'   - `se`: Standard errors of the estimates (optional).
#'   - `z`: Z-scores (optional).
#'   - `pvalue`: P-values (optional).
#'   - `std.all`: Standardized estimates (optional, used if `standardized = TRUE`).
#' @param type Character string indicating the type of parameter (e.g., "Loadings", "Regression").
#' @param op Character string representing the operator (e.g., "=~", "~").
#' @param standardized Logical indicating whether to include standardized estimates in the output (default: `FALSE`).
#'
#' @return Invisibly returns `NULL`. Outputs formatted parameter estimates to the console.
#'
#' @keywords internal
print_link <- function(df, type, op, standardized = FALSE) {

  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))


  cat(type, ":\n")
  if (standardized) {
    cat(sprintf("%-14s %8s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)", "Std.all"))
  } else {
    cat(sprintf("%-14s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))
  }

  latents <- unique(df$lhs)

  for (lat in latents) {

    cat(paste0("  ", lat, " ", op, "\n"))


    group <- df[df$lhs == lat, ]

    for (i in seq_len(nrow(group))) {
      row <- group[i, ]


      rhs <- row$rhs


      est <- sprintf("%.3f", row$est)


      if (is.na(row$se) || row$se == 0) {
        se <- ""
        z_val <- ""
        p_val <- ""
        std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
        stars <- ""
      } else {
        se <- sprintf("%.3f", row$se)
        z_val <- sprintf("%.3f", row$z)
        p_val <- sprintf("%.3f", row$pvalue)
        std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
        stars <- get_stars(row$pvalue)
      }
      if (standardized) {
        cat(sprintf("    %-10s %8s %8s %8s %8s %8s %-3s\n",
                    row$rhs, est, se, z_val, p_val, std_all, stars))
      } else {
        cat(sprintf("    %-10s %8s %8s %8s %8s %-3s\n",
                    row$rhs, est, se, z_val, p_val, stars))
      }
    }
  }
}



#' Print Variance Estimates
#'
#' Formats and prints variance estimates for a given data frame.
#'
#' @param df Data frame containing variance estimates. Expected columns include:
#'   - `lhs`: Names of the variables.
#'   - `est`: Estimated variance values.
#'   - `se`: Standard errors of the estimates (optional).
#'   - `z`: Z-scores (optional).
#'   - `pvalue`: P-values (optional).
#'   - `std.all`: Standardized estimates (optional, used if `standardized = TRUE`).
#' @param type Character string indicating the type of variance (e.g., "Residual Variances").
#' @param standardized Logical indicating whether to include standardized estimates in the output (default: `FALSE`).
#'
#' @return Invisibly returns `NULL`. Outputs formatted variance estimates to the console.
#'
#' @keywords internal
print_variance <- function(df, type, standardized = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  cat("\n", type, ":\n")
  if (standardized) {
    cat(sprintf("%-14s %8s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)", "Std.all"))
  } else {
    cat(sprintf("%-14s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))
  }

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    var_name <- row$lhs
    display_name <- paste0("   .", var_name)

    est <- sprintf("%.3f", row$est)


    if (is.na(row$se) || row$se == 0) {
      se <- ""
      z_val <- ""
      p_val <- ""
      std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
      stars <- ""
    } else {
      se <- sprintf("%.3f", row$se)
      z_val <- sprintf("%.3f", row$z)
      p_val <- sprintf("%.3f", row$pvalue)
      std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
      stars <- get_stars(row$pvalue)
    }

    if (standardized) {
      cat(sprintf("%-14s %8s %8s %8s %8s %8s %-3s\n",
                  display_name, est, se, z_val, p_val, std_all, stars))
    } else {
      cat(sprintf("%-14s %8s %8s %8s %8s %-3s\n",
                  display_name, est, se, z_val, p_val, stars))
    }
  }
}



#' Print Parameter Estimates
#'
#' Formats and prints parameter estimates, including loadings, regression coefficients,
#' residual variances, and effects, based on the provided model estimates.
#'
#' @param estimate List containing the model estimates. Expected elements include:
#'   - `lambda`: Data frame of loadings.
#'   - `residual_variance`: Data frame of residual variances.
#'   - `beta`: Data frame of regression coefficients for endogenous variables.
#'   - `gamma`: Data frame of regression coefficients for exogenous variables.
#'   - `total_effects`: Data frame of total effects (optional, used if `effect = TRUE`).
#'   - `indirect_effects`: Data frame of indirect effects (optional, used if `effect = TRUE`).
#'   - `omega`: Data frame of formative block weights (optional).
#' @param standardized Logical indicating whether to include standardized estimates in the output.
#' @param effect Logical indicating whether to print total and indirect effects.
#'
#' @return Invisibly returns `NULL`. Outputs formatted parameter estimates to the console.
#'
#' @keywords internal

print_estimates <- function(estimate, standardized, effect){

      lambda <- estimate$lambda
      residualvariance <- estimate$residual_variance
      beta <- estimate$beta
      gamma <- estimate$gamma
      total_effects <- estimate$total_effects
      indirect_effects <- estimate$indirect_effects
      omega <- estimate$omega

      cat("\nParameter Estimates:\n")
      print_link(lambda, 'Loadings', '=~', standardized)
      print_link(omega, 'Formative block weight', '<~', standardized)
      print_link(rbind(beta, gamma), 'Regression', '~', standardized)
      print_variance(residualvariance, 'Residual Variances', standardized)
      if (effect){
        print_link(total_effects, 'Total Effects', '~', standardized)
        print_link(indirect_effects, 'Indirect Effects', '~', standardized)
      }
      cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}


#' Print Model Information
#'
#' Formats and prints basic model information, including the estimator type,
#' number of parameters, sample size, degrees of freedom, and fit statistics.
#'
#' @param estimator Character string indicating the estimation method (e.g., "ml", "svd").
#'   The estimator name will be displayed in uppercase.
#' @param len_theta Integer specifying the number of model parameters.
#' @param N Integer specifying the number of observations in the sample.
#' @param dof Integer specifying the degrees of freedom for the model.
#' @param F Numeric value representing the F statistic.
#' @param T_LS Numeric value representing the d_LS (least squares discrepancy) statistic.
#'
#' @return Invisibly returns `NULL`. Outputs formatted model information to the console.
#'
#' @keywords internal
print_model <- function(estimator, len_theta, N, dof, F, T_LS) {
  cat("\n")
  cat(sprintf("%-45s%15s\n", "Estimator", toupper(estimator)))
  cat(sprintf("%-45s%15d\n", "Number of model parameters", len_theta))
  cat(sprintf("%-45s%15d\n", "Number of observations", N))
  cat(sprintf("%-45s%15d\n", "Degrees of freedom", dof))
  cat(sprintf("%-45s%15.3f\n", "F", F))
  cat(sprintf("%-45s%15.3f\n", "d_LS", T_LS))
  cat("\n")
}




#' Print Goodness-of-Fit Statistics for SVD Estimation
#'
#' Formats and prints the results of the Bollen-Stine bootstrap test
#' for models estimated using the SVD method.
#'
#' @param B Integer specifying the number of bootstrap replications performed.
#' @param pvalbs Numeric value representing the Bollen-Stine bootstrap p-value.
#'   This p-value tests the null hypothesis that the model fits the data.
#'
#' @return Invisibly returns `NULL`. Outputs formatted bootstrap test results to the console.
#'
#' @keywords internal

print_gof_svd <- function(B, pvalbs) {

  cat("Bootstrap Test (Bollen Stine):\n\n")
  cat(sprintf("  %-40s%12d\n", "Number of bootstrap replications", B))
  cat(sprintf("  %-40s%12.3f\n", "Bollen Stine bootstrap p-value", pvalbs))
  cat("\n")
}


#' Print Goodness-of-Fit Statistics for ML Estimation
#'
#' Formats and prints comprehensive goodness-of-fit statistics for models
#' estimated using the maximum likelihood (ML) method.
#'
#' @param gof List containing goodness-of-fit statistics. Expected elements include:
#'   - `chi2`: List with `test`, `df`, and `pval` for the user model chi-square test.
#'   - `baseline`: List with `test`, `df`, and `pval` for the baseline model test.
#'   - `cfi`: Numeric value for the Comparative Fit Index.
#'   - `tli`: Numeric value for the Tucker-Lewis Index.
#'   - `RMSEA`: List with `estimate`, `CI_lower`, `CI_upper`, `p_close_fit`, and `p_notclose_fit`.
#'   - `SRMR`: Numeric value for the Standardized Root Mean Square Residual.
#'   - `info_criteria`: List with `loglik_H0`, `loglik_H1`, `AIC`, `BIC`, and `SABIC`.
#'
#' @return Invisibly returns `NULL`. Outputs formatted goodness-of-fit statistics to the console,
#'   including model test statistics, fit indices (CFI, TLI), RMSEA with confidence intervals,
#'   SRMR, and information criteria (AIC, BIC, SABIC).
#'
#' @keywords internal
print_gof_ml <- function(gof) {

    # user test
    testchi2 <- gof$chi2$test
    dfchi2 <- gof$chi2$df
    pvalchi2 <- gof$chi2$pval

    # baseline test

    testbaseline <- gof$baseline$test
    dfbaseline <- gof$baseline$df
    pvalbaseline <- gof$baseline$pval

    #  vs
    cfi <- gof$cfi
    tli <- gof$tli

    # RMSEA and srmr
    rmsea_val <- gof$RMSEA$estimate
    rmsea_ci_lower <- gof$RMSEA$CI_lower
    rmsea_ci_upper <- gof$RMSEA$CI_upper
    p_rmsea_le_005 <- gof$RMSEA$p_close_fit
    p_rmsea_ge_008 <- gof$RMSEA$p_notclose_fit
    srmr_val <- gof$SRMR
    cat("Model Test User Model :\n\n")
    cat(sprintf("  %-40s%12.3f\n", "Test statistic", testchi2))
    cat(sprintf("  %-40s%12d\n", "Degrees of freedom", dfchi2))
    cat(sprintf("  %-40s%12.3f\n", "P-value (Chi-square)", pvalchi2))
    cat("\n")

    cat("Model Test Baseline Model :\n\n")
    cat(sprintf("  %-40s%12.3f\n", "Test statistic", testbaseline))
    cat(sprintf("  %-40s%12d\n", "Degrees of freedom", dfbaseline))
    cat(sprintf("  %-40s%12.3f\n", "P-value", pvalbaseline))
    cat("\n")

    cat("User Model versus Baseline Model:\n\n")
    cat(sprintf("  %-40s%12.3f\n", "Comparative Fit Index (CFI)", cfi))
    cat(sprintf("  %-40s%12.3f\n", "Tucker-Lewis Index (TLI)", tli))
    cat("\n")

    cat("Root Mean Square Error of Approximation:\n\n")
    cat(sprintf("  %-40s%12.3f\n", "RMSEA", rmsea_val))
    cat(sprintf("  %-40s%12.3f\n", "90 Percent confidence interval - lower", rmsea_ci_lower))
    cat(sprintf("  %-40s%12.3f\n", "90 Percent confidence interval - upper", rmsea_ci_upper))
    cat(sprintf("  %-40s%12.3f\n", "P-value H_0: RMSEA <= 0.050", p_rmsea_le_005))
    cat(sprintf("  %-40s%12.3f\n", "P-value H_0: RMSEA >= 0.080", p_rmsea_ge_008))
    cat("\n")

    cat("Standardized Root Mean Square Residual:\n\n")
    cat(sprintf("  %-40s%12.3f\n", "SRMR", srmr_val))
    cat("\n")

    # Information criteria

    loglik_H0 <-  gof$info_criteria$loglik_H0
    loglik_H1 <-  gof$info_criteria$loglik_H1
    AIC <-  gof$info_criteria$AIC
    BIC <-  gof$info_criteria$BIC
    SABIC <-  gof$info_criteria$SABIC

    cat("Loglikelihood and Information Criteria:\n\n")
    cat(sprintf("  %-40s%12.3f\n", "Loglikelihood user model (H0)", loglik_H0))
    cat(sprintf("  %-40s%12.3f\n", "Loglikelihood unrestricted model (H1)", loglik_H1))
    cat("\n")
    cat(sprintf("  %-40s%12.3f\n", "Akaike (AIC)", AIC))
    cat(sprintf("  %-40s%12.3f\n", "Bayesian (BIC)", BIC))
    cat(sprintf("  %-40s%12.3f\n", "Sample-size adjusted BIC (SABIC)", SABIC))
    cat("\n")
}



#' Print Goodness-of-Fit Measures
#'
#' Formats and prints goodness-of-fit measures, reliability coefficients, and R-squared values
#' based on the estimation method and available statistics.
#'
#' @param all_measures Logical indicating whether to print comprehensive goodness-of-fit statistics.
#'   If `TRUE`, calls the appropriate GOF function based on the estimator type.
#' @param estimator Character string indicating the estimation method (e.g., "ml", "svd").
#'   - If `"ml"`, prints ML-specific goodness-of-fit statistics.
#'   - If `"svd"`, prints Bollen-Stine bootstrap test results if available.
#' @param gof List containing goodness-of-fit statistics. Expected elements depend on the estimator:
#'   - For ML: See `print_gof_ml()` documentation.
#'   - For SVD: `bollen_stine` list with `pval`.
#'   - `reliability`: Named vector of reliability coefficients (optional).
#' @param B Integer specifying the number of bootstrap replications (used for SVD estimator).
#' @param R2 Named vector of R-squared values for endogenous variables (optional).
#'
#' @return Invisibly returns `NULL`. Outputs formatted goodness-of-fit measures to the console,
#'   including model-specific fit statistics, reliability coefficients (Dillon), and R-squared values.
#'
#' @keywords internal
print_gof <- function(all_measures, estimator, gof, B, R2) {

    if (all_measures){
      if (estimator == 'ml'){
        print_gof_ml(gof)
      }
      if (estimator == 'svd' && !is.null(gof$bollen_stine)){
        print_gof_svd(B, gof$bollen_stine$pval)
      }
    }

    # reliability

    if (!is.null(gof$reliability)){
      cat("Reliability Coefficients (Dillon):\n")
      print(round(gof$reliability, 3))
      cat("\n")
    }

    # R2

    if (!is.null(R2)){
        cat("R2:\n")
        print(round(R2, 3))
        cat("\n")
    }

}

