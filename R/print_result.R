


get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  return("")
}


print_link <- function(df, type, op, standardized = FALSE) {

  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  # 2. Affichage de l'en-tête
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
    # Affichage du nom de la variable latente (ex: CUSTOMER_E =~)
    cat(paste0("  ", lat, " ", op, "\n"))

    # On prend les lignes correspondant à cette variable
    group <- df[df$lhs == lat, ]

    for (i in seq_len(nrow(group))) {
      row <- group[i, ]

      # Nom de l'indicateur (rhs)
      rhs <- row$rhs

      # Estimation formatée à 3 décimales
      est <- sprintf("%.3f", row$est)

      # Gestion des paramètres fixés : si SE est 0 ou NA, on laisse vide
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
    display_name <- paste0("   .", var_name) # Ajout du point
    # --- Formatage des chiffres ---
    est <- sprintf("%.3f", row$est)

    # Gestion des cas où SE est vide/zéro (paramètres fixés)
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





print_estimates <- function(estimate, standardized, effect){

      lambda <- estimate$lambda
      residualvariance <- estimate$residual_variance
      beta <- estimate$beta
      gamma <- estimate$gamma
      total_effects <- estimate$total_effects
      indirect_effects <- estimate$indirect_effects
      if (!is.null(estimate$omega)){
        omega <- estimate$omega
      }

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






print_gof_svd <- function(B, pvalbs) {

  cat("Bootstrap Test (Bollen Stine):\n\n")
  cat(sprintf("  %-40s%12d\n", "Number of bootstrap replications", B))
  cat(sprintf("  %-40s%12.3f\n", "Bollen Stine bootstrap p-value", pvalbs))
  cat("\n")
}



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

