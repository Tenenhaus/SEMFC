#' Chi-Square Test for Structural Equation Models
#'
#' This function computes the chi-square goodness-of-fit test statistic for
#' structural equation models, along with its degrees of freedom and p-value.
#'
#' @param p Number of observed variables in the model
#' @param q Number of estimated parameters in the model
#' @param r Number of formative blocks
#' @param F The fit function value
#' @param N Sample size
#'
#' @return A list containing:
#'   \item{test}{The chi-square test statistic}
#'   \item{df}{Degrees of freedom for the test}
#'   \item{pval}{The p-value of the chi-square test}
#'
#' @details The chi-square statistic is calculated as (N-1) * F, where F is the
#' fit function value. The degrees of freedom are computed as p(p+1)/2 - q + r
#'
#' @examples
#' \dontrun{
#' chi2_results <- chi2sem(p = 18, q = 61, r = 4, F = 0.428, N = 300)
#' print(chi2_results$pval)
#'}
#' @importFrom stats pchisq
#' @keywords internal
chi2sem <- function(p, q, r, F, N){
  chi2 <- (N-1) * F
  df <- (p * (p+1)/2) - q + r
  p_value <- 1 - pchisq(chi2, df)

  return(list(test = chi2, df = df, pval = p_value))
}

#' Standardized Root Mean Square Residual (SRMR) for Structural Equation Models
#'
#' This function calculates the SRMR, a goodness-of-fit measure for structural equation models,
#' by comparing the empirical covariance matrix with the model-implied covariance matrix.
#'
#' @param S A covariance matrix representing the empirical data.
#' @param Sigma A covariance matrix representing the model-implied data.
#'
#' @return A numeric value representing the SRMR, which is the square root of the mean
#'   of the squared residuals between the lower triangular elements (including the diagonal)
#'   of the standardized empirical and model-implied correlation matrices.
#'
#' @details The function first standardizes the covariance matrices to correlation matrices
#'   using \code{cov2cor}. It then computes the residuals as the difference between the
#'   empirical and model-implied correlation matrices. The SRMR is calculated as the square
#'   root of the mean of the squared residuals from the lower triangular part of the matrix.
#'
#' @examples
#' \dontrun{
#' S <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' Sigma <- matrix(c(1, 0.4, 0.4, 1), ncol = 2)
#' srmr <- srmrsem(S, Sigma)
#' print(srmr)
#'}
#' @keywords internal
srmrsem <- function(S, Sigma) {
  R_emp <- cov2cor(S)
  R_mod <- cov2cor(Sigma)
  resids <- R_emp - R_mod
  SRMR <- sqrt(mean(resids[lower.tri(resids, diag=TRUE)]^2))
  return(SRMR)
}



#' Root Mean Square Error of Approximation (RMSEA) for Structural Equation Models
#'
#' This function calculates the RMSEA for a structural equation model, along with its
#' confidence interval and p-values for close and not-close fit hypotheses.
#'
#' @param chi2 The chi-square test statistic of the model.
#' @param df Degrees of freedom associated with the chi-square test.
#' @param N Sample size used in the model.
#' @param conf.level Confidence level for the confidence interval (default is 0.90).
#' @param close Threshold for close fit hypothesis (default is 0.05).
#' @param notclose Threshold for not-close fit hypothesis (default is 0.08).
#'
#' @return A list containing:
#'   \item{estimate}{The RMSEA estimate.}
#'   \item{CI_lower}{The lower bound of the confidence interval for RMSEA.}
#'   \item{CI_upper}{The upper bound of the confidence interval for RMSEA.}
#'   \item{p_close_fit}{The p-value for the close fit hypothesis (RMSEA <= close).}
#'   \item{p_notclose_fit}{The p-value for the not-close fit hypothesis (RMSEA >= notclose).}
#'
#' @details The RMSEA is calculated as:
#'   \deqn{\sqrt{\max(\chi^2 / df - 1, 0) / (N - 1)}}
#'   Confidence intervals are computed using the non-central chi-square distribution.
#'   P-values for close and not-close fit hypotheses are derived from the non-centrality
#'   parameter of the chi-square distribution.
#'
#' @examples
#' \dontrun{
#' result <- rmseasem(chi2 = 20, df = 10, N = 200)
#' print(result$estimate)
#' }
#' @importFrom stats pchisq uniroot
#' @keywords internal
rmseasem <- function(chi2, df, N, conf.level = 0.90, close = 0.05, notclose = 0.08) {
  RMSEA <- sqrt(max(chi2/df - 1, 0) / (N-1))

  alpha <- 1 - conf.level

  # CI
  find_ncp <- function(p_target) {
    f0 <- pchisq(chi2, df = df, ncp = 0) - p_target

    if (f0 <= 0) return(0)


    upper <- 1e5
    f_upper <- pchisq(chi2, df = df, ncp = upper) - p_target
    max_upper <- 1e10
    while (f_upper > 0 && upper < max_upper) {
      upper <- upper * 10
      f_upper <- pchisq(chi2, df = df, ncp = upper) - p_target
    }
    if (f_upper > 0) return(NA)

    uniroot(
      function(lambda) pchisq(chi2, df = df, ncp = lambda) - p_target,
      lower = 0, upper = upper
    )$root
  }
  lambda.lower <- find_ncp(1 - alpha / 2)
  lambda.upper <- find_ncp(alpha / 2)

  lower <- sqrt(lambda.lower / (df * (N - 1)))
  upper <- sqrt(lambda.upper / (df * (N - 1)))

  # P-value close fit (RMSEA <= close)
  lambda_c <- df * (N - 1) * (close^2)
  p_close_fit <- pchisq(chi2, df = df, ncp = lambda_c, lower.tail = FALSE)
  # P-value not close fit (RMSEA >= notclose)
  lambda_nc <- df * (N - 1) * (notclose^2)
  p_notclose <- pchisq(chi2, df = df, ncp = lambda_nc, lower.tail = TRUE)

  return(
    list(
      estimate         = RMSEA,
      CI_lower      = lower,
      CI_upper      = upper,
      p_close_fit   = p_close_fit,
      p_notclose_fit= p_notclose
    )
  )
}


#' Calculate Goodness-of-Fit Statistics for SEM-ML Estimation
#'
#' @description
#' Computes a comprehensive set of goodness-of-fit indices for structural equation
#' models estimated using Maximum Likelihood (ML) method.
#'
#' @param p Integer. Number of observed variables
#' @param q Integer. Number of free parameters in the model
#' @param r Integer. Number of latent variables
#' @param F Numeric. Value of the fitting function at convergence
#' @param N Integer. Sample size
#' @param S Matrix. Sample covariance matrix
#' @param Sigma Matrix. Model-implied covariance matrix
#'
#' @details
#' This function calculates multiple fit indices commonly used in SEM:
#'
#' * **Chi-square test**: Tests the null hypothesis that the model fits perfectly
#' * **Baseline model**: Independence model where all variables are uncorrelated
#' * **CFI** (Comparative Fit Index): Compares target model to baseline (range: 0-1, >0.95 good)
#' * **TLI** (Tucker-Lewis Index): Non-normed fit index penalizing model complexity
#' * **RMSEA** (Root Mean Square Error of Approximation): Measures approximate fit
#' * **SRMR** (Standardized Root Mean Square Residual): Average residual correlation
#' * **Information criteria** (AIC, BIC, SABIC): For model comparison (lower is better)
#'
#' The baseline model assumes all observed variables are independent (diagonal covariance matrix).
#'
#' @return A list containing:
#' \describe{
#'   \item{chi2}{Chi-square test statistics and p-value}
#'   \item{baseline}{Baseline model chi-square test}
#'   \item{cfi}{Comparative Fit Index}
#'   \item{tli}{Tucker-Lewis Index}
#'   \item{RMSEA}{RMSEA value with confidence interval}
#'   \item{SRMR}{Standardized Root Mean Square Residual}
#'   \item{info_criteria}{List containing loglik_H0, loglik_H1, AIC, BIC, SABIC}
#' }
#'
#' @seealso
#' \code{\link{chi2sem}}, \code{\link{rmseasem}}, \code{\link{srmrsem}}
#'
#' @keywords internal
semML_gof <- function(p, q, r, F, N, S, Sigma) {
  chi2 <- chi2sem(p, q, r, F, N)

  # basline test
  S_baseline <- diag(diag(S))
  F_baseline <- log(det(S_baseline)) + sum(diag(S%*%solve(S_baseline))) - log(det(S)) - NCOL(S)
  baseline <- chi2sem(p, p, 0, F_baseline, N)


  cfi <- 1 - (max(chi2$test - chi2$df, 0)) /
       (max(baseline$test - baseline$df, chi2$test - chi2$df, 0))

  tli <- ( (baseline$test /  baseline$df) - (chi2$test / chi2$df) ) /
         ( (baseline$test /  baseline$df) - 1 )

  # rmsea
  RMSEA <- rmseasem(chi2$test, chi2$df, N)
  SRMR <- srmrsem(S, Sigma)

  # loglik
  loglik_H0 <- -(N/2)*(p*log(2*pi) + log(det(Sigma)) + sum(diag(solve(Sigma) %*% S)))
  if (det(S) == 0) {
    loglik_H1 <- NA
  } else {
    loglik_H1 <- -(N/2)*(p*log(2*pi) + log(det(S)) + sum(diag(solve(S) %*% S)))
  }
  AIC <- -2 * loglik_H0 + 2 * q
  BIC <- -2 * loglik_H0 + q * log(N)
  SABIC <- -2 * loglik_H0 + q * log((N + 2) / 24)

  return(
    list(
      chi2 = chi2,
      baseline = baseline,
      cfi = cfi,
      tli = tli,
      RMSEA = RMSEA,
      SRMR = SRMR,
      info_criteria = list(
        loglik_H0 = loglik_H0,
        loglik_H1 = loglik_H1,
        AIC = AIC,
        BIC = BIC,
        SABIC = SABIC
      )
    )
  )
}



#' Detect Improper Solutions in SEM Estimation
#'
#' Checks for nine types of inadmissible or improper solutions in structural
#' equation models, including out-of-range estimates and non-positive definite
#' matrices.
#'
#' @param fit A list of estimates, including:
#' \itemize{
#'   \item \code{reliability_coef}: Reliability coefficients for each block
#'   \item \code{p_implied}: Implied correlation matrix of latent variables
#'   \item \code{p_tilde}: Corrected correlation matrix of latent variables
#'   \item \code{residual_variance}: List of residual variances
#'   \item \code{std_lambda}: Standardized loadings
#'   \item \code{sigma_implied}: Implied covariance matrix of observed variables
#'   \item \code{r2}: R-squared values for endogenous latent variables
#'   \item \code{psi}: Residual covariance matrix of latent variables
#' }
#'
#' @return Named logical vector of length 8 (or 9 for SVDSEM) indicating presence of each type
#'   of improper solution:
#' \describe{
#'   \item{reliability_coef}{\code{TRUE} if any reliability coefficient is outside (0, 1).}
#'   \item{rho_jh}{\code{TRUE} if any correlation in p_implied is outside (-1, 1).}
#'   \item{p_implied}{\code{TRUE} if p_implied has negative eigenvalues (not positive definite).}
#'   \item{theta_jh}{\code{TRUE} if any residual variance is negative.}
#'   \item{std_lambda}{\code{TRUE} if any standardized loading is outside (-1, 1).}
#'   \item{sigma_implied}{\code{TRUE} if sigma_implied has negative eigenvalues (not positive definite).}
#'   \item{r2}{\code{TRUE} if any R-squared is outside (0, 1).}
#'   \item{psi}{\code{TRUE} if psi has negative eigenvalues (not positive definite).}
#'   \item{p_tilde}{\code{TRUE} if p_tilde has negative eigenvalues (only for SVDSEM).}
#' }
#'
#' @details
#' A value of \code{TRUE} indicates an improper solution for that criterion. An
#' admissible solution should have all values set to \code{FALSE}. Common causes
#' of improper solutions include poor model specification, insufficient sample
#' size, or convergence to a boundary solution.
#'
#' @keywords internal
improper <- function(fit){


  improper_sol <- rep(FALSE, 8)

  # Improper solution 1 : number of out of range reliability coefficient
  improper_sol[1] <- any(fit$reliability_coef < 0 |
                        fit$reliability_coef > 1)
  # Improper solution 2 : number of out of range rho_jk
  improper_sol[2] <- any(abs(fit$p_implied[upper.tri(fit$p_implied)])>1)
  # Improper solution 3 : is p_implied positive definite
  improper_sol[3] <- any(eigen(fit$p_implied)$values<0)
  # Improper solution 4 : is variance of the residual positive
  improper_sol[4] <- any(unlist(fit$residual_variance)<0)
  # Improper solution 5 : number of out of range cor(y_jh, eta_j)
  improper_sol[5] <- any(abs(unlist(fit$std_lambda))>1)
  # Improper solution 6 : is sigma_implied positive definite
  improper_sol[6] <- any(eigen(fit$sigma_implied)$values<0)
  # Improper solution 7 : number of out of range R^2_i
  improper_sol[7] <- any(fit$r2< 0 | fit$r2>1)
  # Improper solution 8 : is psi positive definite
  improper_sol[8] <- any(eigen(fit$psi)$values<0)

  names(improper_sol) <- c("reliability_coef",
                          "rho_jh",
                          "p_implied",
                          "theta_jh",
                          "std_lambda",
                          "sigma_implied",
                          "r2",
                          "psi")

  # Improper solution 9 : is p_tilde positive definite
  if (!is.null(fit$p_tilde)){
    improper_sol[9] <- any(eigen(fit$p_tilde)$values<0)
    names(improper_sol)[9] <- "p_tilde"
  }

  return(improper_sol= improper_sol)
}