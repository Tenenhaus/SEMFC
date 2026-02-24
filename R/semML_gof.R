


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
  loglik_H1 <- -(N/2)*(p*log(2*pi) + log(det(S)) + sum(diag(solve(S) %*% S)))
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