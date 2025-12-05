
#' Chi-Square Test for Structural Equation Models
#'
#' This function computes the chi-square goodness-of-fit test statistic for
#' structural equation models, along with its degrees of freedom and p-value.
#'
#' @param p Number of observed variables in the model
#' @param q Number of estimated parameters in the model
#' @param r Number of formative blocks
#' @param F The fit function value (discrepancy between observed and implied covariance matrices)
#' @param N Sample size
#'
#' @return A list containing:
#'   \item{test}{The chi-square test statistic}
#'   \item{df}{Degrees of freedom for the test}
#'   \item{pval}{The p-value of the chi-square test}
#'
#' @details The chi-square statistic is calculated as (N-1) * F, where F is the
#' fit function value. The degrees of freedom are computed as p(p+1)/2 - q, representing
#' the difference between the number of unique elements in the covariance matrix
#' and the number of estimated parameters. A significant p-value (typically < 0.05)
#' suggests poor model fit.
#'
#' @examples
#' \dontrun{
#' chi2_results <- chi2sem(p = 10, q = 15, F = 0.25, N = 200)
#' print(chi2_results$pval)
#' }
#'
#' @export
chi2sem <- function(p, q, r, F, N){
  chi2 <- (N-1) * F
  df <- (p * (p+1)/2) - q + r
  p_value <- 1 - pchisq(chi2, df)

  return(list(test = chi2, df = df, pval = p_value))
}