
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
#'   \item{p_close_fit}{The p-value for the close fit hypothesis (RMSEA ≤ close).}
#'   \item{p_notclose_fit}{The p-value for the not-close fit hypothesis (RMSEA ≥ notclose).}
#'
#' @details The RMSEA is calculated as:
#'   \deqn{\sqrt{\max(\chi^2 / df - 1, 0) / (N - 1)}}
#'   Confidence intervals are computed using the non-central chi-square distribution.
#'   P-values for close and not-close fit hypotheses are derived from the non-centrality
#'   parameter of the chi-square distribution.
#'
#' @examples
#' result <- rmseasem(chi2 = 20, df = 10, N = 200)
#' print(result$estimate)
#' @importFrom stats pchisq uniroot
#' @export
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

  # P-value close fit (RMSEA ≤ close)
  lambda_c <- df * (N - 1) * (close^2)
  p_close_fit <- pchisq(chi2, df = df, ncp = lambda_c, lower.tail = FALSE)
  # P-value not close fit (RMSEA ≥ notclose)
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
