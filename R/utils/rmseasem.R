

rmseasem <- function(chi2, df, N, conf.level = 0.90, close = 0.05, notclose = 0.08) {
  RMSEA <- sqrt(max(chi2/df - 1, 0) / (N-1))

  alpha <- 1 - conf.level

  # CI
  find_ncp <- function(p_target) {
    uniroot(
      function(lambda) pchisq(chi2, df = df, ncp = lambda) - p_target,
      lower = 0, upper = 1e5
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
