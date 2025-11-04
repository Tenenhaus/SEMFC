

chi2sem <- function(p, q, F, N){
  chi2 <- (N-1) * F
  df <- (p * (p+1)/2) - q
  p_value <- 1 - pchisq(chi2, df)

  return(list(test = chi2, df = df, pval = p_value))
}