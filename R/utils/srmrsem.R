

srmrsem <- function(S, Sigma) {
  R_emp <- cov2cor(S)
  R_mod <- cov2cor(Sigma)
  resids <- R_emp - R_mod
  SRMR <- sqrt(mean(resids[lower.tri(resids, diag=TRUE)]^2))
  return(SRMR)
}

