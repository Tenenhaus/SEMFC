# Function that tranforms the data sets in the way proposed
# by Yuan & Hayashi (2003)

scaleDataSet <- function(data, Sigma, bias = FALSE){
  S <- cov2(data, bias = bias)
  # singular value decomposition
  S1 = svd(S); S2 = svd(Sigma)
  d1 = S1$d; d2 = S2$d
  u1 = S1$u; u2 = S2$u
  v1 = S1$v; v2 = S2$v
  S_half <- u1%*%(diag(d1^(-1/2)))%*%t(v1) #S^(-1/2)
  Sigma_half <- u2%*%(diag(d2^(1/2)))%*%t(v2) # Sigma^(1/2)
  ScaleData <- as.matrix(data)%*%S_half%*%Sigma_half
  colnames(ScaleData) <- colnames(data)
  return(ScaleData = ScaleData)
}
