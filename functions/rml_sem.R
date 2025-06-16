

information_matrix <- function(x, block_sizes, C, mode){
  JAC <- numDeriv::jacobian(s_implied_bis, x = x, block_sizes = block_sizes, C = C, mode =mode ,jac = TRUE)
  full_jac <- sapply(1:NCOL(JAC),
                        function(col){
                          ds_dt <- matrix(0, sum(block_sizes), sum(block_sizes))
                          ds_dt[upper.tri(ds_dt, diag = T)] <- JAC[, col]
                          ds_dt <- ds_dt + t(ds_dt) - diag(diag(ds_dt))
                        }, simplify = FALSE
      )

  nb_param <- length(x)
  Sinv <- solve(s_implied_bis( x = x, block_sizes = block_sizes, C = C, mode =mode, jac = F)$SIGMA_IMPLIED)
  Sinv_full_jac <- lapply(full_jac, function(fj) as.matrix(Sinv %*% fj))
  I_ij <- function(i, j) {
    0.5 * sum(diag(Sinv_full_jac[[i]] %*% Sinv_full_jac[[j]]))
  }
  # Generate only indices of the upper triangular part
  index_upper <- which(upper.tri(matrix(0, nrow = nb_param, ncol = nb_param), diag = TRUE), arr.ind = TRUE)

  # Calculate only the elements of the upper triangular part
  I_upper <- mapply(I_ij, index_upper[, "row"], index_upper[, "col"])

  I <- matrix(0, nrow = nb_param, ncol = nb_param)
  I[upper.tri(I, diag = TRUE)] <- I_upper
  I <- I + t(I) - diag(diag(I)) 

  return(I)
}

Jac_constraints <- function(x, S, block_sizes, C, mode){
  # transpose of the jacobian of constraint function
  H <- t(numDeriv::jacobian(heq1_bis, x = x, S=S, block_sizes=block_sizes, C=C, mode = mode))

  return(H)

}


P_ml <- function(x, S, block_sizes, C, mode){

  I <- information_matrix(x, block_sizes, C, mode)
  H  <- Jac_constraints(x, S, block_sizes, C, mode)
  r  <- ncol(H)
  t <- nrow(H)
  M <- rbind(cbind(I+H%*%t(H), H),
             cbind(t(H), matrix(0, r, r)))
  P <- solve(M)[1:t, 1:t]

  return(P)


}


z_H0 <- function(x, S, X, C, mode, L){
  N <- nrow(X)
  P <- P_ml(x, S, X, C, mode)
  z <- sqrt(N)*(L %*% x)/sqrt(L%*%P%*%t(L))
  pval <- 2*pnorm(z, lower.tail = F)

  return(pval)

}









