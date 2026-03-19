



selection_block <- function(block_sizes, j){

  P <- sum(block_sizes)
  p_j <- block_sizes[j]
  index_end <- cumsum(block_sizes)
  index_start <- c(1, head(index_end, -1) + 1)

  U_j <- sparseMatrix(i = index_start[j]:index_end[j],
                        j = 1:p_j,
                        x = 1,
                        dims = c(P, p_j))


  mask_diagonal <- rep(1, P)
  mask_diagonal[index_start[j]:index_end[j]] <- 0

  P_without_j <- Diagonal(x = mask_diagonal)


  return(c(U_j, P_without_j))
}




jac_u_s <- function(U_j, P_j, S, lambda){

  M_temp <- t(U_j)%*%S%*%P_j

  part1 <- (t(lambda)%*%M_temp) %x% t(U_j)
  part2 <- (t(lambda)%*%t(U_j)) %x% M_temp

  return(part1 + part2)
}


jac_f_u <- function(block_sizes, j, u_j){

  p_j <- block_sizes[j]

  norm2_u_j <- as.numeric(crossprod(u_j))
  norm_u_j <- sqrt(norm2_u_j)

  jac_f_u <- (Diagonal(p_j) - tcrossprod(u_j)/norm2_u_j)/norm_u_j

  return(jac_f_u)
}



jac_f_s <- function(block_sizes, j, S, lambda, U_j, P_j){

  # select_matrix <- selection_block(block_sizes, j)
  # U_j <- select_matrix[[1]]
  # P_j <- select_matrix[[2]]


  A_j <- t(U_j)%*%S%*%P_j%*%S%*%U_j
  u_j <- A_j %*% lambda


  jac_u_s <- jac_u_s(U_j, P_j, S, lambda)
  jac_f_u <- jac_f_u(block_sizes, j, u_j)

  jac_f_s <- jac_f_u%*%jac_u_s

  return(jac_f_s)
}


duplication_matrix <- function(P) {

  n_vech <- P * (P + 1) / 2
  mat_index <- matrix(0, nrow = P, ncol = P)
  mat_index[lower.tri(mat_index, diag = TRUE)] <- 1:n_vech
  mat_index[upper.tri(mat_index)] <- t(mat_index)[upper.tri(mat_index)]
  index_col <- as.vector(mat_index)
  D_P <- sparseMatrix(
    i = 1:(P^2),
    j = index_col,
    x = 1,
    dims = c(P^2, n_vech)
  )

  return(D_P)
}


jac_f_vechs <- function(block_sizes, j, S, lambda, U_j, P_j, D_P = duplication_matrix(sum(block_sizes))){

  jac_f_s <- jac_f_s(block_sizes, j, S, lambda, U_j, P_j)

  jac_f_vechs <- jac_f_s %*% D_P

  return(jac_f_vechs)
}



vcov_vechs <- function(Sigma, D_P){

  D_P_plus <- solve(t(D_P) %*% D_P) %*% t(D_P)

  return(2*D_P_plus %*% kronecker(Sigma, Sigma) %*% t(D_P_plus))

}


vcov_lambda_unit <- function(block_sizes, j, lambda, S, Sigma, U_j, P_j, D_P = duplication_matrix(sum(block_sizes))){

  J <- jac_f_vechs(block_sizes, j, S, lambda, U_j, P_j, D_P)
  V_vechs <- vcov_vechs(Sigma, D_P)

  return(J %*% V_vechs %*% t(J))

}


offdiag <- function(M){

  M_offdiag <- M
  diag(M_offdiag) <- 0

  return(M_offdiag)}


jac_lambda_reflective <- function(J_lambda_unit, lambda, S, U_j){

  # J_lambda_unit <- jac_f_vechs(block_sizes, j, S, lambda,U_j, P_j, D_P)

  S_jj_offdiag <- offdiag(t(U_j)%*%S%*%U_j)
  L_tilde <- offdiag(tcrossprod(lambda))

  A <- as.numeric(t(lambda)%*%S_jj_offdiag%*%lambda)
  B <- as.numeric(1 - as.numeric(crossprod(lambda^2)))
  d_j <- as.numeric(sqrt(A/B))

  v_T <- (1/A)*crossprod(lambda, S_jj_offdiag) + (2/B)*t(lambda^3)
  part1 <- J_lambda_unit + lambda %*% v_T %*% J_lambda_unit


  M <- U_j%*%L_tilde%*%t(U_j)
  vec_Dn_M <- M[lower.tri(M, diag = TRUE)]
  part2 <- (1/A)*lambda%*%t(vec_Dn_M)

  J_lambda <- d_j*(part1 + part2)

  return(J_lambda)
}


jac_lambda_formative <- function(J_lambda_unit, lambda, S, U_j ){

  # J_lambda_unit <- jac_f_vechs(block_sizes, j, S, lambda,U_j, P_j, D_P)
  S_inv_lambda <- solve(t(U_j)%*%S%*%U_j, lambda)
  u <- as.numeric(crossprod(lambda, S_inv_lambda))
  d_j <- 1 / sqrt(u)
  d_j3 <- d_j^3
  part1 <- (d_j * J_lambda_unit) - (d_j3 * lambda%*%crossprod(S_inv_lambda, J_lambda_unit))

  w <- U_j %*% S_inv_lambda
  H <- tcrossprod(w)
  diag_H <- diag(H)
  H <- H * 2
  diag(H) <- diag_H
  vec_Dn_H <- H[lower.tri(H, diag = TRUE)]
  part2 <- (d_j3 / 2) * tcrossprod(lambda, vec_Dn_H)

  J_lambda <- part1 + part2

  return(J_lambda)
}









