



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



jac_f_s <- function(block_sizes, j, S, lambda){

  select_matrix <- selection_block(block_sizes, j)
  U_j <- select_matrix[[1]]
  P_j <- select_matrix[[2]]


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


jac_f_vechs <- function(block_sizes, j, S, lambda, D_P = duplication_matrix(sum(block_sizes))){

  jac_f_s <- jac_f_s(block_sizes, j, S, lambda)

  jac_f_vechs <- jac_f_s %*% D_P

  return(jac_f_vechs)
}



vcov_vechs <- function(Sigma, D_P){

  D_P_plus <- solve(t(D_P) %*% D_P) %*% t(D_P)

  return(2*D_P_plus %*% kronecker(Sigma, Sigma) %*% t(D_P_plus))

}


vcov_svd <- function(block_sizes, j, lambda, S, Sigma, D_P = duplication_matrix(sum(block_sizes))){

  J <- jac_f_vechs(block_sizes, j, S, lambda, D_P)
  V_vechs <- vcov_vechs(Sigma, D_P)

  return(J %*% V_vechs %*% t(J))

}