
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


correlation_duplication_matrix <- function(n, D_n = duplication_matrix(n)) {

  # we remove the columns corresponding to the diagonal elements in vech() to get the correlation duplication matrix
  # Direct computation of the indices of the diagonal elements in vech() for a matrix of size n x n
  j <- 1:n
  indices_diag <- 1 + (j - 1) * n - (j - 1) * (j - 2) / 2
  return(D_n[, -indices_diag, drop = FALSE])
}



commutation_matrix <- function(m, n = m) {
  mn <- m * n
  index_col <- as.vector(t(matrix(1:mn, nrow = m, ncol = n)))
  K_mn <- sparseMatrix(
    i = 1:mn,
    j = index_col,
    x = 1,
    dims = c(mn, mn)
  )
  return(K_mn)
}

diagonal_extraction_matrix <- function(m) {
  L_D <- sparseMatrix(
    i = 1:m,
    j = seq(from = 1, to = m^2, by = m + 1),
    x = 1,
    dims = c(m, m^2)
  )
  return(L_D)
}

elimination_matrix<- function(m) {

  nb_vech <- m * (m + 1) / 2
  mat_index <- matrix(1:(m^2), nrow = m, ncol = m)
  idx_keep <- mat_index[lower.tri(mat_index, diag = TRUE)]
  L_m <- sparseMatrix(
    i = 1:nb_vech,
    j = idx_keep,
    x = 1,
    dims = c(nb_vech, m^2)
  )
  return(L_m)
}
