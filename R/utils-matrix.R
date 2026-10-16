
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








calcul_D_vech_diag <- function(p_j) {

  # Indices des éléments de la diagonale dans l'espace vech
  k <- 1:p_j
  indices_diag <- 1 + (k - 1) * p_j - (k - 1) * (k - 2) / 2

  # Nombre total de lignes de la matrice (taille du vech)
  n_lignes <- p_j * (p_j + 1) / 2

  # Création de la matrice creuse
  # On place des "1" aux lignes correspondant à la diagonale,
  # et aux colonnes 1 à p_j (qui correspondent aux paramètres de variance)
  D_vech_diag_pj <- sparseMatrix(i = indices_diag,
                                 j = 1:p_j,
                                 x = 1,
                                 dims = c(n_lignes, p_j))

  return(D_vech_diag_pj)
}



selection_block <- function(block_sizes, j){

  P <- sum(block_sizes)
  p_j <- block_sizes[j]
  index_end <- cumsum(block_sizes)
  index_start <- c(1, head(index_end, -1) + 1)

  U_j <- sparseMatrix(i = index_start[j]:index_end[j],
                        j = 1:p_j,
                        x = 1,
                        dims = c(P, p_j))


  return(U_j)
}





generate_Pj <- function(block_sizes) {

  P_total <- sum(block_sizes)
  index_end <- cumsum(block_sizes)
  index_start <- c(1, index_end[-length(index_end)] + 1)
  I_P <- Diagonal(P_total)
  list_Pj <- lapply(seq_along(block_sizes), function(j) {
    I_P[, index_start[j]:index_end[j], drop = FALSE]
  })
  return(list_Pj)
}





P_exo_endo <- function(ind_exo, ind_endo) {

  J <- length(ind_exo) + length(ind_endo)

  I_J <- Diagonal(J)

  P_exo <- I_J[, ind_exo, drop = FALSE]
  P_endo <- I_J[, ind_endo, drop = FALSE]
  return(list(P_exo = P_exo, P_endo = P_endo))
}


M_beta_gamma <- function(C) {

  which_exo_endo <- ind_exo_endo(C)
  H <- which_exo_endo$ind_exo
  J <- which_exo_endo$ind_endo
  s_gamma <- as.vector(t(C[H, J, drop = FALSE]))
  s_beta  <- as.vector(t(C[J, J, drop = FALSE]))
  M_gamma <- Diagonal(length(s_gamma))[, s_gamma == 1, drop = FALSE]
  M_beta  <- Diagonal(length(s_beta))[, s_beta == 1, drop = FALSE]

  return(list(M_gamma = M_gamma, M_beta = M_beta))
}




correlation_elimination_matrix <- function(n, L_n = elimination_matrix(n)) {

  # we remove the columns corresponding to the diagonal elements in vech() to get the correlation duplication matrix
  # Direct computation of the indices of the diagonal elements in vech() for a matrix of size n x n
  j <- 1:n
  indices_diag <- 1 + (j - 1) * n - (j - 1) * (j - 2) / 2
  return(L_n[ -indices_diag, , drop = FALSE])
}


extraction_diagonale <- function(n) {

  # Les indices de la diagonale dans un vecteur de taille m_endo^2
  indices_diag <- 1 + (0:(n - 1)) * (n + 1)

  # D_diag : Matrice creuse de dimension (m_endo^2) x (m_endo)
  D_diag <- sparseMatrix(i = indices_diag,
                         j = 1:n,
                         x = 1,
                         dims = c(n^2, n))


  return(D_diag)
}





















#' Least Squares Discrepancy Function
#'
#' Computes the Euclidean distance between two covariance matrices using the
#' least squares discrepancy function.
#'
#' @param A A matrix, typically the observed covariance matrix
#' @param B A matrix, typically the model-implied covariance matrix
#'
#' @return The least squares discrepancy value representing the
#'   Euclidean distance between matrices A and B
#'
#' @details The LS discrepancy function is defined as:
#'   \deqn{d_{LS} = \sqrt{frac{1}{2} \times \text{tr}[(B-A)^2]}}
#'   where tr denotes the trace operator. This measures the overall fit
#'   between the observed and model-implied covariance matrices.
#'
#' @examples
#' \dontrun{
#' S <- diag(2)
#' SIGMA <- diag(2) * 0.5
#' discrepancy <- d_LS(S, SIGMA)
#'}
#'
#' @keywords internal
d_LS <- function(A, B){
  A <- as.matrix(A)
  B <- as.matrix(B)
  d_LS <- sqrt(0.5*sum(diag((B-A)%*%(B-A))))
  return(d_LS = d_LS)
}



#' Compute Reliability for Formative Blocks
#'
#' This function calculates the reliability of formative measurement blocks
#' using the specified metric. Currently, only the "Dillon" metric is supported.
#'
#' @param metric A character string specifying the reliability metric to use.
#'   Currently, only "Dillon" is implemented.
#' @param lambdas A list of factor loadings for each block.
#' @param residual_variances A list of residual variances for each block.
#'
#' @return A named numeric vector containing the reliability values for each block.
#'
#' @details The Dillon reliability metric is calculated as the squared sum of
#'   loadings divided by the sum of the squared loadings and the residual variances.
#'
#' @examples
#' \dontrun{
#' lambdas <- list(block1 = c(0.8, 0.7), block2 = c(0.9, 0.6, 0.5))
#' residuals <- list(block1 = c(0.2, 0.3), block2 = c(0.1, 0.4, 0.5))
#' reliability_values <- reliability(
#'  metric = "Dillon",
#'  lambdas = lambdas,
#'  residual_variances = residuals
#' )
#' print(reliability_values)
#' }
#' @keywords internal

# pour bloc formatif
reliability <- function (metric, lambdas, residual_variances){
  if (metric=='Dillon'){
    # compute sum of loadings and residual variances for each block
    sum_loadings <- sapply(lambdas, function(x) sum(as.numeric(x)))
    sum_residuals <- sapply(residual_variances, function(x) sum(as.numeric(x)))

    #remove formative blocks (blocks without residual variances)
    sum_loadings <- sum_loadings[names(sum_loadings) %in% names(sum_residuals)]

    # compute numerator and denominator for Dillon's rho
    numerator <- sum_loadings^2
    denominator <- numerator + sum_residuals

    # compute reliability results
    results <- numerator / denominator

    return(results)
  }
}

