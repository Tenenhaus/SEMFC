





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





P_exo_endo <- function(m_total, m_exo) {

  # 1. Création de la matrice identité globale (ultra-léger en mémoire)
  I_m <- Diagonal(m_total)

  # 2. P_exo : on prend simplement les m_exo premières colonnes
  P_exo <- I_m[, 1:m_exo, drop = FALSE]

  # 3. P_endo : on prend toutes les colonnes restantes
  P_endo <- I_m[, (m_exo + 1):m_total, drop = FALSE]

  return(list(P_exo = P_exo, P_endo = P_endo))
}


M_beta_gamma <- function(C) {

  which_exo_endo <- ind_exo_endo(C)
  H <- which_exo_endo$ind_exo
  J <- which_exo_endo$ind_endo
  s_gamma <- as.vector(t(C[H, J, drop = FALSE]))
  s_beta  <- as.vector(t(C[J, J, drop = FALSE]))
  M_gamma <- diag(length(s_gamma))[, s_gamma == 1, drop = FALSE]
  M_beta  <- diag(length(s_beta))[, s_beta == 1, drop = FALSE]

  return(list(M_gamma = M_gamma, M_beta = M_beta))
}




correlation_elimination_matrix <- function(n, L_n = elimination_matrix(n)) {

  # we remove the columns corresponding to the diagonal elements in vech() to get the correlation duplication matrix
  # Direct computation of the indices of the diagonal elements in vech() for a matrix of size n x n
  j <- 1:n
  indices_diag <- 1 + (j - 1) * n - (j - 1) * (j - 2) / 2
  return(L_n[ -indices_diag, , drop = FALSE])
}



