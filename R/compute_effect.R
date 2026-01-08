#' Compute Total and Indirect Effects in a SEM Model
#'
#' This function calculates the total and indirect effects for both endogenous
#' and exogenous variables in a Structural Equation Model (SEM).
#'
#' @param BETA A square matrix representing the direct effects between endogenous variables.
#' @param GAMMA A matrix representing the direct effects of exogenous variables on endogenous variables.
#'
#' @return A list containing:
#' \describe{
#'   \item{total_effect}{A matrix combining the total effects of exogenous variables on endogenous variables
#'   and the total effects between endogenous variables.}
#'   \item{indirect_effect}{A matrix combining the indirect effects of exogenous variables on endogenous variables
#'   and the indirect effects between endogenous variables.}
#' }
#'
#' @details
#' The total effects are computed as:
#' \deqn{(I - BETA)^{-1} - I} for endogenous variables, and
#' \deqn{(I - BETA)^{-1} \%*\% GAMMA} for exogenous variables.
#'
#' The indirect effects are derived by subtracting the direct effects (BETA or GAMMA)
#' from the total effects.
#'
#' @examples
#' BETA <- matrix(c(0, 0.5, 0.25, 0), nrow = 2, byrow = TRUE)
#' GAMMA <- matrix(c(-0.3, 0.5, 0, 0, 0, 0, 0.5, 0.25), nrow = 2, byrow = TRUE)
#' result <- compute_effect(BETA, GAMMA)
#' print(result$total_effect)
#' print(result$indirect_effect)
#'
#' @export
compute_effect  <- function(BETA, GAMMA) {

  inv_I_B <- solve(diag(nrow(BETA)) - BETA)
  # total effect endo (endo → endo)
  total_endo <- inv_I_B - diag(nrow(BETA))
  # indirect effects endogenous (endo → endo)
  indirect_endo <- total_endo - BETA

  # total effect exo (exo → endo)
  total_exo <- inv_I_B %*% GAMMA
  # Indirect effect exo (exo → endo)
  indirect_exo <- total_exo - GAMMA


  # [exo → endo | endo → endo]
  total_effect <- cbind(total_exo, total_endo)


  indirect_effect <- cbind(indirect_exo, indirect_endo)

  return(list(total_effect = round(total_effect, 6),
              indirect_effect = round(indirect_effect, 6)))

}







effect_infer <- function(BETA, GAMMA, lengths_parameter, VCOV){

  start <- lengths_parameter[1] + lengths_parameter[2] + 1
  end <- start + lengths_parameter[3] + lengths_parameter[4] -1
  vcov_beta_gamma <- VCOV[start:end, start:end]

  V <- filtering_matrix_effect(BETA)
  inv_I_B <- solve(diag(nrow(BETA)) - BETA)
  K <- kronecker(inv_I_B, t(inv_I_B)) - diag(ncol(BETA)^2)
  J <- t(V) %*% K
  vcov_indirect <- t(J) %*% vcov_beta_gamma %*% J







}



filtering_matrix_effect <- function(M){

  m <- nrow(M)
  n <- ncol(M)
  s <- lengths_parameter[3] + lengths_parameter[4]


  mat_index <- matrix(1:(m * n), m, n)

  # Creation of Vb matrix

  idx_row <- mat_index[M != 0]
  idx_in_theta <- t(mat_index)[t(M) != 0]
  idx_col <- n_gamma + match(idx_row, idx_in_theta)

  # 6. Remplissage de Vg (m*n lignes x s colonnes)
  V <- matrix(0, m * n, s)
  V[cbind(idx_row, idx_col)] <- 1

  return(V)


}








