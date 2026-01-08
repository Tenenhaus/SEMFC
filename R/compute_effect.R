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



filtering_matrix_effect <- function(M, s, start_index) {
  m <- nrow(M)
  n <- ncol(M)
  mat_index <- matrix(1:(m * n), m, n)
  # Creation of V matrix
  idx_row <- mat_index[M != 0]
  idx_in_theta <- t(mat_index)[t(M) != 0]
  # the ordering of idx_in_theta is by row, we need to match the ordering of idx_row
  idx_col <- start_index + match(idx_row, idx_in_theta)
  V <- matrix(0, m * n, s)
  V[cbind(idx_row, idx_col)] <- 1
  return(V)
}





#' Compute Partial Derivative of Endogenous Effects with Respect to Parameters
#'
#' This function calculates the partial derivative of endogenous effects (either total or indirect)
#' with respect to the parameter vector theta. It implements the formula from Bollen (1989),
#' Appendix 8A: Asymptotic Variances of Effects:
#' \deqn{\frac{\partial \text{vec}(\mathbf{I}_{\eta\eta})}{\partial \boldsymbol{\theta}} =
#' \mathbf{V}_B' \left( (\mathbf{I} - \mathbf{B})^{-1} \otimes [(\mathbf{I} - \mathbf{B})^{-1}]' - \mathbf{M} \right)}
#' where \eqn{\mathbf{M} = \mathbf{I}_{m^2}} for indirect effects and \eqn{\mathbf{M} = \mathbf{0}_{m^2 \times m^2}} for total effects.
#'
#' @param BETA A square matrix representing the direct effects between endogenous variables.
#' @param s An integer representing the total number of parameters in BETA and GAMMA combined.
#' @param len_vect_gamma An integer representing the number of non-zero elements in GAMMA.
#' @param effect_type A character string specifying the type of effect: "total" (default) or "indirect".
#'
#' @return A matrix \eqn{\mathbf{J}} containing the partial derivatives of the endogenous effects
#' with respect to the parameters.
#'
#' @details
#' The function uses the Kronecker product to compute the derivative matrix according to the
#' matrix calculus formula. For total effects, the identity matrix term is removed (\eqn{\mathbf{M} = 0}).
#' For indirect effects, the identity matrix is retained (\eqn{\mathbf{M} = \mathbf{I}_{m^2}}).
#'
#' @references
#' Bollen, K. A. (1989). Structural Equations with Latent Variables. Wiley.
#' Appendix 8A: Asymptotic Variances of Effects.
#'
#' @export
partial_derivative_endo_effect <- function(BETA, s, len_vect_gamma, effect_type = "total") {
  Vb <- filtering_matrix_effect(BETA, s, len_vect_gamma)
  m <- nrow(BETA)
  inv_I_B <- solve(diag(m) - BETA)
  M_effect_type <- diag(m^2)
  if (effect_type == "total") {
    M_effect_type <- matrix(0,m^2, m^2)
  }
  K <- kronecker(inv_I_B, t(inv_I_B)) - M_effect_type
  J <- t(Vb) %*% K

  return(J)

}


#' Compute Partial Derivative of Exogenous Effects with Respect to Parameters
#'
#' This function calculates the partial derivative of exogenous effects (either total or indirect)
#' with respect to the parameter vector theta. It implements the formula from Bollen (1989),
#' Appendix 8A: Asymptotic Variances of Effects:
#' \deqn{\frac{\partial \text{vec}(\mathbf{T}_{\xi\eta})}{\partial \boldsymbol{\theta}} =
#' \mathbf{V}_B' \left[ (\mathbf{I} - \mathbf{B})^{-1}\mathbf{\Gamma} \otimes [(\mathbf{I} - \mathbf{B})^{-1}]' \right] +
#' \mathbf{V}_\Gamma' \left[ \mathbf{I}_n \otimes [(\mathbf{I} - \mathbf{B})^{-1} - \mathbf{M}]' \right]}
#' where \eqn{\mathbf{M} = \mathbf{I}_m} for indirect effects and \eqn{\mathbf{M} = \mathbf{0}_{m \times m}} for total effects.
#'
#' @param BETA A square matrix representing the direct effects between endogenous variables.
#' @param GAMMA A matrix representing the direct effects of exogenous variables on endogenous variables.
#' @param s An integer representing the total number of parameters in BETA and GAMMA combined.
#' @param effect_type A character string specifying the type of effect: "total" (default) or "indirect".
#'
#' @return A matrix \eqn{\mathbf{J}} containing the partial derivatives of the exogenous effects
#' with respect to the parameters.
#'
#' @details
#' The function uses the Kronecker product to compute the derivative matrix according to the
#' matrix calculus formula. For total effects, the identity matrix term is removed (\eqn{\mathbf{M} = 0}).
#' For indirect effects, the identity matrix is retained (\eqn{\mathbf{M} = \mathbf{I}_m}).
#'
#' @references
#' Bollen, K. A. (1989). Structural Equations with Latent Variables. Wiley.
#' Appendix 8A: Asymptotic Variances of Effects.
#'
#' @keywords internal
partial_derivative_exo_effect <- function(BETA, GAMMA, s, effect_type = "total") {
  m <- nrow(GAMMA)  # number endo
  n <- ncol(GAMMA) # number exo
  len_vect_gamma <- length(GAMMA[ GAMMA != 0])


  # Beta part
  Vb <- filtering_matrix_effect(BETA, s, len_vect_gamma)
  inv_I_B <- solve(diag(m) - BETA)
  Kb <- kronecker(inv_I_B%*%GAMMA, t(inv_I_B))

  # Gamma part
  Vg <- filtering_matrix_effect(GAMMA, s, 0)
  M_effect_type <- diag(m)
  if (effect_type == "total") {
    M_effect_type <- matrix(0,m, m)
  }
  Kg <- kronecker(diag(n), t(inv_I_B - M_effect_type))

  J <- t(Vb) %*% Kb  + t(Vg) %*% Kg

  return(J)

}


effect_infer <- function(BETA, GAMMA, lengths_parameter, VCOV){

  start <- lengths_parameter[1] + lengths_parameter[2] + 1
  end <- start + lengths_parameter[3] + lengths_parameter[4] -1
  vcov_beta_gamma <- VCOV[start:end, start:end]

  s <- lengths_parameter[3] + lengths_parameter[4]

  len_vect_gamma <- lengths_parameter[3]

  J_endo_total <- partial_derivative_endo_effect(BETA, s, len_vect_gamma, effect_type = "total")
  vcov_endo_total <- t(J_endo_total) %*% vcov_beta_gamma %*% J_endo_total

  J_endo_indirect <- partial_derivative_endo_effect(BETA, s, len_vect_gamma, effect_type = "indirect")
  vcov_endo_indirect <- t(J_endo_indirect) %*% vcov_beta_gamma %*% J_endo_indirect

  J_exo_total <- partial_derivative_exo_effect(BETA, GAMMA, s, effect_type = "total")
  vcov_exo_total <- t(J_exo_total) %*% vcov_beta_gamma %*% J_exo_total

  J_exo_indirect <- partial_derivative_exo_effect(BETA, GAMMA, s, effect_type = "indirect")
  vcov_exo_indirect <- t(J_exo_indirect) %*% vcov_beta_gamma %*% J_exo_indirect

  return(list(vcov_endo_total = vcov_endo_total,
              vcov_endo_indirect = vcov_endo_indirect,
              vcov_exo_total = vcov_exo_total,
              vcov_exo_indirect = vcov_exo_indirect))

}


