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
#' \dontrun{
#' BETA <- matrix(c(0, 0.5, 0.25, 0), nrow = 2, byrow = TRUE)
#' GAMMA <- matrix(c(-0.3, 0.5, 0, 0, 0, 0, 0.5, 0.25), nrow = 2, byrow = TRUE)
#' result <- compute_effect(BETA, GAMMA)
#' print(result$total_effect)
#' print(result$indirect_effect)
#' }
#'
#' @keywords internal
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

#' Create a Filtering Matrix for Effects inference
#'
#' This function generates a filtering matrix \( V \) used to map the indices of non-zero elements
#' in a matrix \( M \) to a parameter vector. The filtering matrix is used in the computation
#' of partial derivatives of effects in structural equation modeling.
#'
#' @param M A matrix whose non-zero elements are to be mapped.
#' @param s An integer representing the number of beta and gamma parameters in the model.
#' @param start_index An integer representing the starting index for mapping the parameters.
#'
#' @return A matrix \eqn{\mathbf{V}_{M}} of dimensions \eqn{(m \times n) \times s}, where \eqn{m} and \eqn{n}
#' are the dimensions of \eqn{M}. This matrix represents:
#' \deqn{\mathbf{V}_{M} = \left[ \mathrm{vec} \frac{\partial \mathbf{M}}{\partial \theta_1},
#' \mathrm{vec} \frac{\partial \mathbf{M}}{\partial \theta_2}, \dots,
#' \mathrm{vec} \frac{\partial \mathbf{M}}{\partial \theta_s} \right]}
#' The matrix contains 1s at positions corresponding to the mapping of non-zero elements in \eqn{M}
#' to the parameter vector, and 0s elsewhere.
#'
#' @details
#' The function works as follows:
#' - It creates an index matrix for \( M \) to identify the positions of non-zero elements.
#' - It constructs the filtering matrix \( V \) with 1s at the mapped positions.
#'
#'
#' @keywords internal

filtering_matrix_effect <- function(M, s, start_index) {
  m <- nrow(M)
  n <- ncol(M)
  mat_index <- matrix(1:(m * n), m, n)
  # Creation of V matrix
  idx_row <- mat_index[M != 0]
  idx_col <- start_index + seq_along(idx_row)

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
#' @keywords internal
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

#' Compute Asymptotic Variance-Covariance Matrices for Effects Inference using delta Method
#'
#' This function calculates the variance-covariance matrices for total and indirect effects
#' (both endogenous and exogenous) in a Structural Equation Model (SEM). It uses the partial
#' derivatives of the effects with respect to the model parameters and the variance-covariance
#' matrix of the parameters using delta Method
#'
#' @param BETA A square matrix representing the direct effects between endogenous variables.
#' @param GAMMA A matrix representing the direct effects of exogenous variables on endogenous variables.
#' @param lengths_parameter A numeric vector specifying the lengths of different parameter groups
#' in the model. The third and fourth elements correspond to the number of parameters in GAMMA
#' and BETA, respectively.
#' @param VCOV A variance-covariance matrix of the model parameters.
#'
#' @return A list containing:
#' \describe{
#'   \item{vcov_endo_total}{Variance-covariance matrix for total effects between endogenous variables.}
#'   \item{vcov_endo_indirect}{Variance-covariance matrix for indirect effects between endogenous variables.}
#'   \item{vcov_exo_total}{Variance-covariance matrix for total effects of exogenous variables on endogenous variables.}
#'   \item{vcov_exo_indirect}{Variance-covariance matrix for indirect effects of exogenous variables on endogenous variables.}
#' }
#'
#' @details
#' The function extracts the relevant portion of the variance-covariance matrix for the parameters
#' associated with BETA and GAMMA. It then computes the variance-covariance matrices for the effects
#' using the partial derivatives of the effects and the extracted variance-covariance matrix.
#'
#'
#' @keywords internal
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


