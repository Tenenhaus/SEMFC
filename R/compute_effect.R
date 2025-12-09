



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
  # indirect effects endogènes (endo → endo)
  indirect_endo <- total_endo - BETA

  # total effect exo (exo → endo)
  total_exo <- inv_I_B %*% GAMMA
  # Indirect effect exo (exo → endo)
  indirect_exo <- total_exo - GAMMA


  # [exo → endo | endo → endo]
  total_effect <- cbind(total_exo, total_endo)


  indirect_effect <- cbind(indirect_exo, indirect_endo)

  return(list(total_effect = total_effect,
              indirect_effect = indirect_effect))

}