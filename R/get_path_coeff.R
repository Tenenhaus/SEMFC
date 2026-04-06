

#' Construct Path Coefficient Matrix from Parameter Vector
#'
#' This function generates a path coefficient matrix for endogenous and exogenous variables
#' based on the parameter vector `x` and the provided relationships between variables.
#'
#' @param x A numeric vector of parameters. Contains path coefficient values to be extracted.
#' @param row_exo_endo A named numeric vector representing the row indices (endogenous variables).
#' @param col_exo_endo A named numeric vector representing the column indices (exogenous or endogenous variables).
#' @param C A matrix defining the relationships between variables (1 indicates a relationship exists).
#' @param initial_start_index An integer specifying the index in `x` where the path coefficient
#'   parameters start.
#'
#' @return A matrix where rows correspond to endogenous variables and columns to exogenous or
#'   endogenous variables, with path coefficients filled in where applicable.
#'
#' @examples
#' \dontrun{
#' x <- c(0.8, 0.5, 0.3, 0.7, 0.8, 0.4)
#' row_exo_endo <- c(X5 = 5, X6 = 6)
#' col_exo_endo <- c(X1 = 1, X2 = 2, X3 = 3, X4 = 4)
#' initial_start_index <- 1
#' C <- matrix(c(0, 0, 0, 0, 1, 0,
#'             0, 0, 0, 0, 1, 0,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 0, 1,
#'             0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)
#' get_path_coeff(x, row_exo_endo, col_exo_endo, C, initial_start_index)}
#' @keywords internal

get_path_coeff <- function(x, row_exo_endo, col_exo_endo, C, initial_start_index) {


  vec_m <- as.vector(t(C[col_exo_endo, row_exo_endo, drop = FALSE]))
  M_select <- diag(length(vec_m))[, vec_m == 1, drop = FALSE]
  len_param <- sum(vec_m)
  end_index <- initial_start_index + len_param - 1
  vec_regression <- M_select%*%x[initial_start_index:end_index]

  Matrix_path <- matrix(vec_regression, nrow = length(row_exo_endo), ncol = length(col_exo_endo))

  return(Matrix_path)
}
