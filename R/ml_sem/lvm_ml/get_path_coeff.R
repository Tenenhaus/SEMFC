

#' Construct Path Coefficient Matrix from Parameter Vector
#'
#' This function generates a path coefficient matrix for endogenous and exogenous variables
#' based on the parameter vector `x` and the provided relationships between variables.
#'
#' @param x A numeric vector of parameters. Contains path coefficient values to be extracted.
#' @param list_linked_exo_endo A list where each element represents the linked exogenous or
#'   endogenous variables for each endogenous variable. Keys correspond to variable names.
#' @param exo_or_endo_variable A named numeric vector representing the exogenous or endogenous
#'   variables. Used to define the structure of the path coefficient matrix.
#' @param initial_start_index An integer specifying the index in `x` where the path coefficient
#'   parameters start.
#'
#' @return A matrix where rows correspond to endogenous variables and columns to exogenous or
#'   endogenous variables, with path coefficients filled in where applicable.
#'
#' @examples
#' x <- c(0.8, 0.5, 0.3, 0.7, 0.8, 0.4)
#' list_linked_exo_endo <- list(c(X1 = 1, X2 = 2), c(X3 = 3, X4 = 4))
#' exo_or_endo_variable <- c(X1 = 1, X2 = 2, X3 = 3, X4 = 4)
#' initial_start_index <- 1
#' get_path_coeff(x, list_linked_exo_endo, exo_or_endo_variable, initial_start_index)

get_path_coeff <- function(x, list_linked_exo_endo, exo_or_endo_variable, initial_start_index) {


  # list of number of linked exos/endos for each endo variable
  # length_exo_endo <- sapply(list_linked_exo_endo, length)
  length_exo_endo <- sapply(list_linked_exo_endo,
                            function(sublist) {
                              ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                                     return(0),
                                     return(length(sublist)))
                            })


  # Initialisation of index for the rows of the path coeff matrix
  start_index <- initial_start_index + c(0, cumsum(head(length_exo_endo, -1)))
  end_index <- start_index + length_exo_endo - 1

  # in the matrix, replace non zero value by x values
  Matrix_row <- mapply(
    function(linked_variable_i, start_ind, end_ind) {
      # Create a vector based on list_linked_exo_endo[i]
      result <- ifelse(names(exo_or_endo_variable) %in% names(linked_variable_i),
                       exo_or_endo_variable, 0)
      result[result != 0] <- x[start_ind:end_ind]
      names(result) <- names(exo_or_endo_variable)
      return(result)
    },
    list_linked_exo_endo,                  # List of linked variables
    start_index,
    end_index,
    SIMPLIFY = FALSE
  )

  Matrix_path <- do.call(rbind, Matrix_row)

  return(Matrix_path)
}
