

#' Calculate parameter vector lengths for SEM model
#'
#' @description
#' Computes the number of parameters for each component of a structural equation model
#' with factors and composites.
#'
#' @param which_exo_endo List containing indices of exogenous/endogenous variables and
#'   their relationships (elements: ind_exo, ind_endo, Hi, Ji)
#' @param block_sizes Numeric vector of block sizes
#' @param mode Character vector specifying "formative" or "reflective" for each block
#' @param dag Logical indicating whether the structural model is recursive (FALSE) or non-recursive (TRUE)
#' @return A numeric vector containing the number of parameters for:
#'   \enumerate{
#'     \item Factor loadings
#'     \item Upper triangular values for exogenous variable covariances
#'     \item Non-zero elements in the gamma matrix
#'     \item Non-zero elements in the beta matrix
#'     \item Upper triangular values for endogenous variable covariances
#'     \item Measurement error covariances/variances (depends on mode)
#'   }
#'
#' @details For formative models, the number of error covariances includes both
#' variances and covariances: (p^2 + p)/2 for each block of size p.
#' For reflective models, only variances are counted: p for each block of size p.
#'
#' @examples
#' \dontrun{
#' which_exo_endo <- list(
#'   Hi = list(
#'     c(LV1 = 1, LV2 = 2),
#'     c(LV3 = 3, LV4 = 4)
#'   ),
#'   Ji = list(
#'     c(LV6 = 6),
#'     c(LV5 = 5)
#'   ),
#'   ind_endo = c(LV5 = 5, LV6 = 6),
#'   ind_exo = c(LV1 = 1, LV2 = 2, LV3 = 3, LV4 = 4)
#' )
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' mode <- c("formative", "formative", "formative", "formative",
#'           "reflective", "reflective")
#' dag <- FALSE
#' lengths <- get_lengths_theta(which_exo_endo, block_sizes, mode, dag)
#' }
#' @keywords internal
get_lengths_theta <- function(which_exo_endo, block_sizes, mode, dag){

  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo

  ##########################################################
  ####### number of parameter for each part ################
  ##########################################################
  number_loadings <- sum(block_sizes)
  number_upper_values_exo <- length(n) * (length(n) - 1) / 2
  # number_non_zero_G <- sum(lengths(which_exo_endo$Hi))
  number_non_zero_G <- sum(length(unlist(which_exo_endo$Hi)[unlist(which_exo_endo$Hi) != 0]))
  # number_non_zero_B <- sum(lengths(which_exo_endo$Ji))
  number_non_zero_B <- sum(length(unlist(which_exo_endo$Ji)[unlist(which_exo_endo$Ji) != 0]))
  number_upper_values_endo <- length(m) * (length(m) - 1) / 2
  if (dag){
    number_upper_values_endo <- 0
  }
  number_cov <- sum(ifelse(mode == "formative", (block_sizes^2 + block_sizes) / 2, block_sizes))

  lengths_parameter <- c(number_loadings,
                         number_upper_values_exo,
                         number_non_zero_G,
                         number_non_zero_B,
                         number_upper_values_endo,
                         number_cov)

  return(lengths_parameter)

}