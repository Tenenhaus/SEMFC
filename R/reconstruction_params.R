
#' Reconstruct Parameters with Empirical Covariance Composite blocks
#'
#' This function reconstructs a complete vector of parameters used for ML
#' by combining residual variances and upper empirical covariance composite.
#'
#' @param params_ml A vector containing the parameters estimated by maximum likelihood.
#' @param mode A vector indicating the mode of each block.
#' @param S_composites A list of empirical covariance composite matrices.
#' @param lengths_cov_parameter A vector containing the lengths of covariance parameters.
#'
#' @return A vector containing the reconstructed parameters used for build the implied covariance matrix.
#'
#' @examples
#' \dontrun{
#' params_ml <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' mode <- c("formative", "reflective", "reflective")
#' S_composites <- list(matrix(1:4, 2, 2), matrix(5:8, 2, 2))
#' lengths_cov_parameter <- c(2, 3, 1)
#' reconstruction_params(params_ml, mode, S_composites, lengths_cov_parameter)
#' }
#' @keywords internal


reconstruction_params <- function(params_ml, mode, S_composites, lengths_cov_parameter, lengths_parameter){

  J <- length(mode)
  S_composites_upper <- lapply(S_composites, function(matrix) matrix[upper.tri(matrix, diag = TRUE)])
  lengths_residual_variance <- lengths_cov_parameter[mode != "formative"]

  len_cov_part <- lengths_parameter[length(lengths_parameter)] - length(unlist(S_composites_upper))
  vect_cov <- tail(params_ml, len_cov_part)

  residual_variance <- split(vect_cov, rep(1:length(lengths_residual_variance), lengths_residual_variance))


  diag_jj <- vector("list", J)
  diag_jj[mode == "formative"] <- S_composites_upper
  diag_jj[mode != "formative"] <- residual_variance

  full_vect_cov <- Reduce("c", diag_jj)

  return(full_vect_cov)


}