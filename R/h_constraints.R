# restrictions for the minimization of the Loglikelihood function



#' @title Equality Constraints for Loglikelihood Minimization (H1)
#'
#' @description
#' Generalized version for multiple blocks with formative/reflective modes.
#'
#' @param x Numeric vector of parameters.
#' @param S Sample covariance matrix.
#' @param model A list containing model specifications with the following elements:
#'   \describe{
#'     \item{block_sizes}{Integer vector specifying the number of indicators in each block.}
#'     \item{mode}{Character vector indicating the measurement mode for each block
#'       ("formative" or "reflective").}
#'   }
#'
#' @return Vector of equality constraints.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 300
#' # Generate sample data (6 blocks with 3 indicators each)
#' X <- matrix(rnorm(n * 18), ncol = 18)
#' S <- cov(X)
#' # Model specification
#' block_sizes <- c(3, 3, 3, 3, 3, 3)
#' mode <- c("formative", "formative", "formative", "formative",
#'           "reflective", "reflective")
#'
#'
#' model <- list(
#'  block_sizes = block_sizes,
#' mode = mode
#' )
#'
#' x  <- rnorm(61)
#'
#' heq1(x, S, model)
#'}
#' @keywords internal
heq1 <- function (x, S, model){

  block_sizes <- model$block_sizes
  mode <- model$mode

  loadings <- get_loadings(x, block_sizes)



 # list of lengths: diagonal covariances for reflective blocks,
 # or lower-triangular covariance values for formative blocks
  lengths_values_cov <- block_sizes
  lengths_values_cov[mode == "formative"] <- (block_sizes[mode == "formative"]^2 + block_sizes[mode == "formative"]) / 2
  # number of parameters for covariance
  total_cov_parameter <- sum(lengths_values_cov)
    # the coefficient are stocked at the end of x
  initial_start_index_cov <- length(x) - total_cov_parameter +1
  end_endex_cov <- initial_start_index_cov + total_cov_parameter - 1

  # part of the vector corresponding to covariance blocks
  extracted_parameters_cov <- x[initial_start_index_cov:end_endex_cov]
  # list of parameters corresponding to each covariance bloc
  list_cov <- split(extracted_parameters_cov,
                    rep(seq_along(lengths_values_cov), lengths_values_cov))

  # list of formative covariance matrices
  S_composite <- lapply(list_cov[mode == "formative"], build_formative_S_diag)
  h <- as.vector(mapply(function(S_composite_i, loadings_i) {
    t(loadings_i) %*% solve(S_composite_i) %*% loadings_i - 1
  }, S_composite, loadings[mode == "formative"], SIMPLIFY = TRUE))

  return(h)

}