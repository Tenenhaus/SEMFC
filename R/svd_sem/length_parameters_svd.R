

# compute the length of theta model parameter


#' Compute Length of SVD-SEM Parameter Vector
#'
#' Calculates the total number of free parameters in a structural equation model
#' estimated via SVD decomposition, accounting for loadings, path coefficients,
#' correlations, residual variances, and covariances.
#'
#' @param block_sizes Integer vector specifying the number of indicators in each block.
#' @param which_exo_endo List containing indices and structure information for
#'   exogenous and endogenous latent variables, including:
#'   \itemize{
#'     \item `ind_exo`: Indices of exogenous latent variables
#'     \item `ind_endo`: Indices of endogenous latent variables
#'     \item `Hi`: List of predecessor exogenous variables for each endogenous variable
#'     \item `Ji`: List of predecessor endogenous variables for each endogenous variable
#'   }
#'
#' @return Integer representing the total number of free parameters in the model.
#'
#' @details
#' The total parameter count is computed as the sum of:
#' \itemize{
#'   \item **Loadings**: Sum of all block sizes
#'   \item **Gamma coefficients**: Number of exogenous → endogenous paths
#'   \item **Beta coefficients**: Number of endogenous → endogenous paths
#'   \item **Exogenous correlations**: \eqn{p(p-1)/2} where p = number of exogenous variables
#'   \item **Endogenous correlations**: \eqn{q(q-1)/2} where q = number of endogenous variables
#'   \item **Residual variances**: Number of indicators in reflective blocks
#'   \item **Composite covariances**: \eqn{\sum_{formative} n_i(n_i+1)/2}
#' }
#'
#' @note This function assumes the `mode` variable is available in the calling environment.
#'
#' @keywords internal
length_parameters_svd <- function (block_sizes, which_exo_endo){
  n <- which_exo_endo$ind_exo
  m <- which_exo_endo$ind_endo

  length_lambda <- sum(block_sizes)
  length_gamma <- sapply(which_exo_endo$Hi,
                       function(sublist) {
                          ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                                 return(0),
                                 return(length(sublist)))
                       })
  length_beta <- sapply(which_exo_endo$Ji,
                     function(sublist) {
                        ifelse((length(sublist) == 1 && sublist[[1]] == 0),
                               return(0),
                               return(length(sublist)))
                     })

  dim_exo <- length(n)
  length_exo <- dim_exo * (dim_exo - 1) / 2

  dim_endo <- length(m)
  length_endo <- dim_endo * (dim_endo - 1) / 2

  length_residual_variance <- sum(block_sizes[mode == 'reflective'])

  length_cov_composite <- sum(unlist(lapply(block_sizes[mode == 'formative'], function(n) n * (n + 1) / 2)))

  len_parameters <- length_lambda +
    length_gamma + length_beta +
    length_exo + length_endo +
    length_residual_variance + length_cov_composite

  return(len_parameters)

}
