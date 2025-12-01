
#' Convert SVD-SEM Parameter  to Vector Form
#'
#' Transforms structural equation model parameters from matrix/list format into
#' a single concatenated parameter vector for optimization or inference. The function
#' handles both reflective and formative measurement models.
#'
#' @param lambda List of loading vectors, one per block. Each element contains the
#'   loadings for the indicators in that block.
#' @param P_EXO Symmetric correlation/covariance matrix of exogenous latent variables.
#' @param G Matrix of gamma coefficients (exogenous → endogenous paths). Rows represent
#'   endogenous variables, columns represent exogenous variables.
#' @param B Matrix of beta coefficients (endogenous → endogenous paths). Rows and columns
#'   represent endogenous variables.
#' @param P_ENDO Symmetric correlation/covariance matrix of endogenous latent variables.
#' @param residual_variance List of residual variances for reflective blocks. Each element
#'   contains the residual variances for the indicators in that block.
#' @param S_composites List of empirical covariance matrices for formative composite blocks.
#'   Each element is a covariance matrix for a formative block.
#' @param mode Character vector indicating the measurement mode for each block:
#'   "reflective" or "formative".
#'
#' @return Unnamed numeric vector containing all free parameters in the following order:
#'   \enumerate{
#'     \item Loadings (all blocks concatenated)
#'     \item Exogenous correlations (upper triangle of P_EXO, excluding diagonal)
#'     \item Gamma coefficients (non-zero elements, row by row)
#'     \item Beta coefficients (non-zero elements, row by row)
#'     \item Endogenous correlations (upper triangle of P_ENDO, excluding diagonal)
#'     \item Residual variances (reflective blocks) or composite covariances
#'         (upper triangle with diagonal, formative blocks)
#'   }
#'
#' @details
#' For formative blocks, the upper triangular part (including diagonal) of the
#' empirical covariance matrix is included. For reflective blocks, only residual
#' variances are included. This ordering must match the structure expected by
#' `length_parameters_svd()` and reconstruction functions.
#'
#' @seealso \code{\link{length_parameters_svd}} for computing the expected vector length.
#'
#' @keywords internal

parameters_svd <- function(lambda,
                           P_EXO,
                           G,
                           B,
                           P_ENDO,
                           residual_variance,
                           S_composites,
                           mode){

  J <- length(mode)

  S_composites_upper <- lapply(S_composites, function(matrix) matrix[upper.tri(matrix, diag = TRUE)])

  # empirical covariance for composite blocks or residual_variance for reflective
  diag_jj <- vector("list", J)
  diag_jj[mode == "formative"] <- S_composites_upper
  diag_jj[mode != "formative"] <- residual_variance


  theta_vect <-
    Reduce("c",
      c(Reduce("c", lambda),
        P_EXO[upper.tri(P_EXO)],
        apply(G, 1, function(row) row[row != 0]),
        apply(B, 1, function(row) row[row != 0]),
        P_ENDO[upper.tri(P_ENDO)],
        Reduce("c", diag_jj)

      )
    )

  return(unname(theta_vect))
}