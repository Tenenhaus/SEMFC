
#' Convert SVD-SEM Parameters to Vector Form
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
#' @param model A list containing the model structure, including:
#'   \describe{
#'     \item{n_blocks}{Number of blocks in the model.}
#'     \item{mode}{Character vector indicating the measurement mode for each block:
#'       \code{"reflective"} or \code{"formative"}.}
#'     \item{dag}{Logical; if \code{TRUE}, endogenous correlations are excluded
#'       (directed acyclic graph assumption).}
#'   }
#'
#' @return Unnamed numeric vector containing all free parameters in the following order:
#'   \enumerate{
#'     \item Loadings (all blocks concatenated)
#'     \item Exogenous correlations (lower triangle of P_EXO, excluding diagonal)
#'     \item Gamma coefficients (non-zero elements, row by row)
#'     \item Beta coefficients (non-zero elements, row by row)
#'     \item Endogenous correlations (lower triangle of P_ENDO, excluding diagonal)
#'     \item Residual variances (reflective blocks) or composite covariances
#'       (lower triangle with diagonal, formative blocks)
#'   }
#'
#' @details
#' For formative blocks, the lower triangular part (including diagonal) of the
#' empirical covariance matrix is included. For reflective blocks, only residual
#' variances are included.
#'
#'
#' @keywords internal

parameters_svd <- function(lambda,
                           P_EXO,
                           G,
                           B,
                           P_ENDO,
                           residual_variance,
                           S_composites,
                           model){

  J <- model$n_blocks
  mode <- model$mode
  dag <- model$dag

  S_composites_lower <- lapply(S_composites, function(mat) mat[lower.tri(mat, diag = TRUE)])

  # empirical covariance for composite blocks or residual_variance for reflective
  diag_jj <- vector("list", J)
  diag_jj[mode == "formative"] <- S_composites_lower
  diag_jj[mode != "formative"] <- residual_variance


  vect_lambda <- Reduce("c", lambda)
  vect_exo <- P_EXO[lower.tri(P_EXO)]
  vect_endo <- P_ENDO[lower.tri(P_ENDO)]
  if (dag){
      vect_endo <- numeric(0)
  }

  vect_beta <- apply(B, 1, function(row) row[row != 0])
  vect_gamma <- apply(G, 1, function(row) row[row != 0])
  vect_cov <- Reduce("c", diag_jj)



  theta_vect <-
    Reduce("c",
      c(vect_lambda,
        vect_exo,
        vect_gamma,
        vect_beta,
        vect_endo,
        vect_cov
      )
    )



  return(unname(theta_vect))
}