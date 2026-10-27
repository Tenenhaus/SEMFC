
#' Compute Information Matrix for ML Estimation
#'
#' Computes the Fisher information matrix using numerical derivatives of the
#' implied covariance matrix with respect to model parameters.
#'
#' @param x Numeric vector of model parameters.
#' @param model List containing the model structure, including block sizes and
#'   other relevant information for the computation.
#'
#' @return Symmetric information matrix of dimension length(x) × length(x).
#'
#' @details
#' The information matrix is computed as:
#' \deqn{I_{ij} = \frac{1}{2} \text{tr}(\Sigma^{-1} \frac{\partial\Sigma}{\partial\theta_i} \Sigma^{-1} \frac{\partial\Sigma}{\partial\theta_j})}
#' Only the upper triangular part is computed for efficiency.
#'
#'
#' @keywords internal
information_matrix <- function(x, model){
  est <- lvm_ml(x, model, jac = FALSE)
  J_Sigma_theta <- compute_J_Sigma_theta(model$block_sizes, model$lengths_theta, model$mode,model$dag,
                                         est$p_implied, model$relation_matrix, est$lambda, est$beta, est$gamma)
  H <- compute_Hessian(est$sigma_implied, J_Sigma_theta)
  I <- 0.5 * H
  return(I)
}



' Compute Projection Matrix for Constrained ML Estimation
#'
#' Computes the projection matrix P for variance-covariance estimation under
#' equality constraints in structural equation models with formative blocks.
#'
#' @param x Numeric vector of model parameters.
#' @param S Sample covariance matrix.
#' @param model List containing the model structure, including:
#'   \itemize{
#'     \item \code{mode}: Character vector indicating the measurement mode for each block.
#'     \item \code{block_sizes}: Integer vector specifying the number of indicators in each block.
#'   }
#'
#' @return Projection matrix P of dimension length(x) × length(x).
#'
#' @details
#' For constrained optimization, the projection matrix is computed as:
#' \deqn{P = [M]_{1:t,1:t}^{-1}}
#' where M is the bordered information matrix:
#' \deqn{M = \begin{bmatrix} I + HH' & H \\ H' & 0 \end{bmatrix}}
#'
#'
#' @importFrom MASS ginv
#'
#' @keywords internal
P_ml <- function(x, S, model){

  mode <- model$mode

  I <- information_matrix(x, model)
  t <- nrow(I)
  r <- length(mode[mode=='formative'])
  H <- matrix(0, r, t)
  if (r>0){
    H  <- grad_heq(x, S, model)
  }



  M <- rbind(cbind(I+t(H)%*%H, t(H)),
             cbind(H, matrix(0, r, r)))

  invM <- tryCatch(
    solve(M),
    error = function(e) {
      warning("Inversion of M failed; using MASS::ginv().")
      MASS::ginv(M)
    }
  )

  P <- invM[1:t, 1:t]

  return(P)


}