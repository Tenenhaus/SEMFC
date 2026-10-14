
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
#' @importFrom numDeriv jacobian
#'
#' @keywords internal
information_matrix <- function(x, model){
  JAC <- numDeriv::jacobian(lvm_ml, x = x, model = model, jac = TRUE)
  full_jac <- sapply(1:NCOL(JAC),
                        function(col){
                          ds_dt <- matrix(0, sum(model$block_sizes), sum(model$block_sizes))
                          ds_dt[upper.tri(ds_dt, diag = T)] <- JAC[, col]
                          ds_dt <- ds_dt + t(ds_dt) - diag(diag(ds_dt))
                        }, simplify = FALSE
      )

  nb_param <- length(x)
  Sinv <- solve(lvm_ml(x = x, model = model, jac = F)$sigma_implied)
  Sinv_full_jac <- lapply(full_jac, function(fj) as.matrix(Sinv %*% fj))
  I_ij <- function(i, j) {
    0.5 * sum(diag(Sinv_full_jac[[i]] %*% Sinv_full_jac[[j]]))
  }
  # Generate only indices of the upper triangular part
  index_upper <- which(upper.tri(matrix(0, nrow = nb_param, ncol = nb_param), diag = TRUE), arr.ind = TRUE)

  # Calculate only the elements of the upper triangular part
  I_upper <- mapply(I_ij, index_upper[, "row"], index_upper[, "col"])

  I <- matrix(0, nrow = nb_param, ncol = nb_param)
  I[upper.tri(I, diag = TRUE)] <- I_upper
  I <- I + t(I) - diag(diag(I))

  return(I)
}



#' Compute Jacobian of Constraint Functions
#'
#' Computes the transpose of the Jacobian matrix for equality constraints
#' in formative measurement models.
#'
#' @param x Numeric vector of model parameters.
#' @param S Sample covariance matrix.
#' @param model List containing the model structure, including block sizes,
#'   measurement modes, and other relevant information.
#'
#' @return Matrix H of dimension length(x) × r, where r is the number of formative blocks.
#'
#' @importFrom numDeriv jacobian
#'
#' @keywords internal
Jac_constraints <- function(x, S, model){
  # transpose of the jacobian of constraint function
  H <- t(numDeriv::jacobian(heq1, x = x, S=S, model = model))

  return(H)

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
  H <- matrix(0, t, r)
  if (r>0){
    H  <- Jac_constraints(x, S, model)
  }



  M <- rbind(cbind(I+H%*%t(H), H),
             cbind(t(H), matrix(0, r, r)))

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