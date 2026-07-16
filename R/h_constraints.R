# restrictions for the minimization of the Loglikelihood function

#' @title Equality Constraints for Loglikelihood Minimization (H0)
#'
#' @description
#' Restrictions for the minimization of the Loglikelihood function.
#'
#' @param x Numeric vector of parameters.
#' @param S Sample covariance matrix.
#'
#' @return Vector of equality constraints.
#'
#' @keywords internal
heq0 <- function(x, S) {
  
  #cov between MVs
  S1 = matrix(c(x[32], x[33], x[35],
                x[33], x[34], x[36],
                x[35], x[36], x[37]), 3, 3)
  #cov between MVs
  S2 = matrix(c(x[38], x[39], x[41],
                x[39], x[40], x[42],
                x[41], x[42], x[43]), 3, 3)
  #cov between MVs
  S3 = matrix(c(x[44], x[45], x[47],
                x[45], x[46], x[48],
                x[47], x[48], x[49]), 3, 3)
  #cov between MVs
  S4 = matrix(c(x[50], x[51], x[53],
                x[51], x[52], x[54],
                x[53], x[54], x[55]), 3, 3)
  
  l1 = x[1:3] ; l2 = x[4:6] ; l3 = x[7:9]
  l4 = x[10:12]
  
  P_EXO = matrix(c(1, x[19], x[20], x[22],
                   x[19], 1, x[21], x[23],
                   x[20], x[21], 1, x[24],
                   x[22], x[23], x[24], 1), 4, 4)
  
  P_ENDO = matrix(c(1, x[31],
                    x[31], 1), 2, 2)
  
  # path coefficients
  G = matrix(c(x[25], x[26], 0, 0,
               0, 0, x[27], x[28]), 2, 4, byrow = TRUE)
  
  # path coefficients
  B = matrix(c(0, x[29],
               x[30], 0), 2, 2, byrow = TRUE)
  
  R = rbind(cbind(P_EXO, P_EXO%*%t(G)%*%t(solve(diag(2) - B))),
            cbind(solve(diag(2) - B)%*%G%*%P_EXO, P_ENDO))
  
  h <- c(rep(0,5))
  h[1] <- t(l1)%*%solve(S1)%*%l1 - 1
  h[2] <- t(l2)%*%solve(S2)%*%l2 - 1
  h[3] <- t(l3)%*%solve(S3)%*%l3 - 1
  h[4] <- t(l4)%*%solve(S4)%*%l4 - 1
  h[5] <- x[27] - x[26]
  
  return(h)
}

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



heq <- function (x, S, model){

  block_sizes <- model$block_sizes
  mode <- model$mode
  r <- sum(mode == "formative")

  loadings <- get_loadings(x, block_sizes)
  loadings_formative <- loadings[mode == "formative"]
  Lambda_formative <- Matrix::bdiag(loadings_formative)



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


  inv_Sigma_formative <- bdiag(lapply(S_composite, solve))

  h <- diag(t(Lambda_formative) %*% inv_Sigma_formative %*% Lambda_formative) - rep(1, r)

  return(h)

}


grad_heq <- function(x, S, model) {
  est <- lvm_ml(x, model, jac = FALSE)
  lambda <- est$lambda
  S_composites <- est$S_composites


  J_H <- compute_gradient_constraint(lambda, S_composites, model$block_sizes, model$lengths_theta, model$mode)

  return(as.matrix(J_H))
}