



one_step_estimator <- function(init_estimate, model, data){
  theta_init <- init_estimate$theta
  lambda <- init_estimate$lambda
  S_diag_composites <- data$S_diag_composites
  block_sizes <- model$block_sizes
  lengths_theta <- model$lengths_theta
  mode <- model$mode
  S <- data$cov_S
  BETA <- init_estimate$beta
  GAMMA <- init_estimate$gamma
  sigma_implied <- init_estimate$sigma_implied
  r <- sum(mode == "formative")


  jac_Sigma_theta <- compute_J_Sigma_theta(block_sizes,
                                           lengths_theta,
                                           mode, model$dag,
                                           init_estimate$p_implied,
                                           model$relation_matrix,
                                           lambda, BETA, GAMMA)

  H <- compute_Hessian(sigma_implied, jac_Sigma_theta)
  grad_F <- as.vector(compute_gradient(sigma_implied, S, jac_Sigma_theta))

  if (r == 0) {
    theta_os <- theta_init - solve(H, grad_F)
    return(theta_os)
  }

  h_eq <- heq(theta_init, S, model)
  J_h_eq <- compute_gradient_constraint(lambda,
                                        S_diag_composites,
                                        block_sizes,
                                        lengths_theta,
                                        mode)
  ktt_mat <- rbind(cbind(H, t(J_h_eq)),
                   cbind(J_h_eq, matrix(0, r, r)))
  ktt_vec <- c(grad_F, h_eq)

  theta_os <- theta_init - solve(ktt_mat, ktt_vec)[seq_along(theta_init)]

  loadings_normalised_formative <- normalization_os(theta_os, model)
  loadings_os <- get_loadings(theta_os, block_sizes)
  loadings_os[mode == "formative"] <- loadings_normalised_formative

  theta_os[seq_along(unlist(lambda))] <- unlist(loadings_os)

  return(theta_os)

}


normalization_os <- function(theta_os, model) {

  block_sizes <- model$block_sizes
  mode <- model$mode
  lengths_theta <- model$lengths_theta
  loadings_os <- get_loadings(theta_os, block_sizes)[mode == "formative"]
  BDIAG <- get_bdiag(theta_os,
                     mode = mode,
                     block_sizes = block_sizes,
                     initial_start_index_cov = cumsum(c(1, head(lengths_theta, -1)))[6])
  S_composites_os <- BDIAG[mode == "formative"]
  loadings_final<- Map(function(x, M) {
    norm_factor <- sqrt(as.numeric(t(x) %*% solve(M, x)))
    return(x / norm_factor)

  }, loadings_os, S_composites_os)

  return(loadings_final)


}