



information_matrix <- function(x, block_sizes, mode, lengths_parameter, which_exo_endo){
  JAC <- numDeriv::jacobian(lvm_ml, x = x, block_sizes = block_sizes, mode =mode ,
                            lengths_parameter = lengths_parameter, which_exo_endo = which_exo_endo, jac = TRUE)
  full_jac <- sapply(1:NCOL(JAC),
                        function(col){
                          ds_dt <- matrix(0, sum(block_sizes), sum(block_sizes))
                          ds_dt[upper.tri(ds_dt, diag = T)] <- JAC[, col]
                          ds_dt <- ds_dt + t(ds_dt) - diag(diag(ds_dt))
                        }, simplify = FALSE
      )

  nb_param <- length(x)
  Sinv <- solve(lvm_ml(x = x, block_sizes = block_sizes, mode =mode,
                       lengths_parameter = lengths_parameter, which_exo_endo = which_exo_endo, jac = F)$SIGMA_IMPLIED)
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

Jac_constraints <- function(x, S, block_sizes, mode, lengths_parameter, which_exo_endo){
  # transpose of the jacobian of constraint function
  H <- t(numDeriv::jacobian(heq1_bis, x = x, S=S, block_sizes=block_sizes, mode = mode,
                            lengths_parameter = lengths_parameter, which_exo_endo = which_exo_endo))

  return(H)

}


P_ml <- function(x, S, block_sizes, mode, lengths_parameter,which_exo_endo){

  I <- information_matrix(x, block_sizes, mode, lengths_parameter,which_exo_endo)
  H  <- Jac_constraints(x, S, block_sizes, mode, lengths_parameter, which_exo_endo)
  r  <- ncol(H)
  t <- nrow(H)
  M <- rbind(cbind(I+H%*%t(H), H),
             cbind(t(H), matrix(0, r, r)))
  P <- solve(M)[1:t, 1:t]

  return(P)


}

formatting_ml_infer <- function(fit, SD, lengths_parameter){

  lambda <- unlist(unname(fit$lambda))
  gamma <- fit$gamma[fit$gamma!=0]
  beta <- fit$beta[fit$beta!=0]

  start_indices_in_x <- cumsum(c(1, head(lengths_parameter, -1)))
  lambda_start_index <- start_indices_in_x[1]
  lambda_end_index <- lambda_start_index + length(lambda) - 1
  gamma_start_index <- start_indices_in_x[3]
  gamma_end_index <- gamma_start_index + length(gamma) - 1
  beta_start_index <- start_indices_in_x[4]
  beta_end_index <- beta_start_index + length(beta) - 1

  sd_lambda <- SD[lambda_start_index: lambda_end_index]
  sd_gamma <- SD[gamma_start_index: gamma_end_index]
  sd_beta <- SD[beta_start_index: beta_end_index]

  z_lambda <- lambda/sd_lambda
  z_gamma<- gamma/sd_gamma
  z_beta <- beta/sd_beta


  table_lambda <- data.frame(Estimate = lambda,
                             std = sd_lambda,
                             z_score = z_lambda,
                             pval = unlist(lapply(z_lambda, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))
  )

  table_gamma <- data.frame(Estimate = gamma,
                            std = sd_gamma,
                            z_score = z_gamma,
                            pval = unlist((lapply(z_gamma, function (z) 2*pnorm(abs(z), lower.tail = FALSE))))
  )

  table_beta <- data.frame(Estimate = beta,
                           std = sd_beta,
                           z_score = z_beta,
                           pval = unlist(lapply(z_beta, function (z) 2*pnorm(abs(z), lower.tail = FALSE)))
  )



  out <- list(
    lambda = table_lambda,
    gamma = table_gamma,
    beta = table_beta

  )

  return(out)


}



mlSEM_infer <- function(x, S, block_sizes, mode, lengths_parameter, N, fit,which_exo_endo){

  P_ml <-P_ml(x, S, block_sizes, mode, lengths_parameter, which_exo_endo)
  VCOV <- P_ml/N
  SD <- sqrt(diag(VCOV))

  table <- formatting_ml_infer(fit, SD, lengths_parameter)

  out <- list(
    estimate = table,
    VCOV = VCOV
  )

  return(out)


}







z_H0 <- function(x, S, X, C, mode, L){
  N <- nrow(X)
  P <- P_ml(x, S, X, C, mode)
  z <- sqrt(N)*(L %*% x)/sqrt(L%*%P%*%t(L))
  pval <- 2*pnorm(z, lower.tail = F)

  return(pval)

}















