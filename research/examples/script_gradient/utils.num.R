information_matrix_old <- function(x, model){
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


Jac_constraints <- function(x, S, model){
  # transpose of the jacobian of constraint function
  H <- t(numDeriv::jacobian(heq1, x = x, S=S, model = model))

  return(H)

}

P_ml_old <- function(x, S, model){

  mode <- model$mode

  I <- information_matrix_old(x, model)
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


f <- function (x) { return(F1(x, S, model)) }
grad <- function(x){ return(grad_F1(x, S, model)) }
h_eq <- function(x){ return(heq(x, S, model)) }
grad_h_eq <- function(x){ return(grad_heq(x, S, model)) }
hess <- function(x){ return(hess_F1(x, S, model)) }