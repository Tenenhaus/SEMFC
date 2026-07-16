



lvm <- function(R, C){
  gr <- igraph::graph_from_adjacency_matrix(C)
  which_exo_endo <- ind_exo_endo(C)
  rownames(C) <- colnames(C) <- rownames(R)

  H <- which_exo_endo$ind_exo
  J <- which_exo_endo$ind_endo


  if(any(colSums(C[, J, drop = FALSE]) == 0)){
    bad_idx <- which(colSums(C[, J, drop = FALSE]) == 0)
    stop(paste("The structural model is not properly specified for endo variable(s):",
               paste(colnames(C[, J, drop = FALSE])[bad_idx], collapse = ", ")))
  }
  R_proj <- R
  if(!igraph::is_dag(gr)){
    R_proj <-  R[, H, drop = F] %*% solve(R[H, H, drop = F]) %*% t(R[, H, drop = F])
  }
  I <- diag(length(J))
  S <- C[, J, drop = F]
  s <- as.vector(S)
  S <- diag(ncol(C)*length(J))[, s == 1]
  gb_vec <- solve(t(S)%*%kronecker(I, R_proj)%*%S)%*%(t(S)%*%as.vector(R_proj[, J, drop = F]))

  GAMMA_BETA <- matrix(S%*%gb_vec, ncol(C), length(J))
  dimnames(GAMMA_BETA) <- dimnames(C[, J, drop = F])

  BETA <- t(GAMMA_BETA[J, , drop = F])
  GAMMA <- t(GAMMA_BETA[H, , drop = F])


  PI <- solve(diag(NROW(BETA))- BETA)
  if(!igraph::is_dag(gr)){
    PSI <- (diag(NROW(BETA))- BETA)%*%R[J,J]%*%t(diag(NROW(BETA))- BETA) -
      GAMMA%*%R[H,H]%*%t(GAMMA)
  } else {
    D <- diag(diag(ncol(BETA)) - (PI%*%GAMMA%*%R[H, H]%*%t(GAMMA)%*%t(PI)))
    diag_PSI <- drop(solve(PI*PI)%*%D)
    PSI <- if(length(diag_PSI) == 1) {
      matrix(diag_PSI, 1, 1)
    } else {
      diag(diag_PSI)
    }

  }


  dimnames(PSI) <- dimnames(BETA)

  R2 <- 1-diag(PSI)

  R_LVM <- matrix(0, NCOL(C), NCOL(C))

  if(!igraph::is_dag(gr)){
    R_LVM[H, H] <- R[H, H]
    R_LVM[J, J] <- R[J, J]
    R_LVM[H, J] <- R[H, H]%*%t(GAMMA)%*%t(PI)
    R_LVM[J, H] <- PI%*%GAMMA%*%R[H, H]
  }else{
    R_LVM[H, H] <- R[H, H]
    R_LVM[H, J] <- R[H, H]%*%t(GAMMA)%*%t(PI)
    R_LVM[J, H] <- PI%*%GAMMA%*%R[H, H]
    R_LVM[J, J] <- PI%*%(GAMMA%*%R[H, H]%*%t(GAMMA) + PSI)%*%t(PI)
  }

  dimnames(R_LVM) <- dimnames(C)


  return(list(gr = gr,
              BETA = BETA, GAMMA = GAMMA,
              PSI = PSI,
              R2 = R2,
              P_EXO = R_LVM[H, H],
              P_ENDO = R_LVM[J, J],
              R_LVM = R_LVM))



}