

source('R/SEMFC/sem_f_c.R')




get_lambdas <- function(A, a, mode, bias = FALSE){

  J <- length(A)

  a <- lapply(a, function(x) {if (x[1]>0) {x<-x} else {x<--x}})
  names(a) <- names(A)

  #Compute disattenuation
  d <- rep(0, J)

  for(j in seq_len(J)){
    d[j] <- correction(A[[j]], a[[j]],
                       mode = mode[j],
                       bias = bias)$d
  }

  names(d) = names(A)

  #... and apply the correction
  lambda <- mapply("*", a, d, SIMPLIFY = FALSE)

  return(list(lambda = lambda, d = d))

}


get_R <- function(A, C, a, d, bias = FALSE){

  J <- length(A)

  Phat <- matrix(0, J, J)
  for (j in seq_len(J)){
    for (k in j:J){
          Phat[j, k] <- t(a[[j]])%*%cov2(A[[j]], A[[k]], bias = bias)%*%a[[k]]/(d[j]*d[k])
    }
  }

  Phat <- Phat + t(Phat)
  diag(Phat) <- 1
  colnames(Phat) <- rownames(Phat) <- names(A)

  #Call lvm() for the structural model
  lv <- lvm(Phat, C)

  return(lv$R_LVM)


}





create_us <- function(R, lambdas) {


  # u_1 <- c(R[1,2]*lambdas[[2]], R[1,3]*lambdas[[3]], R[1,4]*lambdas[[4]], R[1,5]*lambdas[[5]], R[1,6]*lambdas[[6]])
  # u_2 <- c(R[2,1]*lambdas[[1]], R[2,3]*lambdas[[3]], R[2,4]*lambdas[[4]], R[2,5]*lambdas[[5]], R[2,6]*lambdas[[6]])
  # u_3 <- c(R[3,1]*lambdas[[1]], R[3,2]*lambdas[[2]], R[3,4]*lambdas[[4]], R[3,5]*lambdas[[5]], R[3,6]*lambdas[[6]])
  # u_4 <- c(R[4,1]*lambdas[[1]], R[4,2]*lambdas[[2]], R[4,3]*lambdas[[3]], R[4,5]*lambdas[[5]], R[4,6]*lambdas[[6]])
  # u_5 <- c(R[5,1]*lambdas[[1]], R[5,2]*lambdas[[2]], R[5,3]*lambdas[[3]], R[5,4]*lambdas[[4]], R[5,6]*lambdas[[6]])
  # u_6 <- c(R[6,1]*lambdas[[1]], R[6,2]*lambdas[[2]], R[6,3]*lambdas[[3]], R[6,4]*lambdas[[4]], R[6,5]*lambdas[[5]])
  #
  #
  # us <- list(u_1, u_2, u_3, u_4, u_5, u_6)
  #

  us <- vector("list", 6)
  for (j in 1:6) {
    u_j <- c()
    for (i in setdiff(1:6, j)) {
      u_j <- c(u_j, R[j, i] * lambdas[[i]])
    }
    us[[j]] <- u_j
  }
  return(us)
}



soft = function(x, d) {
  return(sign(x) * pmax(0, abs(x) - d))
}


l2n = function(vec) {
  a <- sqrt(sum(vec ^ 2))
  if (a == 0) a <- 0.05
  return(a)
}




BinarySearch = function(argu, sumabs) {
  if (l2n(argu) == 0 || sum(abs(argu / l2n(argu))) <= sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu)) - 1e-05
  iter <- 1
  while (iter < 150) {
    su <- soft(argu, (lam1 + lam2) / 2)
    if (sum(abs(su / l2n(su))) < sumabs) {
      lam2 <- (lam1 + lam2) / 2
    } else {
      lam1 <- (lam1 + lam2) / 2
    }
    if ((lam2 - lam1) < 1e-06) return((lam1 + lam2) / 2)
    iter <- iter + 1
  }
  warning("Didn\'t quite converge")
  return((lam1 + lam2) / 2)
}







newpmd <- function (Y, C, mode, sumabsv, niter, tol){

  model<- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen=FALSE,
                       len_seq = NULL,
                       nfold=NULL,
                       niter=NULL,
                       sparse_val=NULL)
  model$fit_svd()

  lambdas <- model$svd_parameters$lambda
  R <- model$svd_parameters$P_IMPLIED



  us <- create_us(R, lambdas)
  vs <- vector("list", 6)
  for (x in 1:length(Y)){
    u <- us[[x]]
    data <- t(t(Y[[x]])%*%Reduce("cbind", Y[-x]))
    argv <- t(u) %*% data
    lamv <- BinarySearch(argv, sumabsv[[x]])
    sv <- soft(argv, lamv)
    v <- matrix(sv / l2n(sv), ncol = 1)
    vs[[x]] <- v
  }


  iter <- 0
  diff <- tol + 1

  while (iter < niter & diff > tol){
    v_old <- vs

    fit <- get_lambdas(Y, vs, mode, bias = FALSE)
    lambdas_k <- fit$lambda
    d_k <- fit$d
    R_k <- get_R(Y, C, vs, d_k, bias = FALSE)
    us <- create_us(R_k, lambdas_k)

    for (x in 1:length(Y)){
      u <- us[[x]]
      data <- t(t(Y[[x]])%*%Reduce("cbind", Y[-x]))
      argv <- t(u) %*% data
      lamv <- BinarySearch(argv, sumabsv[[x]])
      sv <- soft(argv, lamv)
      v <- matrix(sv / l2n(sv), ncol = 1)
      vs[[x]] <- v
    }

    # Compute difference for convergence check
    diff <- max(sapply(seq_along(vs), function(i) max(abs(vs[[i]] - v_old[[i]]))))
    iter <- iter + 1
  }

  return(vs)

}