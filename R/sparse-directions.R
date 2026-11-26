# Directions parcimonieuses de SVD-SEM : SPC bloc par bloc.
# c[j] est une borne L1 absolue dans [1, sqrt(p_j)].

fit_spc_directions <- function(blocks, sparsity) {
  X <- lapply(blocks, function(x) scale(as.matrix(x), center = TRUE, scale = FALSE))
  n <- nrow(X[[1]])
  if (is.null(names(blocks))) names(blocks) <- paste0("LV", seq_along(blocks))

  a <- lapply(seq_along(X), function(j) {
    if (ncol(X[[j]]) == 1) return(1)

    # Covariances entre le bloc j et les autres blocs.
    B <- crossprod(do.call(cbind, X[-j]), X[[j]]) / (n - 1)
    u <- drop(PMA::SPC(B, sumabsv = sparsity[j], K = 1,
                       center = FALSE, niter = 1000,
                       trace = FALSE, compute.pve = FALSE)$v)
    u <- u / sqrt(sum(u^2))
    if (u[which.max(abs(u))] < 0) u <- -u
    u
  })

  names(a) <- names(blocks)
  a
}
