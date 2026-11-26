# Ajustement a regularisations fixees (fold ou echantillon complet).

fit_sparse_semfc <- function(blocks, model, sparsity) {
  p <- vapply(blocks, ncol, integer(1))
  idx <- split(seq_len(sum(p)), rep(seq_along(p), times = p))
  S <- cov(do.call(cbind, blocks))
  a <- fit_spc_directions(blocks, sparsity)

  fit <- reconstruct_sparse_semfc(S, a, model, idx)
  fit$sparsity <- sparsity
  fit
}
