# Real integration smoke test for the 7 sparse SEMFC tuning methods.
# Run from a directory containing the seven sparse-*.R files and svd_lvm.R.
# Requires PMA, RGCCA, and igraph to be installed.
# This test intentionally uses few folds, grid points and permutations.

needed <- c("PMA", "RGCCA", "igraph")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R dependencies: ", paste(missing, collapse = ", "))
}


set.seed(123)
n <- 240L
z1 <- rnorm(n)
z2 <- .65 * z1 + sqrt(1-.65^2) * rnorm(n)
z3 <- .6 * z2 + sqrt(1-.6^2) * rnorm(n)
mk <- function(z, prefix) {
  out <- cbind(z + rnorm(n, sd = .35), .9 * z + rnorm(n, sd = .45),
               .7 * z + rnorm(n, sd = .55))
  colnames(out) <- paste0(prefix, 1:3)
  out
}
L <- list(A = mk(z1, "a"), B = mk(z2, "b"), C = mk(z3, "c"))
C <- matrix(0, 3, 3, dimnames = list(names(L), names(L)))
C[1, 2] <- C[2, 3] <- 1
model <- list(relation_matrix = C,
              mode = rep("reflective", 3))
folds <- make_cv_folds(n, nfolds = 3L, seed = 5L)
methods <- list(
  M1 = c("spc", "spc_cv"),
  M2 = c("sgcca", "permutation"),
  M3 = c("spc", "spectral"),
  M4 = c("spc", "targeted"),
  M5 = c("spc", "global"),
  M6 = c("sgcca", "targeted"),
  M7 = c("sgcca", "global")
)
for (m in names(methods)) {
  cat("Testing", m, "...\n")
  typ <- methods[[m]]
  out <- tune_sparse_semfc(
    blocks = L, model = model, estimator = typ[1L],
    criterion = typ[2L], grid_size = 3L,
    nfolds = 3L, folds = folds, seed = 5L, eta = if (m %in%
      c("M4", "M5", "M6", "M7")) 1 else 0,
    baseline_control = if (m == "M2") list(n_perms = 3L) else
      list(niter = 10L),
    search_control = list(n_starts = 2L, max_passes = 3L)
  )
  stopifnot(identical(out$method, m), all(is.finite(out$sparsity)))
  fit <- fit_sparse_semfc(L, model, estimator = typ[1L],
                          sparsity = out)
  stopifnot(is.matrix(fit$sigma), all(is.finite(fit$sigma)),
            nrow(fit$sigma) == sum(vapply(L, ncol, integer(1))))
}
message("Seven-method integration smoke test passed.")
