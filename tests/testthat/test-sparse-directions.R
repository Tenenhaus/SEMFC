# Minimal smoke test. Run after source("sparse-directions.R") and installing
# packages PMA and RGCCA.
set.seed(123)
n <- 150L
z <- matrix(rnorm(n * 3), ncol = 3)
blocks <- list(
  A = tcrossprod(z[, 1], c(0.9, 0.7, 0.0)) + matrix(rnorm(n * 3, sd = 0.4), n, 3),
  B = tcrossprod(z[, 1] + 0.2 * z[, 2], c(0.8, 0.3, 0.0, 0.2)) +
      matrix(rnorm(n * 4, sd = 0.4), n, 4),
  C = tcrossprod(z[, 1] + 0.3 * z[, 3], c(0.7, 0.0, 0.6)) +
      matrix(rnorm(n * 3, sd = 0.4), n, 3)
)
penalties <- c(1.6, 1.8, 1.5)

for (f in list(fit_spc_directions, fit_sgcca_directions)) {
  ans <- f(blocks, sparsity = penalties)
  stopifnot(length(ans$directions) == length(blocks))
  for (j in seq_along(blocks)) {
    u <- ans$directions[[j]]
    stopifnot(length(u) == ncol(blocks[[j]]))
    stopifnot(abs(sum(u^2) - 1) < 1e-6)
    stopifnot(sum(abs(u)) <= penalties[j] + 1e-4)
  }
  print(lapply(ans$directions, round, 3))
}
