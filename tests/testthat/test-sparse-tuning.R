# Standalone tests for sparse-tuning.R / sparse-fit.R.
# Run from the folder containing the files:
#   Rscript test-sparse-tuning.R
# These tests do not require PMA, RGCCA, igraph, or SEMFC internals;
# they use deterministic stubs ONLY for the integration boundaries.
# Actual package integration must be tested separately in R.



throws <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")

# 1. Separable exact grid minimization and constant reference across targets.
g <- list(A = c(1, 1.5, 2), B = c(1, 1.5, 2), C = c(1, 1.5, 2))
reference <- c(A = 2, B = 2, C = 2)
opt <- c(A = 1, B = 1.5, C = 2)
seen <- list()
obj <- function(c, target) {
  seen[[length(seen) + 1L]] <<- list(c = c, target = target)
  j <- target
  list(mean_loss = (c[j] - opt[j])^2, se_loss = 0.01)
}
sel <- select_sparsity(obj, g, reference, rule = "min")
stopifnot(identical(unname(sel$sparsity), unname(opt)),
          sel$n_evaluations == sum(lengths(g)))
for (entry in seen) {
  stopifnot(isTRUE(all.equal(unname(entry$c[-entry$target]),
                             unname(reference[-entry$target]))))
}

# 2. 1-SE chooses smallest feasible c (sparsest candidate).
flat <- select_sparsity(function(c, target) {
  list(mean_loss = if (c[target] == 1) 1.08 else 1,
       se_loss = 0.1)
}, g, reference, rule = "1se")
stopifnot(all(flat$sparsity == 1))
stopifnot(throws(select_sparsity(function(c, target) 1,
                                 g, rule = "1se")))

# 3. Coordinate search: coupled example with a coordinatewise trap.
coupled_grids <- list(x = c(1, 2), y = c(1, 2))
objective <- function(c) {
  if (c["x"] == 1 && c["y"] == 1) return(3)
  if (c["x"] == 2 && c["y"] == 2) return(1)
  4
}
local <- coordinate_search(objective, coupled_grids,
                           initial = c(x = 1, y = 1), n_starts = 1)
global_found <- coordinate_search(objective, coupled_grids,
                                  n_starts = 3)
stopifnot(local$value == 3, global_found$value == 1,
          all(global_found$sparsity == 2),
          global_found$n_evaluations <= 4,
          all(diff(global_found$history$loss) <= 1e-10))

# 4. Same fold assignment for all methods + routing of 5 new CV procedures.
set.seed(21)
n <- 40
z <- rnorm(n)
L <- list(
  A = cbind(a1 = z + rnorm(n, sd = .3), a2 = z + rnorm(n, sd = .3),
            a3 = z + rnorm(n, sd = .3), a4 = z + rnorm(n, sd = .3)),
  B = cbind(b1 = z + rnorm(n, sd = .3), b2 = z + rnorm(n, sd = .3),
            b3 = z + rnorm(n, sd = .3), b4 = z + rnorm(n, sd = .3)),
  C = cbind(c1 = z + rnorm(n, sd = .3), c2 = z + rnorm(n, sd = .3),
            c3 = z + rnorm(n, sd = .3), c4 = z + rnorm(n, sd = .3))
)
connection <- matrix(0, 3, 3, dimnames = list(names(L), names(L)))
connection[1, 2] <- connection[2, 3] <- 1
model <- list(relation_matrix = connection,
              mode = c("reflective", "reflective", "reflective"))
folds <- make_cv_folds(n, nfolds = 4L, seed = 5)
mock_calls <- list()
evaluate_cv <- function(sparsity, cv_data, estimator, criterion,
                        model, target = NULL, eta = 0,
                        normalization = "block", control = list(),
                        strict = FALSE, on_error = "inf") {
  mock_calls[[length(mock_calls) + 1L]] <<- list(
    c = sparsity, target = target, criterion = criterion,
    estimator = estimator)
  target_opt <- c(A = 1, B = 1.5, C = 2)
  if (criterion %in% c("spectral", "targeted")) {
    score <- (sparsity[target] - target_opt[target])^2
  } else {
    score <- sum((sparsity - target_opt)^2) +
      .4 * (sparsity[1L] - sparsity[2L])^2
  }
  list(mean_loss = unname(score), se_loss = 0.01)
}

methods <- list(
  M3 = c("spc", "spectral"),
  M4 = c("spc", "targeted"),
  M5 = c("spc", "global"),
  M6 = c("sgcca", "targeted"),
  M7 = c("sgcca", "global")
)
out <- lapply(methods, function(m) tune_sparse_semfc(
  blocks = L, model = model, estimator = m[1L], criterion = m[2L],
  grids = g, folds = folds, rule = "min",
  search_control = list(n_starts = 3, max_passes = 10)
))
stopifnot(identical(names(out), names(methods)),
          all(vapply(out, function(z) all(z$sparsity %in% c(1,1.5,2)), logical(1))),
          all(vapply(out, function(z) identical(z$fold_id, folds$fold_id), logical(1))),
          all(out$M3$sparsity == opt), all(out$M4$sparsity == opt),
          all(out$M6$sparsity == opt),
          identical(out$M5$method, "M5"), identical(out$M7$method, "M7"),
          throws(tune_sparse_semfc(L, model, "sgcca", "spectral")),
          throws(tune_sparse_semfc(L, model, "spc", "global", rule = "1se")))

# 5. Final-fit orchestration: downstream covariance must only use given blocks.
fit_spc_directions <- function(blocks, sparsity, control = list()) {
  U <- lapply(blocks, function(b) rep(1 / sqrt(ncol(b)), ncol(b)))
  list(directions = U, sparsity = sparsity, method = "spc")
}
fit_sgcca_directions <- function(blocks, sparsity, control = list()) {
  out <- fit_spc_directions(blocks, sparsity, control)
  out$method <- "sgcca"
  out
}
reconstruct_sparse_semfc <- function(S, directions, model,
                                     block_indices = NULL, compute_ls = FALSE,
                                     strict = FALSE) {
  list(sigma = S, parameters = list(), admissible = TRUE,
       diagnostics = list(min_eigenvalue_sigma = min(eigen(S)$values)))
}
final <- fit_sparse_semfc(L, model, estimator = "spc",
                          sparsity = out$M3)
stopifnot(isTRUE(all.equal(final$sigma,
                           stats::cov(do.call(cbind, L)))),
          all(final$sparsity == out$M3$sparsity),
          final$admissible,
          throws(fit_sparse_semfc(L, model, estimator = "sgcca",
                                  sparsity = out$M3)))
message("All sparse tuning and final-fit orchestration tests passed.")
