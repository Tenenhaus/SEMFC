# Standalone smoke tests; no package dependencies.
# Run from the directory containing sparse-cv-data.R and this test file:
# source("sparse-cv-data.R")
# source("test-sparse-cv-data.R")

set.seed(2026)
n <- 100L
z <- rnorm(n)
L <- list(
  factor_a = cbind(a1 = z + rnorm(n), a2 = 2 * z + rnorm(n)),
  factor_b = cbind(b1 = z + rnorm(n), b2 = z + rnorm(n), b3 = rnorm(n)),
  composite_c = matrix(z + rnorm(n), ncol = 1L,
                       dimnames = list(NULL, "c1"))
)

# Reproducibility and RNG state preservation.
state_before <- .Random.seed
folds <- make_cv_folds(n = n, nfolds = 5L, seed = 42)
stopifnot(identical(state_before, .Random.seed))
stopifnot(identical(folds$fold_id,
                    make_cv_folds(n = n, nfolds = 5L, seed = 42)$fold_id))
stopifnot(all(tabulate(folds$fold_id, nbins = 5L) == 20L))

cv_data <- prepare_cv_data(L, folds)
stopifnot(length(cv_data$folds) == 5L)
stopifnot(identical(unname(cv_data$block_sizes), c(2L, 3L, 1L)))
stopifnot(identical(lapply(cv_data$block_indices, length),
                    list(factor_a = 2L, factor_b = 3L, composite_c = 1L)))

for (f in cv_data$folds) {
  Xt <- do.call(cbind, f$train_blocks)
  Xv <- do.call(cbind, lapply(L, function(x) x[f$valid_idx, , drop = FALSE]))
  stopifnot(isTRUE(all.equal(f$S_train, cov(Xt))))
  stopifnot(isTRUE(all.equal(f$S_valid, cov(Xv))))
  stopifnot(length(intersect(f$train_idx, f$valid_idx)) == 0L)
  stopifnot(identical(f$S_train, t(f$S_train)))
  stopifnot(identical(f$S_valid, t(f$S_valid)))
}

# No data leakage: modifying validation rows leaves training covariance unchanged.
f1 <- cv_data$folds[[1L]]
L_changed <- L
L_changed[[1L]][f1$valid_idx, 1L] <-
  L_changed[[1L]][f1$valid_idx, 1L] +
  seq_along(f1$valid_idx) * 100
cv_changed <- prepare_cv_data(L_changed, folds)
stopifnot(isTRUE(all.equal(f1$S_train, cv_changed$folds[[1L]]$S_train)))
stopifnot(!isTRUE(all.equal(f1$S_valid, cv_changed$folds[[1L]]$S_valid)))

# Grouped folds: every subject belongs to a single validation fold.
subjects <- rep(seq_len(25L), each = 4L)
folds_grouped <- make_cv_folds(n, 5L, seed = 42, groups = subjects)
stopifnot(all(vapply(split(folds_grouped$fold_id, subjects),
                     function(x) length(unique(x)) == 1L, logical(1))))
prepare_cv_data(L, folds_grouped)

# The training data can be passed straight to either direction estimator:
# fit_spc_directions(cv_data$folds[[1L]]$train_blocks, sparsity = c(1.3, 1.5, 1))
# fit_sgcca_directions(cv_data$folds[[1L]]$train_blocks, sparsity = c(1.3, 1.5, 1))

message("All CV preparation tests passed.")
