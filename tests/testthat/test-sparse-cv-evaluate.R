# Standalone orchestration tests for sparse-cv-evaluate.R.
# Run in the directory holding sparse-cv-data.R, sparse-cv-loss.R and
# sparse-cv-evaluate.R: Rscript test-sparse-cv-evaluate.R
# No PMA, RGCCA, igraph or Matrix required: estimation and reconstruction
# are replaced by deterministic stand-ins ONLY in this test environment.



set.seed(123)
n <- 60L
z1 <- stats::rnorm(n)
z2 <- 0.7 * z1 + stats::rnorm(n, sd = 0.6)
z3 <- 0.8 * z2 + stats::rnorm(n, sd = 0.6)
L <- list(
  first = cbind(a1 = z1 + rnorm(n, sd = .2),
                a2 = .8 * z1 + rnorm(n, sd = .4)),
  second = cbind(b1 = z2 + rnorm(n, sd = .3),
                 b2 = z2 + rnorm(n, sd = .3)),
  third = cbind(c1 = z3 + rnorm(n, sd = .4),
                c2 = z3 + rnorm(n, sd = .4))
)
cv <- prepare_cv_data(L, make_cv_folds(n, nfolds = 5L, seed = 22L))
c <- setNames(rep(1.2, 3L), names(L))
C <- matrix(0, 3L, 3L, dimnames = list(names(L), names(L)))
C[1, 2] <- 1
C[2, 3] <- 1
model <- list(relation_matrix = C,
              mode = c(first = "reflective", second = "reflective",
                       third = "formative"))

# Deterministic estimator stub: direction depends on TRAINING blocks only.
fit_spc_directions <- function(blocks, sparsity, control = list()) {
  if (isTRUE(control$fail)) stop("forced estimator failure")
  list(directions = lapply(blocks, function(x) {
    u <- rep(1, ncol(x)); u / sqrt(sum(u^2))
  }), method = "spc")
}
fit_sgcca_directions <- function(blocks, sparsity, control = list()) {
  ans <- fit_spc_directions(blocks, sparsity, control)
  ans$method <- "sgcca"
  ans
}

# Stub for the SEM reconstruction; records that the training S is passed in.
seen_S <- list()
reconstruct_sparse_semfc <- function(S, directions, model,
                                     block_indices = NULL,
                                     compute_ls = FALSE, strict = FALSE) {
  seen_S[[length(seen_S) + 1L]] <<- S
  list(sigma = S, admissible = TRUE)
}

# 1. Spectral block-target score agrees with an explicit fold-by-fold loop.
result <- evaluate_cv(c, cv, estimator = "spc", criterion = "spectral",
                      target = "first", details = TRUE)
manual <- vapply(cv$folds, function(f) {
  u <- rep(1 / sqrt(2), 2)
  cv_loss(f$S_train, f$S_valid, cv$block_indices,
          directions = list(first = u, second = u, third = u),
          criterion = "spectral", target = 1L)
}, numeric(1))
stopifnot(
  isTRUE(all.equal(unname(result$fold_losses), unname(manual))),
  isTRUE(all.equal(result$mean_loss, mean(manual))),
  isTRUE(all.equal(result$se_loss, stats::sd(manual) / sqrt(5))),
  result$n_failed == 0L,
  length(result$fold_details) == 5L
)

# 2. Complete separable objective is the sum of independent block scores.
all_spectral <- evaluate_cv(c, cv, "spc", "spectral")
by_block <- vapply(seq_along(L), function(j) {
  evaluate_cv(c, cv, "spc", "spectral", target = j)$mean_loss
}, numeric(1))
stopifnot(isTRUE(all.equal(all_spectral$mean_loss, sum(by_block))))

# 3. SGCCA estimator path works for a bilateral diagnostic.
sgcca <- evaluate_cv(c, cv, estimator = "sgcca", criterion = "bilateral")
stopifnot(is.finite(sgcca$mean_loss), sgcca$n_failed == 0L)

# 4. SEM paths use precisely S_train; compare with direct loss.
seen_S <- list()
global <- evaluate_cv(c, cv, "spc", "global", model = model,
                      eta = 1, keep_fits = TRUE)
stopifnot(
  global$n_failed == 0L,
  length(seen_S) == length(cv$folds),
  length(global$fold_fits) == length(cv$folds),
  all(global$fold_admissible)
)
for (ell in seq_along(cv$folds)) {
  stopifnot(isTRUE(all.equal(seen_S[[ell]], cv$folds[[ell]]$S_train)))
  expected <- cv_loss(
    cv$folds[[ell]]$S_train, cv$folds[[ell]]$S_valid,
    cv$block_indices, criterion = "global", sigma_hat = seen_S[[ell]],
    eta = 1, reflective = c(TRUE, TRUE, FALSE)
  )
  stopifnot(isTRUE(all.equal(global$fold_losses[[ell]], expected)))
}

targeted <- evaluate_cv(c, cv, "sgcca", "targeted", model = model,
                        target = "second", eta = 1)
stopifnot(is.finite(targeted$mean_loss), targeted$target == 2L)

# 5. Folding failures: Inf, not an average over only successful folds.
failed <- evaluate_cv(c, cv, "spc", "spectral", control = list(fail = TRUE),
                      on_error = "inf")
stopifnot(is.infinite(failed$mean_loss), failed$n_failed == 5L,
          all(is.infinite(failed$fold_losses)), all(!is.na(failed$fold_errors)),
          is.na(failed$se_loss))

throws <- function(expr) inherits(
  tryCatch(force(expr), error = identity), "error")
stopifnot(
  throws(evaluate_cv(c, cv, "spc", "spectral",
                     control = list(fail = TRUE), on_error = "stop")),
  throws(evaluate_cv(c, cv, "spc", "targeted", model = model)),
  throws(evaluate_cv(c, cv, "spc", "global")),
  throws(evaluate_cv(c, cv, "spc", "global", model = model, eta = -1)),
  throws(evaluate_cv(rep(10, 3), cv, "spc", "spectral")),
  throws(evaluate_cv(c, cv, "spc", "global", model = model, target = 1))
)

# 6. Named sparsities are properly reordered.
reordered <- evaluate_cv(rev(c), cv, "spc", "spectral", target = 1)
stopifnot(identical(names(reordered$sparsity), names(L)),
          isTRUE(all.equal(unname(reordered$sparsity), unname(c))))

message("All evaluate_cv() orchestration tests passed.")
