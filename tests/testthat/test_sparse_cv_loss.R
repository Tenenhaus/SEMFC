# Base-R smoke tests for sparse-cv-loss.R; no PMA / RGCCA required.
# Run from the directory containing sparse-cv-data.R and sparse-cv-loss.R:
# source("sparse-cv-data.R")
# source("sparse-cv-loss.R")
# source("test-sparse-cv-loss.R")

set.seed(37)
n <- 80L
z <- rnorm(n)
blocks <- list(
  factor_a = cbind(a1 = z + rnorm(n), a2 = 2*z + rnorm(n)),
  factor_b = cbind(b1 = z + rnorm(n), b2 = -z + rnorm(n), b3 = rnorm(n)),
  composite_c = matrix(z + rnorm(n), ncol = 1L)
)
cv_data <- prepare_cv_data(blocks, make_cv_folds(n, nfolds = 4L, seed = 42L))
fold <- cv_data$folds[[1L]]
S_train <- fold$S_train
S_valid <- fold$S_valid
idx <- cv_data$block_indices

U <- list(
  factor_a = c(1, 0),
  factor_b = c(0, 1, 0),
  composite_c = 1
)

# Spectral loss is separable, and both signs of a unit direction are equivalent.
fj <- vapply(seq_along(idx), function(j) {
  cv_loss(S_train, S_valid, idx, U, "spectral", target = j)
}, numeric(1))
Fsep <- cv_loss(S_train, S_valid, idx, U, "spectral")
stopifnot(isTRUE(all.equal(Fsep, sum(fj))))
U_flip <- U
U_flip[[1L]] <- -U_flip[[1L]]
stopifnot(isTRUE(all.equal(Fsep, cv_loss(S_train, S_valid, idx,
                                        U_flip, "spectral"))))
stopifnot(isTRUE(all.equal(fj[2L], cv_loss(S_train, S_valid, idx,
                                           U, "spectral", target = "factor_b"))))

# The one-sided separable loss upper-bounds the bilateral rank-one loss.
Fbil <- cv_loss(S_train, S_valid, idx, U, "bilateral")
stopifnot(Fbil <= Fsep + 1e-10)

# Explicit matrix implementation, for the same fold (no shortcut involved).
pair <- function(j, k) {
  A <- S_valid[idx[[j]], idx[[k]], drop = FALSE]
  B <- S_train[idx[[j]], idx[[k]], drop = FALSE]
  P <- tcrossprod(U[[j]])
  sum((A - P %*% B)^2) / (length(idx[[j]]) * length(idx[[k]]))
}
stopifnot(isTRUE(all.equal(fj[1L], pair(1, 2) + pair(1, 3))))

# SEM criteria: same full sigma_hat reused for targeted and global.
# Deliberately use the training covariance as a stand-in for a prediction;
# this tests the loss module, not the SEM reconstruction module.
sigma_hat <- S_train
eta <- 0.6
refl <- c(TRUE, TRUE, FALSE)
Fglobal <- cv_loss(S_train, S_valid, idx, criterion = "global",
                   sigma_hat = sigma_hat, eta = eta,
                   reflective = refl)
targeted <- vapply(seq_along(idx), function(j) {
  cv_loss(S_train, S_valid, idx, criterion = "targeted",
          target = j, sigma_hat = sigma_hat, eta = eta,
          reflective = refl)
}, numeric(1))

# Global interblock score counts each unordered pair once. Sum of targeted
# interblock scores counts each twice, but within-block scores only once.
detail <- cv_loss(S_train, S_valid, idx, criterion = "global",
                  sigma_hat = sigma_hat, eta = eta,
                  reflective = refl, details = TRUE)
stopifnot(isTRUE(all.equal(Fglobal, detail$loss)))
stopifnot(isTRUE(all.equal(sum(targeted),
  2 * sum(detail$between) + sum(detail$within))))
stopifnot(detail$within["composite_c"] == 0)

# Manual check of the reflective within-block contribution.
i <- idx[[1L]]
D <- S_valid[i, i, drop = FALSE] - sigma_hat[i, i, drop = FALSE]
diag(D) <- 0
manual_within <- eta * sum(D^2) / (length(i) * (length(i)-1))
stopifnot(isTRUE(all.equal(unname(detail$within["factor_a"]), manual_within)))

# Renaming / reordering directions does not change the result.
stopifnot(isTRUE(all.equal(
  cv_loss(S_train, S_valid, idx, U[3:1], "spectral"), Fsep)))

# Invalid inputs should fail informatively.
expect_error <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
stopifnot(expect_error(cv_loss(S_train, S_valid, idx, criterion = "global")))
stopifnot(expect_error(cv_loss(S_train, S_valid, idx, U,
                               criterion = "spectral", eta = 1)))
stopifnot(expect_error(cv_loss(S_train, S_valid, idx,
                               criterion = "targeted", target = 1,
                               sigma_hat = sigma_hat, eta = 1)))
message("All sparse CV loss tests passed.")
