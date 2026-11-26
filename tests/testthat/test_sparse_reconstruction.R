# Regression tests for sparse-reconstruction.R against the supplied svdSEM().
# Run from the directory containing:
#   sparse-reconstruction.R, sparse-cv-data.R, sparse-cv-loss.R,
#   svd_lvm.R, svd_estimation.R
#
# Rscript test-sparse-reconstruction.R
# Requires igraph and Matrix for the original svdSEM implementation.

source("svd_lvm.R")
source("svd_estimation.R")
source("sparse-reconstruction.R")
source("sparse-cv-data.R")
source("sparse-cv-loss.R")

if (!requireNamespace("igraph", quietly = TRUE) ||
    !requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install igraph and Matrix to run the regression tests.")
}

# Minimal stand-ins for package functions not included in the supplied files.
# Only used here (scale = FALSE, bias = FALSE).
if (!exists("scale2", mode = "function")) {
  scale2 <- function(x, bias = FALSE, scale = TRUE) {
    base::scale(as.matrix(x), center = TRUE, scale = scale)
  }
}
if (!exists("cov2", mode = "function")) {
  cov2 <- function(x, y = NULL, bias = FALSE) {
    x <- as.matrix(x)
    ans <- if (is.null(y)) stats::cov(x) else stats::cov(x, as.matrix(y))
    if (bias) ans <- ans * ((nrow(x) - 1) / nrow(x))
    ans
  }
}
if (!exists("d_LS", mode = "function")) {
  d_LS <- function(S, sigma) sum((S - sigma)^2)
}

set.seed(125)
n <- 450L
z1 <- rnorm(n)
z2 <- 0.7 * z1 + sqrt(1 - 0.7^2) * rnorm(n)
z3 <- 0.6 * z2 + sqrt(1 - 0.6^2) * rnorm(n)
A <- list(
  A = cbind(a1 = z1 + rnorm(n, sd = .3),
            a2 = 0.8 * z1 + rnorm(n, sd = .4),
            a3 = 1.2 * z1 + rnorm(n, sd = .5)),
  B = cbind(b1 = z2 + rnorm(n, sd = .5),
            b2 = 1.1 * z2 + rnorm(n, sd = .4),
            b3 = 0.6 * z2 + rnorm(n, sd = .6)),
  C = cbind(c1 = z3 + rnorm(n, sd = .5),
            c2 = 0.75 * z3 + rnorm(n, sd = .6))
)
C <- matrix(0, 3, 3, dimnames = list(names(A), names(A)))
C[1, 2] <- 1
C[2, 3] <- 1
mode <- c("reflective", "formative", "reflective")
model <- list(relation_matrix = C, mode = mode)
S <- stats::cov(do.call(cbind, A))

old <- svdSEM(A, C, mode = mode, scale = FALSE, bias = FALSE)
new <- reconstruct_sparse_semfc(S, old$a, model, compute_ls = TRUE)

check_equal <- function(x, y, tol = 1e-7) {
  isTRUE(all.equal(unname(x), unname(y), tolerance = tol,
                   check.attributes = FALSE))
}
stopifnot(
  check_equal(new$sigma, old$sigma_implied),
  check_equal(new$d, old$d),
  check_equal(new$p_tilde, old$p_tilde),
  check_equal(new$p_implied, old$p_implied),
  check_equal(new$beta, old$beta),
  check_equal(new$gamma, old$gamma),
  check_equal(new$T_LS, old$T_LS)
)
for (j in seq_along(A)) {
  stopifnot(check_equal(new$lambda[[j]], old$lambda[[j]]))
}

# The formative diagonal block is left equal to its empirical covariance;
# for reflective blocks, the model reconstructs the observed variances.
idx <- split(seq_len(ncol(S)),
             rep(seq_along(A), times = vapply(A, ncol, integer(1))))
stopifnot(check_equal(new$sigma[idx[[2L]], idx[[2L]], drop = FALSE],
                      S[idx[[2L]], idx[[2L]], drop = FALSE]))
stopifnot(check_equal(diag(new$sigma), diag(S)))

# Flipping the sign of a latent direction must not change the implied Sigma.
flipped <- old$a
flipped[[1L]] <- -flipped[[1L]]
flipped_sem <- reconstruct_sparse_semfc(S, flipped, model)
stopifnot(check_equal(flipped_sem$sigma, new$sigma))

# Reconstruct from the fold training covariance only and pass to cv_loss().
cv <- prepare_cv_data(A, make_cv_folds(n, nfolds = 5L, seed = 8L))
fold <- cv$folds[[1L]]
train_fit <- svdSEM(fold$train_blocks, C, mode = mode,
                    scale = FALSE, bias = FALSE)
sem <- reconstruct_sparse_semfc(
  S = fold$S_train,
  directions = train_fit$a,
  model = model,
  block_indices = cv$block_indices
)
stopifnot(check_equal(sem$sigma, train_fit$sigma_implied))

loss_target <- cv_loss(
  S_train = fold$S_train, S_valid = fold$S_valid,
  block_indices = cv$block_indices, criterion = "targeted",
  target = 1L, sigma_hat = sem$sigma, eta = 1,
  reflective = mode == "reflective"
)
loss_global <- cv_loss(
  S_train = fold$S_train, S_valid = fold$S_valid,
  block_indices = cv$block_indices, criterion = "global",
  sigma_hat = sem$sigma, eta = 1,
  reflective = mode == "reflective"
)
stopifnot(is.finite(loss_target), is.finite(loss_global),
          loss_target >= 0, loss_global >= 0)

# Alter directions to check the module actually uses the supplied directions.
modified <- old$a
modified[[1L]] <- c(1, 0, 0)
other_sem <- reconstruct_sparse_semfc(S, modified, model)
stopifnot(!check_equal(other_sem$sigma, new$sigma))

# Incorrect direction dimensions and missing structural model must fail.
expect_error <- function(expr) inherits(
  tryCatch(force(expr), error = identity), "error")
stopifnot(expect_error(reconstruct_sparse_semfc(S, old$a[1:2], model)))
stopifnot(expect_error(reconstruct_sparse_semfc(S, old$a, list(mode = mode))))

message("All sparse SEMFC reconstruction regression tests passed.")
