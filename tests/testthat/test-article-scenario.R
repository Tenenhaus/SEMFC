# Base-R scenario test; no PMA/RGCCA/igraph required.
# Run from extracted benchmark folder:
#   Rscript tests/test-article-scenario.R
args <- commandArgs(trailingOnly = FALSE)
flag <- grep("^--file=", args, value = TRUE)
root <- if (length(flag)) {
  dirname(dirname(normalizePath(sub("^--file=", "", flag[[1L]]))))
} else getwd()
source(file.path(root, "article-scenario.R"))
s <- article_scenario()
stopifnot(
  length(s$lambda) == 6L,
  identical(unname(s$block_sizes), c(60L, 3L, 3L, 3L, 3L, 3L)),
  length(s$strong) == 20L, length(s$null) == 20L,
  length(s$weak) == 20L,
  all(s$lambda[[1]][s$strong] == 0.9),
  all(s$lambda[[1]][s$null] == 0),
  all(s$lambda[[1]][s$weak] == 0.2),
  all(unlist(s$lambda[-1L]) == 0.7),
  identical(dim(s$Sigma), c(75L, 75L)),
  identical(dim(s$R), c(6L, 6L)),
  all(s$min_eigenvalues > 0),
  max(abs(diag(s$Sigma) - 1)) < 1e-10,
  s$C[1, 5] == 1, s$C[2, 5] == 1,
  s$C[3, 6] == 1, s$C[4, 6] == 1,
  s$C[5, 6] == 1, s$C[6, 5] == 1,
  sum(s$C) == 6
)

# Gaussian samples have the intended block dimensions and covariance.
L <- simulate_article_blocks(12000L, s, seed = 123L)
stopifnot(
  length(L) == 6L,
  identical(vapply(L, ncol, integer(1)), s$block_sizes),
  all(vapply(L, nrow, integer(1)) == 12000L),
  max(abs(stats::cov(do.call(cbind, L)) - s$Sigma)) < 0.08
)
message("PASS: article population, indices, structural graph and simulation.")
