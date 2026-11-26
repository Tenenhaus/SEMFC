# M4 : chaque regularisation c_j est choisie independamment.
# Les autres blocs restent toujours aux memes regularisations de reference.

tune_sparse_semfc <- function(blocks, model, grids = NULL, grid_size = 9,
                              nfolds = 5, eta = 1, seed = 123,
                              reference = NULL) {
  p <- vapply(blocks, ncol, integer(1))
  J <- length(blocks)
  idx <- split(seq_len(sum(p)), rep(seq_along(p), times = p))

  if (is.null(grids)) {
    grids <- lapply(p, function(pj) seq(1, sqrt(pj), length.out = grid_size))
  }
  if (is.null(reference)) reference <- sqrt(p)

  folds <- make_cv_folds(nrow(blocks[[1]]), nfolds, seed)
  cv_data <- prepare_cv_data(blocks, folds)

  selected <- reference
  curves <- vector("list", J)
  for (j in seq_len(J)) {
    scores <- vapply(grids[[j]], function(cj) {
      candidate <- reference    # Toujours la meme reference, pas 'selected'.
      candidate[j] <- cj
      evaluate_cv(candidate, cv_data, model, target = j,
                  block_indices = idx, eta = eta)
    }, numeric(1))

    selected[j] <- grids[[j]][which.min(scores)]
    curves[[j]] <- data.frame(sparsity = grids[[j]], loss = scores)
  }

  names(selected) <- names(blocks)
  names(curves) <- names(blocks)
  list(sparsity = selected, curves = curves, reference = reference)
}

# M1 : PMA::SPC.cv appliquee separément a chaque matrice de
# covariances croisees. Meme estimateur SPC final que M4.
# La CV de SPC.cv porte sur les entrees de B, pas sur les observations.

tune_spc_cv <- function(blocks, grids = NULL, grid_size = 9,
                        nfolds = 5, rule = "1se", seed = 123,
                        niter = 40) {
  p <- vapply(blocks, ncol, integer(1))
  J <- length(blocks)
  n <- nrow(blocks[[1]])

  if (is.null(grids)) {
    grids <- lapply(p, function(pj) seq(1, sqrt(pj), length.out = grid_size))
  }

  X <- lapply(blocks, function(x)
    scale(as.matrix(x), center = TRUE, scale = FALSE))

  selected <- numeric(J)
  cv_fits <- vector("list", J)

  for (j in seq_len(J)) {
    if (p[j] == 1) {
      selected[j] <- 1
      next
    }

    B <- crossprod(do.call(cbind, X[-j]), X[[j]]) / (n - 1)
    set.seed(seed + j)
    cv <- PMA::SPC.cv(
      x = B, sumabsvs = grids[[j]],
      nfolds = nfolds, center = FALSE,
      niter = niter, trace = FALSE
    )

    selected[j] <- if (rule == "1se") cv$bestsumabsv1se else cv$bestsumabsv
    cv_fits[[j]] <- cv
  }

  names(selected) <- names(blocks)
  names(cv_fits) <- names(blocks)
  list(sparsity = selected, cv_fits = cv_fits)
}
