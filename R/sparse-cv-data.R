# Les memes folds sont utilises pour tous les candidats de regularisation.

make_cv_folds <- function(n, nfolds = 5, seed = 123) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  lapply(seq_len(nfolds), function(k) {
    list(train = which(fold_id != k), valid = which(fold_id == k))
  })
}

prepare_cv_data <- function(blocks, folds) {
  lapply(folds, function(f) {
    A_train <- lapply(blocks, function(x) x[f$train, , drop = FALSE])
    A_valid <- lapply(blocks, function(x) x[f$valid, , drop = FALSE])
    list(
      train_blocks = A_train,
      S_train = cov(do.call(cbind, A_train)),
      S_valid = cov(do.call(cbind, A_valid))
    )
  })
}
