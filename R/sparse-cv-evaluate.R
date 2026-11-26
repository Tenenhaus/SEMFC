# Score CV d'UNE configuration c pour UN bloc cible j.

evaluate_cv <- function(sparsity, cv_data, model, target,
                        block_indices, eta = 1) {
  losses <- vapply(cv_data, function(f) {
    # Le SEM est estime uniquement sur les donnees d'apprentissage.
    a <- fit_spc_directions(f$train_blocks, sparsity)
    fit <- reconstruct_sparse_semfc(f$S_train, a, model, block_indices)

    cv_loss(f$S_valid, fit$sigma, block_indices, target,
            model$mode, eta)
  }, numeric(1))

  mean(losses)
}
