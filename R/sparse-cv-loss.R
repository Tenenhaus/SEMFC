# Perte M4 : blocs croises touchant j + hors-diagonale du bloc j reflectif.

cv_loss <- function(S_valid, sigma_hat, block_indices, target,
                    mode, eta = 1) {
  idx <- block_indices
  j <- target
  ij <- idx[[j]]
  pj <- length(ij)
  loss <- 0

  for (k in setdiff(seq_along(idx), j)) {
    ik <- idx[[k]]
    D <- S_valid[ij, ik, drop = FALSE] - sigma_hat[ij, ik, drop = FALSE]
    loss <- loss + sum(D^2) / (pj * length(ik))
  }

  if (mode[j] == "reflective" && pj > 1) {
    D <- S_valid[ij, ij, drop = FALSE] - sigma_hat[ij, ij, drop = FALSE]
    diag(D) <- 0
    loss <- loss + eta * sum(D^2) / (pj * (pj - 1))
  }

  loss
}

# Evaluation commune M1/M4 sur un echantillon TEST independant.
# Toutes les paires de blocs sont comptees une seule fois.
sem_test_loss <- function(S_test, sigma_hat, block_indices, mode, eta = 1) {
  loss <- 0
  J <- length(block_indices)

  for (j in seq_len(J)) {
    ij <- block_indices[[j]]
    pj <- length(ij)

    for (k in seq_len(J)) {
      if (k <= j) next
      ik <- block_indices[[k]]
      D <- S_test[ij, ik, drop = FALSE] - sigma_hat[ij, ik, drop = FALSE]
      loss <- loss + sum(D^2) / (pj * length(ik))
    }

    if (mode[j] == "reflective" && pj > 1) {
      D <- S_test[ij, ij, drop = FALSE] - sigma_hat[ij, ij, drop = FALSE]
      diag(D) <- 0
      loss <- loss + eta * sum(D^2) / (pj * (pj - 1))
    }
  }

  loss
}
