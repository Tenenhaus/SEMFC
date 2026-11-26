# Reconstruction SVD-SEM/SEMFC a partir de directions donnees.
# Reprend la correction des loadings et l'estimation structurelle par lvm().

reconstruct_sparse_semfc <- function(S, directions, model, block_indices) {
  a <- directions
  C <- model$relation_matrix
  mode <- model$mode
  idx <- block_indices
  J <- length(a)
  p <- nrow(S)

  d <- numeric(J)
  lambda <- vector("list", J)
  omega <- vector("list", J)

  for (j in seq_len(J)) {
    Sjj <- S[idx[[j]], idx[[j]], drop = FALSE]
    u <- a[[j]]

    if (mode[j] == "reflective") {
      W <- tcrossprod(u)
      diag(W) <- 0
      WS <- W * Sjj
      pos <- WS > 0
      d[j] <- if (any(pos)) sqrt(sum(WS[pos]) / sum(W[pos]^2)) else 1
    } else {
      v <- solve(Sjj, u)
      d[j] <- 1 / sqrt(drop(crossprod(u, v)))
      omega[[j]] <- d[j] * v
    }
    lambda[[j]] <- d[j] * u
  }

  # Correlations latentes estimees bloc par bloc.
  Phat <- diag(1, J)
  for (j in seq_len(J - 1)) {
    for (k in (j + 1):J) {
      Sjk <- S[idx[[j]], idx[[k]], drop = FALSE]
      rho <- drop(crossprod(a[[j]], Sjk %*% a[[k]])) / (d[j] * d[k])
      Phat[j, k] <- Phat[k, j] <- rho
    }
  }
  dimnames(Phat) <- list(names(a), names(a))

  # Modele structurel existant de SEMFC, y compris liens reciproques.
  lv <- lvm(Phat, C)

  Lambda <- matrix(0, p, J)
  E <- matrix(0, p, p)
  for (j in seq_len(J)) {
    ij <- idx[[j]]
    Lambda[ij, j] <- lambda[[j]]
    Sjj <- S[ij, ij, drop = FALSE]

    if (mode[j] == "reflective") {
      E[ij, ij] <- diag(diag(Sjj) - lambda[[j]]^2, nrow = length(ij))
    } else {
      E[ij, ij] <- Sjj - tcrossprod(lambda[[j]])
    }
  }

  Sigma <- Lambda %*% lv$R_LVM %*% t(Lambda) + E
  names(lambda) <- names(a)
  names(omega) <- names(a)

  list(sigma = Sigma, directions = a, lambda = lambda, d = d,
       omega = omega[mode == "formative"], p_tilde = Phat,
       p_implied = lv$R_LVM, beta = lv$BETA,
       gamma = lv$GAMMA, psi = lv$PSI)
}
