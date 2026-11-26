# M4 : sparse SPC + validation croisee SEMFC ciblee
# Dependances : PMA et la fonction lvm() de semfc (svd_lvm.R).
# A : liste des blocs (n x p_j) ; C : matrice des relations structurelles.
# c_j est une borne L1 absolue, comprise entre 1 et sqrt(p_j).

# 1. Directions SPC, independamment pour chaque bloc
m4_directions <- function(A, c) {
  X <- lapply(A, function(x) scale(x, center = TRUE, scale = FALSE))
  n <- nrow(X[[1]])

  a <- lapply(seq_along(X), function(j) {
    if (ncol(X[[j]]) == 1) return(1)

    # Covariances entre le bloc j et les autres blocs
    B <- crossprod(do.call(cbind, X[-j]), X[[j]]) / (n - 1)
    u <- drop(PMA::SPC(B, sumabsv = c[j], K = 1,
                       center = FALSE, niter = 1000,
                       trace = FALSE, compute.pve = FALSE)$v)
    u <- u / sqrt(sum(u^2))
    if (u[which.max(abs(u))] < 0) u <- -u
    u
  })
  names(a) <- names(A)
  a
}

# 2. Reconstruction SEMFC a partir de la covariance d'apprentissage
m4_reconstruct <- function(S, a, C, mode, idx) {
  J <- length(a)
  p <- nrow(S)
  d <- numeric(J)
  lambda <- vector("list", J)
  omega <- vector("list", J)

  for (j in seq_len(J)) {
    ij <- idx[[j]]
    Sjj <- S[ij, ij, drop = FALSE]
    u <- a[[j]]

    if (mode[j] == "reflective") {
      # Correction identique a celle de svdSEM()
      W <- tcrossprod(u)
      diag(W) <- 0
      WS <- W * Sjj
      pos <- WS > 0
      d[j] <- if (any(pos)) sqrt(sum(WS[pos]) / sum(W[pos]^2)) else 1
    } else {
      # Normalisation des composites formatifs
      v <- solve(Sjj, u)
      d[j] <- 1 / sqrt(drop(crossprod(u, v)))
      omega[[j]] <- d[j] * v
    }
    lambda[[j]] <- d[j] * u
  }

  # Correlations latentes non contraintes
  Phat <- diag(1, J)
  for (j in seq_len(J - 1)) {
    for (k in (j + 1):J) {
      Sjk <- S[idx[[j]], idx[[k]], drop = FALSE]
      rho <- drop(crossprod(a[[j]], Sjk %*% a[[k]])) / (d[j] * d[k])
      Phat[j, k] <- Phat[k, j] <- rho
    }
  }
  dimnames(Phat) <- list(names(a), names(a))

  # Structure SEMFC : prend en charge les relations reciproques
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
  list(sigma_implied = Sigma, a = a, lambda = lambda, d = d,
       omega = omega[mode == "formative"], p_tilde = Phat,
       p_implied = lv$R_LVM, beta = lv$BETA,
       gamma = lv$GAMMA, psi = lv$PSI)
}

# 3. Estimation a c fixe (utilisee dans chaque fold et sur tout l'echantillon)
fit_m4_semfc <- function(A, C, mode, c) {
  p <- vapply(A, ncol, integer(1))
  idx <- split(seq_len(sum(p)), rep(seq_along(p), times = p))
  S <- cov(do.call(cbind, A))
  a <- m4_directions(A, c)
  fit <- m4_reconstruct(S, a, C, mode, idx)
  fit$c <- c
  fit
}

# 4. Perte ciblee : seules les covariances impliquant le bloc j
m4_loss <- function(S_valid, Sigma, idx, j, mode, eta = 1) {
  ij <- idx[[j]]
  pj <- length(ij)
  loss <- 0

  for (k in setdiff(seq_along(idx), j)) {
    ik <- idx[[k]]
    D <- S_valid[ij, ik, drop = FALSE] - Sigma[ij, ik, drop = FALSE]
    loss <- loss + sum(D^2) / (pj * length(ik))
  }

  if (mode[j] == "reflective" && pj > 1) {
    D <- S_valid[ij, ij, drop = FALSE] - Sigma[ij, ij, drop = FALSE]
    diag(D) <- 0
    loss <- loss + eta * sum(D^2) / (pj * (pj - 1))
  }
  loss
}

# 5. M4 : chaque c_j varie, les autres restent a la MEME reference c0
#    La reference par defaut est c0_j = sqrt(p_j) (SPC non penalise).
tune_m4_semfc <- function(A, C, mode, G = 9, nfolds = 5,
                          eta = 1, seed = 123, reference = NULL) {
  p <- vapply(A, ncol, integer(1))
  J <- length(A)
  idx <- split(seq_len(sum(p)), rep(seq_along(p), times = p))
  grids <- lapply(p, function(pj) seq(1, sqrt(pj), length.out = G))
  if (is.null(reference)) reference <- sqrt(p)

  # Memes folds pour tous les candidats et tous les blocs
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = nrow(A[[1]])))
  folds <- lapply(seq_len(nfolds), function(f) {
    train <- which(fold_id != f)
    valid <- which(fold_id == f)
    list(
      A_train = lapply(A, function(x) x[train, , drop = FALSE]),
      S_valid = cov(do.call(cbind, lapply(A, function(x)
        x[valid, , drop = FALSE])))
    )
  })

  selected <- reference
  curves <- vector("list", J)
  for (j in seq_len(J)) {
    scores <- numeric(length(grids[[j]]))
    for (i in seq_along(grids[[j]])) {
      candidate <- reference      # surtout PAS selected
      candidate[j] <- grids[[j]][i]

      scores[i] <- mean(vapply(folds, function(f) {
        fit <- fit_m4_semfc(f$A_train, C, mode, candidate)
        m4_loss(f$S_valid, fit$sigma_implied, idx, j, mode, eta)
      }, numeric(1)))
    }
    selected[j] <- grids[[j]][which.min(scores)]
    curves[[j]] <- data.frame(c = grids[[j]], loss = scores)
  }

  names(selected) <- names(A)
  names(curves) <- names(A)
  list(c = selected, curves = curves, reference = reference)
}
