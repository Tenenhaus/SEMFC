
make_test_semfc_data <- function(N = 300, seed = 123, measurement = c("mixed", "reflective"), recursive = FALSE) {
  measurement <- match.arg(measurement)
  set.seed(seed)

  if (measurement == "mixed") {
    mode <- c(rep("formative", 4), rep("reflective", 2))
  } else {
    mode <- rep("reflective", 6)
  }

  # Structural relation matrix
  C <- matrix(0, 6, 6)

  C[1, 5] <- 1
  C[2, 5] <- 1
  C[3, 6] <- 1
  C[4, 6] <- 1
  C[5, 6] <- 1
  C[6, 5] <- 1

  BETA <- matrix(c(
    0,   0.25,
    0.5, 0
  ), 2, 2, byrow = TRUE)

  GAMMA <- matrix(c(
    -0.30, 0.50, 0,    0,
     0,    0,   0.50, 0.25
  ), 2, 4, byrow = TRUE)

  PHI <- matrix(0.5, 4, 4)
  diag(PHI) <- 1

  R22 <- matrix(c(
    1, sqrt(0.5),
    sqrt(0.5), 1
  ), 2, 2)

  if (recursive) {
    BETA[1, 2] <- 0
    C[6, 5] <- 0

    R22 <- matrix(
      c(1,     0.575,
        0.575, 1),
      nrow = 2,
      byrow = TRUE
    )
  }

  A <- diag(2) - BETA
  A_inv <- solve(A)

  PSI <- A %*% R22 %*% t(A) -
    GAMMA %*% PHI %*% t(GAMMA)

  R <- rbind(
    cbind(
      PHI,
      PHI %*% t(GAMMA) %*% t(A_inv)
    ),
    cbind(
      A_inv %*% GAMMA %*% PHI,
      A_inv %*%
        (GAMMA %*% PHI %*% t(GAMMA) + PSI) %*%
        t(A_inv)
    )
  )

  # Within-block covariance matrices
  S_form <- matrix(c(
    1,   0.3, 0.4,
    0.3, 1,   0.5,
    0.4, 0.5, 1
  ), 3, 3, byrow = TRUE)

  S_refl <- matrix(0.49, 3, 3)
  diag(S_refl) <- 1


  if (measurement == "mixed") {

    normalize_weight <- function(w) {
      as.vector(
        w / sqrt(drop(t(w) %*% S_form %*% w))
      )
    }

    w1 <- w2 <- normalize_weight(rep(1, 3))
    w3 <- w4 <- normalize_weight(1:3)

    # Theoretical loadings
    lambda <- list(
      as.vector(S_form %*% w1),
      as.vector(S_form %*% w2),
      as.vector(S_form %*% w3),
      as.vector(S_form %*% w4),
      rep(0.7, 3),
      rep(0.7, 3)
    )

    omega <- list(w1, w2, w3, w4)

    sigma_blocks <- c(
      rep(list(S_form), 4),
      rep(list(S_refl), 2)
    )

  } else {

    omega <- NULL

    lambda <- rep(list(rep(0.7, 3)), 6)

    sigma_blocks <- rep(list(S_refl), 6)



  }

  # Construct the population covariance matrix
  SIGMA <- matrix(0, 18, 18)

  indices <- split(1:18, rep(1:6, each = 3))

  for (j in 1:6) {
    SIGMA[indices[[j]], indices[[j]]] <- sigma_blocks[[j]]
  }

  for (j in 1:5) {
    for (k in (j + 1):6) {
      block <- R[j, k] *
        tcrossprod(lambda[[j]], lambda[[k]])

      SIGMA[indices[[j]], indices[[k]]] <- block
      SIGMA[indices[[k]], indices[[j]]] <- t(block)
    }
  }

  # Generate data with exact sample covariance
  X <- MASS::mvrnorm(
    n = N,
    mu = rep(0, 18),
    Sigma = SIGMA,
    empirical = TRUE
  )

  colnames(X) <- paste0(
    "X",
    rep(1:6, each = 3),
    rep(1:3, 6)
  )

  Y <- setNames(
    lapply(indices, function(idx) {
      X[, idx, drop = FALSE]
    }),
    paste0("LV", 1:6)
  )



  dimnames(C) <- list(names(Y), names(Y))

 # Construct the theoretical parameter vector
  # Identify formative and reflective blocks
  formative_blocks <- which(mode == "formative")
  reflective_blocks <- which(mode == "reflective")

  # Within-block covariance parameters for formative blocks
  within_cov_blocks <- unlist(
    lapply(formative_blocks, function(j) {
      idx <- indices[[j]]
      S_jj <- SIGMA[idx, idx, drop = FALSE]
      S_jj[lower.tri(S_jj, diag = TRUE)]
    }),
    use.names = FALSE
  )

  # Measurement error variances for reflective blocks
  error_variances <- unlist(
    lapply(reflective_blocks, function(j) {
      idx <- indices[[j]]
      diag(SIGMA[idx, idx, drop = FALSE]) -
        as.numeric(lambda[[j]])^2
    }),
    use.names = FALSE
  )

  # Structural coefficients
  beta_free <- if (recursive) {
    BETA[2, 1]
  } else {
    c(BETA[2, 1], BETA[1, 2])
  }

  # Additional latent correlation in nonrecursive models
  rho_free <- if (recursive) {
    numeric(0)
  } else {
    R[5, 6]
  }

  # Construct the theoretical parameter vector
  theta_true <- c(
    unlist(lambda, use.names = FALSE),
    PHI[lower.tri(PHI)],
    GAMMA[1, 1:2],
    GAMMA[2, 3:4],
    beta_free,
    rho_free,
    within_cov_blocks,
    error_variances
  )

  return(list(
    theta_true = theta_true,
    data = Y,
    X = X,
    relation_matrix = C,
    mode = mode,
    SIGMA = SIGMA,
    lambda = lambda,
    omega = omega,
    BETA = BETA,
    GAMMA = GAMMA,
    PHI = PHI,
    PSI = PSI,
    R = R
  ))
}
