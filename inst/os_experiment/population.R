make_population_J <- function(J, q = 3) {

  stopifnot(J >= 6)

  m <- J - 2              # number of formative blocks
  stopifnot(m %% 2 == 0)

  k <- m / 2              #ifnot(m %% 2 == 0)

  k <- m / 2              # formative blocks per endogenous factor

  mode <- c(
    rep("formative", m),
    rep("reflective", 2)
  )

  # ============================================================
  # 1. Structural model
  # ============================================================

  BETA <- matrix(
    c(
      0,   0.25,
      0.5, 0
    ),
    2, 2,
    byrow = TRUE
  )

  R22 <- matrix(
    c(
      1, sqrt(1 / 2),
      sqrt(1 / 2), 1
    ),
    2, 2
  )

  # ------------------------------------------------------------
  # Exogenous correlation matrix
  # ------------------------------------------------------------

  PHI <- matrix(0.5, m, m)
  diag(PHI) <- 1


  # ============================================================
  # 2. GAMMA
  # ============================================================

  GAMMA <- matrix(0, 2, m)

  # Extend the original coefficient patterns
  g1 <- rep(c(-0.30, 0.50), length.out = k)
  g2 <- rep(c( 0.50, 0.25), length.out = k)

  # Preserve the structural signal of the original J = 6 model
  target1 <- 0.19
  target2 <- 0.4375

  PHI1 <- PHI[1:k, 1:k, drop = FALSE]
  PHI2 <- PHI[(k + 1):m, (k + 1):m, drop = FALSE]

  g1 <- g1 * sqrt(
    target1 /
      drop(t(g1) %*% PHI1 %*% g1)
  )

  g2 <- g2 * sqrt(
    target2 /
      drop(t(g2) %*% PHI2 %*% g2)
  )

  GAMMA[1, 1:k] <- g1
  GAMMA[2, (k + 1):m] <- g2


  # ============================================================
  # 3. Disturbance covariance
  # ============================================================

  PSI <-
    (diag(2) - BETA) %*%
    R22 %*%
    t(diag(2) - BETA) -
    GAMMA %*% PHI %*% t(GAMMA)

  # Check positive definiteness
  stopifnot(
    min(
      eigen(
        PSI,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    ) > 0
  )


  # ============================================================
  # 4. Complete latent correlation matrix R
  # ============================================================

  IB_inv <- solve(diag(2) - BETA)

  R <- rbind(
    cbind(
      PHI,
      PHI %*% t(GAMMA) %*% t(IB_inv)
    ),
    cbind(
      IB_inv %*% GAMMA %*% PHI,
      IB_inv %*%
        (GAMMA %*% PHI %*% t(GAMMA) + PSI) %*%
        t(IB_inv)
    )
  )


  # ============================================================
  # 5. Formative blocks
  # ============================================================

  SIGMA_form_base <- matrix(
    c(
      1,   0.3, 0.4,
      0.3, 1,   0.5,
      0.4, 0.5, 1
    ),
    3, 3,
    byrow = TRUE
  )

  SIGMA_form <- replicate(
    m,
    SIGMA_form_base,
    simplify = FALSE
  )


  # ============================================================
  # 6. Formative weights
  # ============================================================

  omega <- lapply(seq_len(m), function(j) {

    # Recover exactly the original model for J = 6:
    #
    # eta1, eta2 : (1,1,1)
    # eta3, eta4 : (1,2,3)

    if (j <= k) {
      w <- rep(1, q)
    } else {
      w <- seq_len(q)
    }

    Sj <- SIGMA_form[[j]]

    w / drop(
      sqrt(
        t(w) %*% Sj %*% w
      )
    )
  })


  # ============================================================
  # 7. Formative loadings
  # ============================================================

  lambda_form <- lapply(
    seq_len(m),
    function(j) {
      drop(
        SIGMA_form[[j]] %*% omega[[j]]
      )
    }
  )


  # ============================================================
  # 8. Reflective blocks
  # ============================================================

  lambda_ref <- list(
    rep(0.7, q),
    rep(0.7, q)
  )

  SIGMA_ref_base <- matrix(0.49, q, q)
  diag(SIGMA_ref_base) <- 1

  SIGMA_ref <- list(
    SIGMA_ref_base,
    SIGMA_ref_base
  )


  # ============================================================
  # 9. Complete loading list
  # ============================================================

  lambda <- c(
    lambda_form,
    lambda_ref
  )

  LAMBDA <- Matrix::bdiag(
    lapply(
      lambda,
      function(x) matrix(x, ncol = 1)
    )
  )


  # ============================================================
  # 10. Observed covariance matrix
  # ============================================================

  SIGMA_blocks <- c(
    SIGMA_form,
    SIGMA_ref
  )

  SIGMA <- Matrix::bdiag(SIGMA_blocks)

  index_end <- cumsum(lengths(lambda))

  index_start <- c(
    1,
    index_end[-length(index_end)] + 1
  )

  range_index <- lapply(
    seq_len(J),
    function(j) {
      index_start[j]:index_end[j]
    }
  )

  for (j in 1:(J - 1)) {

    for (i in (j + 1):J) {

      row_index <- range_index[[j]]
      col_index <- range_index[[i]]

      lj <- lambda[[j]]
      li <- lambda[[i]]

      SIGMA[row_index, col_index] <-
        R[j, i] * lj %*% t(li)

      SIGMA[col_index, row_index] <-
        t(SIGMA[row_index, col_index])
    }
  }


  # ============================================================
  # 11. True parameter vector
  # ============================================================

  true_lambda <- unlist(lambda)

  true_phi <- PHI[upper.tri(PHI)]

  true_gamma <- c(
    GAMMA[1, 1:k],
    GAMMA[2, (k + 1):m]
  )

  true_beta <- c(
    BETA[2, 1],
    BETA[1, 2]
  )

  true_endo_cov <- R[m + 1, m + 2]

  true_sigma_form <- unlist(
    lapply(
      SIGMA_form,
      function(S) {
        S[lower.tri(S, diag = TRUE)]
      }
    )
  )

  true_theta_ref <- unlist(
    lapply(
      lambda_ref,
      function(l) {
        1 - l^2
      }
    )
  )

  true_param_with_S <- c(
    true_lambda,
    true_phi,
    true_gamma,
    true_beta,
    true_endo_cov,
    true_sigma_form,
    true_theta_ref
  )


  # ============================================================
  # Output
  # ============================================================

  list(
    J = J,
    m = m,
    k = k,
    q = q,
    mode = mode,
    BETA = BETA,
    GAMMA = GAMMA,
    PHI = PHI,
    PSI = PSI,
    R = R,
    omega = omega,
    lambda = lambda,
    LAMBDA = LAMBDA,
    SIGMA_blocks = SIGMA_blocks,
    SIGMA = as.matrix(SIGMA),
    true_param_with_S = true_param_with_S
  )
}


generate_mixed_sample <- function(N, J, q, SIGMA, empirical = FALSE) {

  # ============================================================
  # Checks
  # ============================================================

  stopifnot(J >= 4)
  stopifnot(q >= 1)

  n_formative <- J - 2

  # We want to split formative blocks equally between
  # the two reflective endogenous variables
  stopifnot(n_formative %% 2 == 0)

  k <- n_formative / 2
  p <- J * q

  stopifnot(
    nrow(SIGMA) == p,
    ncol(SIGMA) == p
  )


  # ============================================================
  # Variable names
  # ============================================================

  variable_names <- unlist(
    lapply(
      seq_len(J),
      function(j) {
        paste0("X", j, seq_len(q))
      }
    )
  )


  # ============================================================
  # Generate sample
  # ============================================================

  X <- MASS::mvrnorm(
    n = N,
    mu = rep(0, p),
    Sigma = SIGMA,
    empirical = empirical
  )

  colnames(X) <- variable_names


  # ============================================================
  # Split observed variables into J blocks
  # ============================================================

  Y <- lapply(
    seq_len(J),
    function(j) {

      ind <- ((j - 1) * q + 1):(j * q)

      X[
        ,
        ind,
        drop = FALSE
      ]
    }
  )

  names(Y) <- paste0("LV", seq_len(J))


  # ============================================================
  # Structural relation matrix C
  # ============================================================

  C <- matrix(
    0,
    nrow = J,
    ncol = J
  )

  # Indices of the two reflective endogenous blocks
  endo1 <- J - 1
  endo2 <- J


  # First half of formative variables -> first endogenous LV
  C[
    seq_len(k),
    endo1
  ] <- 1


  # Second half of formative variables -> second endogenous LV
  C[
    (k + 1):n_formative,
    endo2
  ] <- 1


  # Reciprocal relation between endogenous variables
  C[endo1, endo2] <- 1
  C[endo2, endo1] <- 1

  rownames(C) <- colnames(C) <- names(Y)


  # ============================================================
  # Block modes
  # ============================================================

  mode <- c(
    rep("formative", n_formative),
    rep("reflective", 2)
  )


  # ============================================================
  # Output
  # ============================================================

  list(
    X = X,
    Y = Y,
    C = C,
    mode = mode,
    J = J,
    q = q,
    p = p,
    n_formative = n_formative,
    n_reflective = 2
  )
}


generate_sem_model_J <- function(J, q = 3) {

  m <- J - 2
  k <- m / 2

  eta_ref1 <- m + 1
  eta_ref2 <- m + 2

  lines <- character()

  # Reflective blocks
  lines <- c(
    lines,
    paste0(
      "eta", eta_ref1, " =~ ",
      paste0(
        "X", eta_ref1, seq_len(q),
        collapse = " + "
      )
    ),
    paste0(
      "eta", eta_ref2, " =~ ",
      paste0(
        "X", eta_ref2, seq_len(q),
        collapse = " + "
      )
    )
  )

  # Formative blocks
  for (j in seq_len(m)) {

    lines <- c(
      lines,
      paste0(
        "eta", j, " <~ ",
        paste0(
          "X", j, seq_len(q),
          collapse = " + "
        )
      )
    )
  }

  # First structural equation
  lines <- c(
    lines,
    paste0(
      "eta", eta_ref1,
      " ~ ",
      paste0(
        "eta", seq_len(k),
        collapse = " + "
      ),
      " + eta", eta_ref2
    )
  )

  # Second structural equation
  lines <- c(
    lines,
    paste0(
      "eta", eta_ref2,
      " ~ ",
      paste0(
        "eta", (k + 1):m,
        collapse = " + "
      ),
      " + eta", eta_ref1
    )
  )

  # Residual covariance
  lines <- c(
    lines,
    paste0(
      "eta", eta_ref1,
      " ~~ eta", eta_ref2
    )
  )

  paste(lines, collapse = "\n")
}