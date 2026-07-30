# Configuration dataset model_mixed
# Source: inst/model/model_mixed.R

# Créer les données pour model_mixed
create_model_mixed_dataset <- function() {
  library(MASS)
  library(mvtnorm)
  set.seed(20091979)
  N <- 300

  # Load SIGMA from data_simulation_mixed.R
  source(system.file("simulations", "data_simulation_mixed.R", package = "SEMFC"))

  # Generate dataset
  X <- mvrnorm(N, rep(0, 18), SIGMA, empirical = FALSE)
  colnames(X) <- paste("X", rep(1:6, each = 3), rep(1:3, 6), sep = "")

  # Create list of latent variables
  A <- list(
    LV1 = X[, 1:3],
    LV2 = X[, 4:6],
    LV3 = X[, 7:9],
    LV4 = X[, 10:12],
    LV5 = X[, 13:15],
    LV6 = X[, 16:18]
  )

  return(A)
}

# Configuration de model_mixed
create_model_mixed_config <- function() {
  # Récupérer les données générées
  A <- create_model_mixed_dataset()

  # Relation matrix (C) - Path structure
  C <- matrix(c(0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 1, 0), 6, 6, byrow = TRUE)

  colnames(C) <- rownames(C) <- names(A)

  # Mode: formative for first 4, reflective for last 2
  mode <- c(rep("formative", 4), rep("reflective", 2))
  names(mode) <- names(A)


  sem <-  '
    # latent variable definitions
    LV1 =~ X11+X12+X13
    LV2 =~ X21+X22+X23
    LV3 =~ X31+X32+X33
    LV4 =~ X41+X42+X43
    LV5 =~ X51+X52+X53
    LV6 =~ X61+X62+X63

    # Regressions
    LV5 ~ LV1 + LV2 + LV6
    LV6 ~ LV3 + LV4 + LV5

    # residual covariances
    LV5 ~~ LV6
    '

  return(list(data = A, relation_matrix = C, mode = mode, sem = sem))
}

