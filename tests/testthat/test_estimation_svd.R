
check_svd_sem <- function(measurement = measurement, recursive = recursive) {
  setup <- make_test_semfc_data(measurement = measurement, recursive = recursive)
  x <- setup$theta_true

  # Verify the exact empirical covariance
  expect_equal(
    unname(cov(setup$X)),
    unname(setup$SIGMA),
    tolerance = 1e-10
  )

  # Initialize the model
  model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    mode = setup$mode,
    estimator = "svd"
  )



  # Fit with statistical inference
  expect_no_error(model$fit(infer = TRUE, B = 50))
  # Print the summary of the fitted model
  expect_no_error(model$summary(standardized = TRUE, all_measures = TRUE, effect = TRUE))

  expect_no_error(model$parameterEstimates(standardized = TRUE))

  expect_no_error(model$get_gof('all'))

  # Retrieve estimated parameters
  theta_hat <- model$get_estimate("theta")
  lambda_hat <- model$get_estimate("lambda")
  beta_hat <- model$get_estimate("beta")
  gamma_hat <- model$get_estimate("gamma")
  psi_hat <- model$get_estimate("psi")
  Sigma_hat <- as.matrix(
    model$get_estimate("sigma_implied")
  )

  expect_equal(
    unname(as.numeric(theta_hat)),
    unname(x),
    tolerance = 1e-10
  )

  expect_type(lambda_hat, "list")
  expect_length(lambda_hat, 6)

  # Compare with theoretical loadings
  for (j in seq_along(setup$lambda)) {
    expect_equal(
      as.numeric(lambda_hat[[j]]),
      setup$lambda[[j]],
      tolerance = 1e-10,
      info = paste("Incorrect loadings for block", j)
    )
  }


  # Check structural coefficients
  expect_equal(
    unname(as.matrix(beta_hat)),
    unname(setup$BETA),
    tolerance = 1e-10
  )

  expect_equal(
    unname(as.matrix(gamma_hat)),
    unname(setup$GAMMA),
    tolerance = 1e-10
  )

  expect_equal(
    unname(as.matrix(psi_hat)),
    unname(setup$PSI),
    tolerance = 1e-10
  )


  # Indices of the four formative blocks
  block_indices <- list(
    1:3,
    4:6,
    7:9,
    10:12
  )
  if (measurement == 'mixed'){
    # Check identification constraints
    for (j in 1:4) {

      idx <- block_indices[[j]]

      Sigma_jj <- Sigma_hat[idx, idx]
      lambda_j <- as.numeric(lambda_hat[[j]])

      constraint <- drop(
        crossprod(
          lambda_j,
          solve(Sigma_jj, lambda_j)
        )
      )

      expect_equal(
        constraint,
        1,
        tolerance = 1e-10,
        info = paste(
          "Identification constraint failed for block", j
        )
      )
    }
  }

  # Check the implied covariance matrix

  expect_equal(
    unname(Sigma_hat),
    unname(setup$SIGMA),
    tolerance = 1e-10
  )
}

test_that("SVD-SEM estimation - nonrecursive - mixed", {
  check_svd_sem(measurement = "mixed", recursive = FALSE)
})


test_that("SVD-SEM estimation - recursive - mixed", {
  check_svd_sem(measurement = "mixed", recursive = TRUE)
})


test_that("SVD-SEM estimation - nonrecursive - reflective", {
  check_svd_sem(measurement = "reflective",recursive = FALSE)
})


test_that("SVD-SEM estimation - recursive - reflective", {
  check_svd_sem(measurement = "reflective",recursive = TRUE)
})