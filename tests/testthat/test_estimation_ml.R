check_ml_sem <- function(recursive) {
  setup <- make_test_semfc_data(recursive = recursive)


  x <- setup$theta_true



  # Initialize the RML-SEM model
  model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    mode = setup$mode,
    estimator = "ml"
  )

  expect_true(inherits(model, "SemFC"))

  # Fit without statistical inference
  expect_no_error(model$fit(infer = FALSE))

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

  # Check theoretical loadings
  for (j in 1:6) {
    expect_equal(
      as.numeric(lambda_hat[[j]]),
      as.numeric(setup$lambda[[j]]),
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

  # Check the implied covariance matrix
  expect_equal(
    unname(Sigma_hat),
    unname(setup$SIGMA),
    tolerance = 1e-10
  )

}
test_that("RML-SEM estimation - nonrecursive", {
  check_ml_sem(recursive = FALSE)
})

test_that("RML-SEM estimation - recursive", {
  check_ml_sem(recursive = TRUE)
})