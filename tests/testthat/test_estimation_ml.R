test_that("RML-SEM estimation", {

  setup <- make_test_semfc_data()


  # Construct the theoretical parameter vector
  sigma_blocks <- unlist(lapply(
    list(1:3, 4:6, 7:9, 10:12),
    function(idx) {
      S_jj <- setup$SIGMA[idx, idx]
      S_jj[lower.tri(S_jj, diag = TRUE)]
    }
  ))

  x <- c(
    unlist(setup$lambda, use.names = FALSE),
    setup$PHI[lower.tri(setup$PHI)],
    setup$GAMMA[1, 1:2],
    setup$GAMMA[2, 3:4],
    setup$BETA[2, 1],
    setup$BETA[1, 2],
    setup$R[5, 6],
    sigma_blocks,
    1 - setup$lambda[[5]]^2,
    1 - setup$lambda[[6]]^2
  )



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

})