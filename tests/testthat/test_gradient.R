test_that("Analytical Fisher information matches numerical information", {

  setup <- make_test_semfc_data()

  init <- get_parameter_model_sem(
    data = setup$data,
    mode = setup$mode,
    relation_matrix = setup$relation_matrix,
    bias = FALSE,
    scale = FALSE
  )

  model <- init$model

  # Construct the theoretical parameter vector
  sigma_blocks <- unlist(lapply(
    list(1:3, 4:6, 7:9, 10:12),
    function(idx) {
      S_jj <- setup$SIGMA[idx, idx]
      S_jj[upper.tri(S_jj, diag = TRUE)]
    }
  ))

  x <- c(
    unlist(setup$lambda, use.names = FALSE),
    setup$PHI[upper.tri(setup$PHI)],
    setup$GAMMA[1, 1:2],
    setup$GAMMA[2, 3:4],
    setup$BETA[1, 2],
    setup$BETA[2, 1],
    setup$R[5, 6],
    sigma_blocks,
    1 - setup$lambda[[5]]^2,
    1 - setup$lambda[[6]]^2
  )

  # Analytical Fisher information
  I_analytical <- information_matrix(x, model)

  # Numerical Jacobian of the upper triangular covariance
  JAC <- numDeriv::jacobian(
    lvm_ml,
    x = x,
    model = model,
    jac = TRUE
  )

  p <- sum(model$block_sizes)
  nb_param <- length(x)

  # Reconstruct symmetric covariance derivative matrices
  full_jac <- lapply(seq_len(ncol(JAC)), function(j) {

    D <- matrix(0, p, p)
    D[upper.tri(D, diag = TRUE)] <- JAC[, j]

    D + t(D) - diag(diag(D))
  })

  # Inverse of the implied covariance matrix
  Sigma <- lvm_ml(x, model, jac = FALSE)$sigma_implied
  Sinv <- solve(Sigma)

  # Numerical Fisher information
  Sinv_full_jac <- lapply(
    full_jac,
    function(D) Sinv %*% D
  )

  I_ij <- function(i, j) {
    0.5 * sum(diag(
      Sinv_full_jac[[i]] %*% Sinv_full_jac[[j]]
    ))
  }

  index_upper <- which(
    upper.tri(
      matrix(0, nb_param, nb_param),
      diag = TRUE
    ),
    arr.ind = TRUE
  )

  I_upper <- mapply(
    I_ij,
    index_upper[, "row"],
    index_upper[, "col"]
  )

  I_numerical <- matrix(0, nb_param, nb_param)
  I_numerical[upper.tri(I_numerical, diag = TRUE)] <- I_upper
  I_numerical <- I_numerical + t(I_numerical) -
    diag(diag(I_numerical))

  # Compare analytical and numerical information matrices
  expect_equal(
    unname(as.matrix(I_analytical)),
    unname(I_numerical),
    tolerance = 1e-10
  )




  # Compute the sample covariance matrix
  S <- cov(setup$X)

  # Analytical constraint Jacobian
  H_analytical <- as.matrix(
    grad_heq(x, S, model)
  )

  # Numerical constraint Jacobian
  H_numerical <- numDeriv::jacobian(
    func = heq,
    x = x,
    S = S,
    model = model
  )

  # Check dimensions
  expect_equal(dim(H_analytical), c(4,length(x)))
  expect_equal(dim(H_analytical), dim(H_numerical))

  # Compare analytical and numerical Jacobians
  expect_equal(
    unname(H_analytical),
    unname(H_numerical),
    tolerance = 1e-10
  )

})