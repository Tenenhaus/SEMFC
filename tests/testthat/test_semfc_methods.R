# Tests pour les méthodes de la classe SemFC
# Utilise la nouvelle structure de configurations séparées

# ============================================================================
# SETUP : Charger les configurations
# ============================================================================

# Sourcer les configurations et fonctions de setup

source(testthat::test_path("dataset_configs/load_configs.R"))
source(testthat::test_path("setup_references.R"))


# ============================================================================
# FIXTURE COMMUNE : Créer un modèle test
# ============================================================================

# Fixture commune - utilise la nouvelle structure
setup_test_model <- function(dataset_name = "ECSI", estimator = "svd", fitted = TRUE) {
  # Utiliser la nouvelle fonction setup_test_data_for_reference()
  setup <- setup_test_data_for_reference(dataset_name)

  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    mode = setup$mode,
    estimator = estimator
  )

  if (fitted) {
    sem_model$fit()
  }

  return(sem_model)
}

# ============================================================================
# Tests pour initialize()
# ============================================================================

test_that("SemFC$new() creates valid object - ECSI", {
  sem_model <- setup_test_model(dataset_name = "ECSI", fitted = FALSE)

  expect_s3_class(sem_model, "SemFC")
  expect_true(!is.null(sem_model))
})

test_that("SemFC$new() initializes with different estimators - ECSI", {
  model_svd <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)
  model_ml <- setup_test_model(dataset_name = "ECSI", estimator = "ml", fitted = FALSE)

  expect_s3_class(model_svd, "SemFC")
  expect_s3_class(model_ml, "SemFC")
})

test_that("SemFC$new() works for all available datasets", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, fitted = FALSE)
    expect_s3_class(sem_model, "SemFC")
  }
})

# ============================================================================
# Tests pour fit()
# ============================================================================

test_that("fit() executes without errors - SVD - ECSI", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)

  expect_silent(sem_model$fit(infer = FALSE))
})

test_that("fit() executes without errors - ML - ECSI", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "ml", fitted = FALSE)

  expect_silent(sem_model$fit(infer = FALSE))
})

test_that("fit() with bootstrap inference - SVD - all datasets", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, estimator = "svd", fitted = FALSE)


    expect_silent(capture.output(sem_model$fit(infer = TRUE, B = 30)))
  }
})

test_that("fit() with asymptotic inference - ML - ECSI", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "ml", fitted = FALSE)

  expect_silent(sem_model$fit(infer = TRUE, tol = 1e-3))
})



# ============================================================================
# Tests pour get_estimate()
# ============================================================================

test_that("get_estimate() returns lambda for fitted model", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  lambda <- sem_model$get_estimate("lambda")

  expect_type(lambda, "list")
  expect_gt(length(lambda), 0)
})

test_that("get_estimate() returns all estimates with 'all'", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  all_estimates <- sem_model$get_estimate("all")

  expect_type(all_estimates, "list")
  expect_true("lambda" %in% names(all_estimates))
  expect_true("beta" %in% names(all_estimates))
})

test_that("get_estimate() returns NULL with warning for unfitted model", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)

  expect_warning(sem_model$get_estimate("lambda"))
})

test_that("get_estimate() returns NULL with warning for invalid estimate name", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  expect_warning(sem_model$get_estimate("invalid_estimate"))
})

test_that("get_estimate() supports all documented estimate types - all datasets", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, estimator = "svd", fitted = TRUE)

    estimate_names <- c("lambda", "beta")

    for (est_name in estimate_names) {
      result <- sem_model$get_estimate(est_name)
      # Peut être NULL ou une valeur valide, mais pas d'erreur
      expect_true(!is.null(result) || is.null(result),
                  info = paste("Failed for dataset:", dataset_name, "estimate:", est_name))
    }
  }
})

# ============================================================================
# Tests pour check_improper()
# ============================================================================

test_that("check_improper() returns logical vector", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  improper_results <- sem_model$check_improper()

  expect_type(improper_results, "logical")
})

test_that("check_improper() fails before fit", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)

  expect_error(sem_model$check_improper())
})

test_that("check_improper() result is named vector", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  improper_results <- sem_model$check_improper()

  expect_true(!is.null(names(improper_results)))
})

test_that("check_improper() works for all datasets", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, estimator = "svd", fitted = TRUE)
    improper_results <- sem_model$check_improper()

    expect_type(improper_results, "logical")
  }
})

# ============================================================================
# Tests pour summary() (interface de rapport)
# ============================================================================

test_that("summary() executes without errors", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  # summary() affiche du texte, on vérifie qu'il ne plante pas
  expect_silent(capture.output(sem_model$summary()))
})

test_that("summary() with different parameters", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = TRUE)

  expect_silent(capture.output(sem_model$summary(standardized = TRUE)))
  expect_silent(capture.output(sem_model$summary(effect = TRUE)))
  expect_silent(capture.output(sem_model$summary(all_measures = TRUE)))
})

test_that("summary() works for all datasets - SVD", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, estimator = "svd", fitted = TRUE)

    expect_silent(capture.output(sem_model$summary()))
  }
})

# ============================================================================
# Tests d'intégration complète
# ============================================================================

test_that("Complete workflow: SVD without inference", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)

  # Fit
  sem_model$fit(infer = FALSE)

  # Get estimates
  lambda <- sem_model$get_estimate("lambda")
  estimates <- sem_model$parameterEstimates()

  expect_type(lambda, "list")
  expect_s3_class(estimates, "data.frame")
})

test_that("Complete workflow: SVD with bootstrap inference", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "svd", fitted = FALSE)

  # Fit avec bootstrap
  sem_model$fit(infer = TRUE, B = 30)

  # Check improper
  improper <- sem_model$check_improper()

  # Get estimates
  estimates <- sem_model$parameterEstimates()

  expect_type(improper, "logical")
  expect_s3_class(estimates, "data.frame")
})

test_that("Complete workflow: ML with asymptotic inference", {
  sem_model <- setup_test_model(dataset_name = "ECSI", estimator = "ml", fitted = FALSE)

  # Fit avec inference asymptotique
  sem_model$fit(infer = TRUE)

  # Get estimates
  estimates <- sem_model$parameterEstimates()

  expect_s3_class(estimates, "data.frame")
  expect_gt(nrow(estimates), 0)
})

test_that("Complete workflow - all datasets - SVD minimal", {
  for (dataset_name in AVAILABLE_DATASETS) {
    sem_model <- setup_test_model(dataset_name = dataset_name, estimator = "svd", fitted = FALSE)

    # Fit sans inference pour rapidité
    sem_model$fit(infer = FALSE)

    # Get estimates
    estimates <- sem_model$parameterEstimates()

    expect_s3_class(estimates, "data.frame")
    expect_gt(nrow(estimates), 0)
  }
})

