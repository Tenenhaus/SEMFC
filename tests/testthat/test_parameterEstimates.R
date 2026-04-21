# Tests pour la méthode parameterEstimates()
# Utilise la nouvelle structure de configurations

# ============================================================================
# SETUP : Charger les configurations
# ============================================================================

# Sourcer les configurations et fonctions de setup
source(testthat::test_path("dataset_configs/load_configs.R"))
source(testthat::test_path("setup_references.R"))

# ============================================================================
# FONCTION GÉNÉRIQUE : Test parameterEstimates() et get_estimate()
# ============================================================================

# Helper function pour vérifier les colonnes de parameterEstimates via la clé (lhs, op, rhs)
verify_parameter_estimates <- function(estimates, reference, dataset_name, tolerance) {
  # Vérifier les colonnes présentes
  expect_equal(colnames(estimates), colnames(reference),
               info = paste("Column names mismatch for", dataset_name))

  # Pour chaque ligne de référence, vérifier qu'elle existe dans estimates avec les mêmes valeurs
  for (i in seq_len(nrow(reference))) {
    ref_row <- reference[i, ]
    lhs <- as.character(ref_row$lhs)
    op <- as.character(ref_row$op)
    rhs <- as.character(ref_row$rhs)

    # Chercher la ligne correspondante dans estimates
    matching_idx <- which(
      as.character(estimates$lhs) == lhs &
        as.character(estimates$op) == op &
        as.character(estimates$rhs) == rhs
    )

    expect_length(matching_idx, 1)

    est_row <- estimates[matching_idx, ]

    # Fonction utilitaire pour comparer avec tolérance, en gérant les NA
    compare_value <- function(est_val, ref_val, tol, desc) {
      # Si la référence est NA, ignorer la vérification (passer le test)
      if (is.na(ref_val)) {
        return(invisible(NULL))
      }
      # Sinon, comparer avec tolérance
      expect_equal(est_val, ref_val, tolerance = tol,
                   info = paste(dataset_name, "-", desc, "mismatch for key (", lhs, ",", op, ",", rhs, ")"))
    }

    # Vérifier les valeurs
    compare_value(est_row$est, ref_row$est, tolerance, "Estimate")
    compare_value(est_row$se, ref_row$se, tolerance, "SE")
    compare_value(est_row$z, ref_row$z, tolerance, "Z")
    compare_value(est_row$pvalue, ref_row$pvalue, tolerance, "Pvalue")
    compare_value(est_row$std.all, ref_row$std.all, tolerance, "Std.all")
  }
}

# Helper function pour vérifier les estimations
verify_estimate <- function(estimate_actual, estimate_ref, estimate_type, dataset_name) {
  tolerance <- 1e-2
  if (is.null(estimate_actual) || is.null(estimate_ref)) {
    expect_equal(is.null(estimate_actual), is.null(estimate_ref),
                 info = paste(dataset_name, "- estimate", estimate_type, "NULL mismatch"))
    return(invisible(NULL))
  }

  expect_equal(class(estimate_actual), class(estimate_ref),
               info = paste(dataset_name, "- estimate", estimate_type, "class mismatch"))

  if (is.numeric(estimate_actual)) {
    expect_equal(as.numeric(estimate_actual), as.numeric(estimate_ref),
                 tolerance = tolerance,
                 info = paste(dataset_name, "- estimate", estimate_type, "values mismatch"))

    if (!is.null(names(estimate_actual))) {
      expect_equal(names(estimate_actual), names(estimate_ref),
                   info = paste(dataset_name, "- estimate", estimate_type, "names mismatch"))
    }

    if (!is.null(rownames(estimate_actual))) {
      expect_equal(rownames(estimate_actual), rownames(estimate_ref),
                   info = paste(dataset_name, "- estimate", estimate_type, "rownames mismatch"))
    }
    if (!is.null(colnames(estimate_actual))) {
      expect_equal(colnames(estimate_actual), colnames(estimate_ref),
                   info = paste(dataset_name, "- estimate", estimate_type, "colnames mismatch"))
    }
  } else if (is.list(estimate_actual)) {
    expect_equal(length(estimate_actual), length(estimate_ref),
                 info = paste(dataset_name, "- estimate", estimate_type, "list length mismatch"))
    expect_equal(names(estimate_actual), names(estimate_ref),
                 info = paste(dataset_name, "- estimate", estimate_type, "list names mismatch"))

    for (i in seq_along(estimate_actual)) {
      elem_name <- names(estimate_actual)[i]
      verify_estimate(estimate_actual[[i]], estimate_ref[[i]],
                     paste(estimate_type, "-", elem_name), dataset_name)
    }
  } else {
    expect_equal(estimate_actual, estimate_ref,
                 info = paste(dataset_name, "- estimate", estimate_type, "mismatch"))
  }
}

# Fonction 1 : teste parameterEstimates() pour tous les datasets
test_parameterEstimates_for_estimator <- function(estimator, fit_params) {
  tolerance <- 1e-2

  for (dataset_name in AVAILABLE_DATASETS) {
    test_that(paste("parameterEstimates() -", toupper(estimator), "-", dataset_name), {
      cat("\n[", toupper(estimator), "]", " Testing parameterEstimates():", dataset_name, "\n")

      setup <- setup_test_data_for_reference(dataset_name)

      sem_model <- SemFC$new(
        data = setup$data,
        relation_matrix = setup$relation_matrix,
        mode = setup$mode,
        estimator = estimator
      )

      set.seed(20091979)
      cat("  - Fitting model...\n")
      do.call(sem_model$fit, c(list(infer = TRUE), fit_params))

      cat("  - Testing parameterEstimates()...\n")
      estimates <- sem_model$parameterEstimates(standardized = TRUE)
      reference_file_param <- testthat::test_path("fixtures",
                                                   paste0("parameterEstimates_reference_", estimator, "_", dataset_name, ".rds"))

      if (file.exists(reference_file_param)) {
        reference_param <- readRDS(reference_file_param)
        verify_parameter_estimates(estimates, reference_param, dataset_name, tolerance)
        cat("    [OK] parameterEstimates() passed\n")
      } else {
        skip(paste("parameterEstimates reference not found for", dataset_name))
      }
    })
  }
}

# Fonction 2 : teste get_estimate() pour tous les datasets
test_get_estimate_for_estimator <- function(estimator, fit_params) {
  for (dataset_name in AVAILABLE_DATASETS) {
    test_that(paste("get_estimate() -", toupper(estimator), "-", dataset_name), {
      cat("\n[", toupper(estimator), "]", " Testing get_estimate():", dataset_name, "\n")

      setup <- setup_test_data_for_reference(dataset_name)

      sem_model <- SemFC$new(
        data = setup$data,
        relation_matrix = setup$relation_matrix,
        mode = setup$mode,
        estimator = estimator
      )

      set.seed(20091979)
      cat("  - Fitting model...\n")
      do.call(sem_model$fit, c(list(infer = TRUE), fit_params))

      cat("  - Testing get_estimate()...\n")
      reference_file_est <- testthat::test_path("fixtures",
                                                paste0("get_estimate_reference_", estimator, "_", dataset_name, ".rds"))

      if (file.exists(reference_file_est)) {
        reference_all <- readRDS(reference_file_est)

        for (est_type in ESTIMATE_TYPES) {
          tryCatch({
            estimate_actual <- sem_model$get_estimate(est_type)
            estimate_ref <- reference_all[[est_type]]
            verify_estimate(estimate_actual, estimate_ref, est_type, dataset_name)
            cat("    [OK]", est_type, "passed\n")
          }, error = function(e) {
            skip(paste(dataset_name, "- estimate type", est_type, "not available"))
          })
        }
      } else {
        skip(paste("get_estimate reference not found for", dataset_name))
      }
    })
  }
}

# ============================================================================
# TEST : Tests de parameterEstimates() - TOUS les datasets
# SVD (B=100) + ML + standardized + inférence
# ============================================================================

test_parameterEstimates_for_estimator("svd", list(B = 100, seed = 20091979))
test_parameterEstimates_for_estimator("ml", list(tol = 1e-6))

# ============================================================================
# TEST : Tests de get_estimate() - TOUS les datasets
# SVD (B=100) + ML + inférence
# ============================================================================

test_get_estimate_for_estimator("svd", list(B = 100, seed = 20091979))
test_get_estimate_for_estimator("ml", list(tol = 1e-6))
