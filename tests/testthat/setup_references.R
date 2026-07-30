# Script pour créer les fichiers de référence RDS
# À exécuter une seule fois pour générer les données de référence
# Maintenant supporte PLUSIEURS datasets !

# Créer le dossier fixtures s'il n'existe pas
dir.create("tests/testthat/fixtures", showWarnings = FALSE)

# Charger les configurations des datasets et variables globales
source(testthat::test_path("dataset_configs/load_configs.R"))

# ============================================================================
# FIXTURE FLEXIBLE POUR LES RÉFÉRENCES
# Utilise les configurations externes pour une meilleure lisibilité
# ============================================================================



# Charger les données du dataset source (si NULL, les données sont générées par la config)
load_dataset_data <- function(dataset_name) {
  data_source <- DATASET_DATA_SOURCE[[dataset_name]]

  if (is.null(data_source)) {
    return(NULL)  # Les données seront générées par la config
  }

  data(list = data_source, envir = environment())
  return(get(data_source, envir = environment()))
}

setup_test_data_for_reference <- function(dataset_name = "ECSI") {
  # Charger les données si nécessaire
  dataset <- load_dataset_data(dataset_name)

  # Récupérer la configuration
  config <- get_dataset_config(dataset_name, dataset)

  return(list(
    data = config$data,
    relation_matrix = config$relation_matrix,
    mode = config$mode,
    dataset_name = dataset_name
  ))
}

# ============================================================================
# Générer les références pour parameterEstimates() ET get_estimate()
# Estimateur: "svd" (B=100) ou "ml" (tol=1e-6)
# ============================================================================
create_reference <- function(dataset_name = "ECSI", estimator = "svd") {
  setup <- setup_test_data_for_reference(dataset_name)

  # Créer le modèle
  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    mode = setup$mode,
    estimator = estimator
  )
  set.seed(20091979)
  # Fit avec paramètres spécifiques à l'estimateur
  if (estimator == "svd") {
    sem_model$fit(infer = TRUE, B = 100, seed = 20091979)
  } else if (estimator == "ml") {
    sem_model$fit(infer = TRUE, tol = 1e-6)
  } else {
    stop("Estimator must be 'svd' or 'ml'")
  }

  # ============================================================
  # 1. Sauvegarder parameterEstimates
  # ============================================================
  estimates <- sem_model$parameterEstimates(standardized = TRUE)

  output_path_est <- paste0(
    "tests/testthat/fixtures/parameterEstimates_reference_",
    estimator,
    "_",
    dataset_name,
    ".rds"
  )
  saveRDS(estimates, output_path_est)
  cat(" parameterEstimates created:", output_path_est, "\n")

  # ============================================================
  # 2. Sauvegarder get_estimate() pour tous les types
  # ============================================================
  all_estimates <- list()
  for (est_type in ESTIMATE_TYPES) {
    tryCatch({
      all_estimates[[est_type]] <- sem_model$get_estimate(est_type)
    }, error = function(e) {
      # Silencieusement ignorer si le type n'existe pas
      NULL
    }, warning = function(w) {
      NULL
    })
  }

  output_path_get <- paste0(
    "tests/testthat/fixtures/get_estimate_reference_",
    estimator,
    "_",
    dataset_name,
    ".rds"
  )
  saveRDS(all_estimates, output_path_get)
  cat(" get_estimate created:", output_path_get, "\n")

  return(list(
    parameterEstimates = estimates,
    get_estimate = all_estimates
  ))
}

# ============================================================================
# Générer TOUS les fichiers de référence
# ============================================================================
create_all_references <- function(estimators = c("svd", "ml")) {
  cat(" Création de toutes les références...\n\n")

  for (dataset_name in AVAILABLE_DATASETS) {
    cat(sprintf("%-30s: ", dataset_name))

    tryCatch({
      if ("svd" %in% estimators) {
        create_reference(dataset_name, estimator = "svd")
        cat("SVD ok ")
      }

      if ("ml" %in% estimators) {
        create_reference(dataset_name, estimator = "ml")
        cat("ML ok")
      }
      cat("\n")

    }, error = function(e) {
      cat("Error:", as.character(e), "\n")
    })
  }

  cat("\n Done !\n")
}
