# Helper pour charger les configurations des datasets
# Sourçe automatiquement tous les fichiers de configuration
# Configuration globale des tests SEMFC
# À sourcer dans tous les scripts de test et setup



# ============================================================================
# TYPES D'ESTIMATIONS POUR get_estimate()
# ============================================================================
ESTIMATE_TYPES <- c(
  "lambda",
  "std_lambda",
  "omega",
  "beta",
  "gamma",
  "p_exo",
  "p_endo",
  "psi",
  "p_implied",
  # "sigma_implied",
  "r2",
  "T_LS",
  # "theta",
  "effect"
)


# ============================================================================
# DICTIONNAIRE CENTRAL DE CONFIGURATION DES DATASETS
# ============================================================================
DATASETS_CONFIG <- list(
  "ECSI" = list(
    data_source = "ECSI",
    config_func = "create_ecsi_config",
    is_generated = FALSE
  ),
  "ECSI_full" = list(
    data_source = "ECSI",
    config_func = "create_ecsi_full_config",
    is_generated = FALSE
  ),
  "BergamiBagozzi2000" = list(
    data_source = "BergamiBagozzi2000",
    config_func = "create_bergami_config",
    is_generated = FALSE
  ),
  "ITFlex" = list(
    data_source = "ITFlex",
    config_func = "create_itflex_config",
    is_generated = FALSE
  ),
  "PoliticalDemocracy" = list(
    data_source = "PoliticalDemocracy",
    config_func = "create_politicaldemocracy_config",
    is_generated = FALSE
  ),
  "Russett" = list(
    data_source = "Russett",
    config_func = "create_russett_config",
    is_generated = FALSE
  ),
  "LancelotMiltgenetal2016" = list(
    data_source = "LancelotMiltgenetal2016",
    config_func = "create_lancelot_config",
    is_generated = FALSE
  ),
  "model_mixed" = list(
    data_source = NULL,
    config_func = "create_model_mixed_config",
    is_generated = TRUE
  ),
  "model_reflective" = list(
    data_source = NULL,
    config_func = "create_model_reflective_config",
    is_generated = TRUE
  )
)

# Dériver les variables utiles du dictionnaire
AVAILABLE_DATASETS <- names(DATASETS_CONFIG)
DATASET_DATA_SOURCE <- setNames(
  sapply(DATASETS_CONFIG, function(x) x$data_source),
  AVAILABLE_DATASETS
)
GENERATED_DATASETS <- names(DATASETS_CONFIG)[sapply(DATASETS_CONFIG, function(x) x$is_generated)]


# Dossier des configurations
config_dir <- testthat::test_path("dataset_configs")

# Sourcer tous les fichiers de configuration (sauf celui-ci!)
.load_dataset_configs <- function() {
  config_files <- list.files(
    config_dir,
    pattern = "\\.R$",
    full.names = TRUE
  )

  # Chemin complet du fichier courant
  this_file <- normalizePath(testthat::test_path("dataset_configs/load_configs.R"))

  for (file in config_files) {
    # Ignorer load_configs.R lui-même pour éviter la récursion
    if (normalizePath(file) != this_file) {
      source(file, local = FALSE)
    }
  }
}

# Fonction utilitaire pour obtenir la config d'un dataset
get_dataset_config <- function(dataset_name, dataset = NULL) {
  if (!dataset_name %in% AVAILABLE_DATASETS) {
    stop("Dataset '", dataset_name, "' not configured.")
  }

  config <- DATASETS_CONFIG[[dataset_name]]
  func <- get(config$config_func)

  # Pour les datasets générés, on ne passe pas dataset en paramètre
  if (config$is_generated) {
    return(func())
  } else {
    return(func(dataset))
  }
}

# Charger les configs au démarrage
.load_dataset_configs()
