# Script de rapport comparative - SVD vs ML
# Générer un rapport comparatif pour tous les datasets
#
# Utilisation:
#   source("scripts/generate_comparative_report.R")

# ============================================================================
# SETUP - Charger les configurations et fonctions utiles
# ============================================================================
library(devtools)
load_all()
library(lavaan)


pkg_root <- getwd()

# Sourcer les configurations
source(file.path(pkg_root, "tests/testthat/dataset_configs/load_configs.R"))

# Sourcer les fonctions utiles
source(file.path(pkg_root, "inst/report/report_utils.R"))

# Créer le dossier de rapport
report_dir <- file.path(pkg_root, "inst/report/reports")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# GÉNÉRER LE RAPPORT COMPLET
# ============================================================================

# Ouvrir un fichier de sortie pour le rapport complet
report_file <- file.path(report_dir, paste0("rapport_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
file_connection <- file(report_file, "w")
sink(file_connection, append = TRUE)

# Afficher le titre du rapport
cat("\n")
cat("=", strrep("=", 78), "\n", sep = "")
cat("RAPPORT COMPARATIF SVD vs ML\n")
cat(strrep("=", 80), "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Répertoire:", pkg_root, "\n")
cat("Datasets:", length(AVAILABLE_DATASETS), "| Estimations:", length(ESTIMATE_TYPES), "\n\n")

# Générer le rapport pour chaque dataset
all_results <- list()
for (dataset_name in AVAILABLE_DATASETS) {
  tryCatch({
    all_results[[dataset_name]] <- generate_dataset_report(dataset_name)
  }, error = function(e) {
    cat("[ERREUR] Pour", dataset_name, ":", as.character(e), "\n\n")
  })
}

cat("\n")
cat(strrep("=", 80), "\n")
cat("FIN DU RAPPORT\n")
cat(strrep("=", 80), "\n")
cat("Report generated at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Dossier de rapports:", report_dir, "\n\n")

# Fermer la redirection de la sortie
sink()
close(file_connection)


cat("\n")
cat("*** report saved at:", report_file, "***\n")
cat("\n")


