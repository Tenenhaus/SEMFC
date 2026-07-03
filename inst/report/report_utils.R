# Utilitaires pour la génération de rapports comparatifs
# Fonctions d'aide pour le script generate_comparative_report.R

# ============================================================================
# FONCTION UTILITAIRE - Traiter les estimations
# ============================================================================

process_estimate <- function(est) {
  # Si c'est une liste de vecteurs, la dérouler
  if (is.list(est) && !is.data.frame(est) && !is.matrix(est)) {
    # Vérifier si tous les éléments sont des vecteurs atomiques
    all_vectors <- all(sapply(est, function(x) is.atomic(x) && length(x) > 0))
    if (all_vectors) {
      est <- unlist(est)
    }
  }
  return(est)
}

# ============================================================================
# FONCTION - Fit et récupérer les estimations
# ============================================================================

fit_and_get_estimates <- function(config, estimator, fit_params = list()) {
  fit_label <- paste0(toupper(estimator), " (",
    if (estimator == "svd") paste0("B=", fit_params$B) else paste0("tol=", fit_params$tol),
    ")")
  cat("[INFO] Fitting ", fit_label, "...\n")

  sem <- SemFC$new(
    data = config$data,
    relation_matrix = config$relation_matrix,
    mode = config$mode,
    estimator = estimator
  )
  set.seed(20091979)
  # Fit avec les paramètres appropriés
  if (estimator == "svd") {
    sem$fit(infer = TRUE, B = fit_params$B, verbose = FALSE, seed = 20091979)
  } else if (estimator == "ml") {
    sem$fit(infer = TRUE, tol = fit_params$tol)
  }

  # Récupérer les estimations
  param_est <- sem$parameterEstimates(standardized = TRUE)
  estimates <- list()
  for (est_type in ESTIMATE_TYPES) {
    tryCatch({
      est <- sem$get_estimate(est_type)

      # Si c'est effect, séparer les deux matrices
      if (est_type == "effect") {
        estimates[["total_effect"]] <- est$total_effect
        estimates[["indirect_effect"]] <- est$indirect_effect
      } else {
        estimates[[est_type]] <- process_estimate(est)
      }
    }, error = function(e) {
      NULL
    }, warning = function(w) {
      NULL
    })
  }

  return(list(sem = sem, param = param_est, estimates = estimates))
}

# ============================================================================
# FONCTION - Afficher les parameterEstimates
# ============================================================================

display_parameter_estimates <- function(param_est_svd, param_est_ml) {
  cat("\n", strrep("-", 80), "\n")
  cat("PARAMETER ESTIMATES - SVD (B=100)\n")
  cat(strrep("-", 80), "\n")
  print(param_est_svd)

  cat("\n", strrep("-", 80), "\n")
  cat("PARAMETER ESTIMATES - ML (tol=1e-8)\n")
  cat(strrep("-", 80), "\n")
  print(param_est_ml)
}

# ============================================================================
# FONCTION - Comparer les estimations pour un type donné
# ============================================================================

compare_estimate_type <- function(est_type, est_svd, est_ml) {
  cat("\n[", est_type, "]\n")


  if (is.numeric(est_svd) && is.numeric(est_ml) && !is.matrix(est_svd)) {
    if (length(est_svd) == length(est_ml)) {

      row_names <- names(est_svd)
      if (is.null(row_names)) {
        row_names <- rownames(est_svd)
      }
      if (is.null(row_names)) {
        row_names <- paste0("V", seq_along(est_svd))
      }

      comparison <- data.frame(
        SVD = as.numeric(est_svd),
        ML = as.numeric(est_ml),
        Difference = round(as.numeric(est_svd) - as.numeric(est_ml),3)
      )
      rownames(comparison) <- row_names
      print(comparison)
    }
  }
  # Matrices
  else if (is.matrix(est_svd) && is.matrix(est_ml)) {
    cat("SVD:\n")
    print(est_svd)
    cat("\nML:\n")
    print(est_ml)
    cat("\nDifference:\n")
    print(round(est_svd - est_ml, 3))
  }
  # Listes
  else if (is.list(est_svd) || is.list(est_ml)) {
    svd_names <- if (is.null(names(est_svd))) "NULL" else paste(names(est_svd), collapse = ", ")
    ml_names <- if (is.null(names(est_ml))) "NULL" else paste(names(est_ml), collapse = ", ")
    cat("SVD names:", svd_names, "\n")
    cat("ML names:", ml_names, "\n")
  }
}

# ============================================================================
# FONCTION - Comparer tous les types d'estimations
# ============================================================================

compare_all_estimates <- function(get_est_svd, get_est_ml) {
  cat("\n", strrep("-", 80), "\n")
  cat("COMPARAISON DES ESTIMATIONS (SVD vs ML)\n")
  cat(strrep("-", 80), "\n")

  # Types à comparer (exclure "effect" car il est séparé en total_effect et indirect_effect)
  est_types_to_compare <- ESTIMATE_TYPES[ESTIMATE_TYPES != "sigma_implied" & ESTIMATE_TYPES != "effect"]

  for (est_type in est_types_to_compare) {
    if (!is.null(get_est_svd[[est_type]]) || !is.null(get_est_ml[[est_type]])) {
      est_svd <- get_est_svd[[est_type]]
      est_ml <- get_est_ml[[est_type]]
      compare_estimate_type(est_type, est_svd, est_ml)
    }
  }

  # Comparer total_effect et indirect_effect
  for (effect_type in c("total_effect", "indirect_effect")) {
    if (!is.null(get_est_svd[[effect_type]]) || !is.null(get_est_ml[[effect_type]])) {
      est_svd <- get_est_svd[[effect_type]]
      est_ml <- get_est_ml[[effect_type]]
      compare_estimate_type(effect_type, est_svd, est_ml)
    }
  }
}

# ============================================================================
# FONCTION UTILITAIRE - Charger les données d'un dataset
# ============================================================================

load_dataset_data <- function(dataset_name) {
  data_source <- DATASET_DATA_SOURCE[[dataset_name]]

  if (is.null(data_source)) {
    return(NULL)  # Les données seront générées par la config
  }

  data(list = data_source, envir = environment())
  return(get(data_source, envir = environment()))
}

# ============================================================================
# FONCTION PRINCIPALE - Générer le rapport pour un dataset
# ============================================================================

generate_dataset_report <- function(dataset_name) {
  cat("\n", strrep("=", 140), "\n")
  cat("DATASET:", dataset_name, "\n")
  cat(strrep("=", 140), "\n\n")

  # Charger les données si nécessaire
  dataset <- load_dataset_data(dataset_name)

  # Récupérer la configuration
  config <- get_dataset_config(dataset_name, dataset)

  # Fit et récupérer les estimations
  svd_results <- fit_and_get_estimates(config, "svd", list(B = 100))
  ml_results <- fit_and_get_estimates(config, "ml", list(tol = 1e-8))

  # Afficher les parameterEstimates
  display_parameter_estimates(svd_results$param, ml_results$param)

  # si config$sem est non NULL
    if (!is.null(config$sem)) {
      lavaan_ml <- sem(config$sem,
                       data = data.frame(Reduce("cbind", config$data)),
                       estimator = "ML",
                       likelihood = "wishart")
      lavaan_results <- parameterEstimates(lavaan_ml, standardized = TRUE)

      cat("\n", strrep("-", 80), "\n")
      cat(" Comparison with lavaan \n")
      cat(strrep("-", 80), "\n")

      print(merge_estimates_with_lavaan(ml_results$param, lavaan_results))

    }

  # Comparer les estimations
  compare_all_estimates(svd_results$estimates, ml_results$estimates)

  cat("\n")
  return(list(
    svd = list(param = svd_results$param, estimates = svd_results$estimates),
    ml = list(param = ml_results$param, estimates = ml_results$estimates)
  ))
}

# ============================================================================
# FONCTION - Fusionner les estimations ML et lavaan
# ============================================================================

merge_estimates_with_lavaan <- function(ml_param, lavaan_param) {
  # Créer des clés de fusion basées sur lhs, op, rhs
  ml_param$merge_key <- paste(ml_param$lhs, ml_param$op, ml_param$rhs, sep = "___")
  lavaan_param$merge_key <- paste(lavaan_param$lhs, lavaan_param$op, lavaan_param$rhs, sep = "___")

  # Garder seulement les colonnes nécessaires de lavaan (référence)
  lavaan_subset <- lavaan_param[, c("merge_key", "lhs", "op", "rhs", "std.all")]
  colnames(lavaan_subset)[colnames(lavaan_subset) == "std.all"] <- "std.all.lavaan"
  lavaan_subset$std.all.lavaan <- round(lavaan_subset$std.all.lavaan, 3)

  # Garder seulement les colonnes nécessaires de ML
  ml_subset <- ml_param[, c("merge_key", "std.all")]
  colnames(ml_subset)[colnames(ml_subset) == "std.all"] <- "std.all.ml"
  ml_subset$std.all.ml <- round(ml_subset$std.all.ml, 3)

  # Fusionner par la clé avec all.x = TRUE pour garder l'ordre de lavaan
  merged <- merge(lavaan_subset, ml_subset, by = "merge_key", all.x = TRUE, sort = FALSE)

  # Garder seulement les colonnes utiles
  merged <- merged[, c("lhs", "op", "rhs", "std.all.lavaan", "std.all.ml")]

  # Ajouter une colonne de différence
  merged$difference <- round(merged$std.all.ml - merged$std.all.lavaan, 3)

  return(merged)
}


merge_estimates <- function(ml_param, svd_param, lavaan_param) {
  # Créer des clés de fusion basées sur lhs, op, rhs
  ml_param$merge_key <- paste(ml_param$lhs, ml_param$op, ml_param$rhs, sep = "___")
  svd_param$merge_key <- paste(svd_param$lhs, svd_param$op, svd_param$rhs, sep = "___")
  lavaan_param$merge_key <- paste(lavaan_param$lhs, lavaan_param$op, lavaan_param$rhs, sep = "___")

  # Garder seulement les colonnes nécessaires de lavaan (référence)
  lavaan_subset <- lavaan_param[, c("merge_key", "lhs", "op", "rhs", "std.all")]
  colnames(lavaan_subset)[colnames(lavaan_subset) == "std.all"] <- "std.all.lavaan"
  lavaan_subset$std.all.lavaan <- round(lavaan_subset$std.all.lavaan, 3)

  # Garder seulement les colonnes nécessaires de ML
  ml_subset <- ml_param[, c("merge_key", "std.all")]
  colnames(ml_subset)[colnames(ml_subset) == "std.all"] <- "std.all.ml"
  ml_subset$std.all.ml <- round(ml_subset$std.all.ml, 3)

  # Garder seulement les colonnes nécessaires de SVD
  svd_subset <- svd_param[, c("merge_key", "std.all")]
  colnames(svd_subset)[colnames(svd_subset) == "std.all"] <- "std.all.svd"
  svd_subset$std.all.svd <- round(svd_subset$std.all.svd, 3)

  # Fusionner par la clé avec all.x = TRUE pour garder l'ordre de lavaan
  merged <- merge(lavaan_subset, ml_subset, by = "merge_key", all.x = TRUE, sort = FALSE)
  merged <- merge(merged, svd_subset, by = "merge_key", all.x = TRUE, sort = FALSE)

  # Garder seulement les colonnes utiles (sans la colonne différence)
  merged <- merged[, c("lhs", "op", "rhs", "std.all.lavaan", "std.all.ml", "std.all.svd")]

  return(merged)
}