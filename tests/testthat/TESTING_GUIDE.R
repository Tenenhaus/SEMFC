# ============================================================================
# EXEMPLE D'UTILISATION DU SYSTÈME DE CONFIGURATIONS
# ============================================================================

# Localisation correcte du fichier testthat
# Si vous exécutez ce script depuis la racine du projet:

# 1️⃣ CHARGER LES CONFIGURATIONS
# ================================================
source("tests/testthat/dataset_configs/load_configs.R")

# 2️⃣ CRÉER LES RÉFÉRENCES POUR LES TESTS
# ================================================
source("tests/testthat/setup_references.R")

# Créer les références pour SVD et ML sur tous les datasets:
create_all_references(estimators = c("svd", "ml"))

# Ou seulement SVD:
create_all_references(estimators = "svd")

# Ou seulement ML:
create_all_references(estimators = "ml")

# 3️⃣ UTILISER LES CONFIGURATIONS DANS VOS TESTS
# ================================================
# Dans vos tests, vous pouvez faire:

library(testthat)

test_that("parameterEstimates retourne la même chose", {
  # Charger la config ECSI
  setup <- setup_test_data_for_reference("ECSI")

  # Créer et fitter le modèle
  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    estimator = "svd"
  )
  sem_model$fit()

  # Récupérer les estimations
  actual <- sem_model$parameterEstimates(standardized = FALSE)

  # Charger la référence
  expected <- readRDS("tests/testthat/fixtures/parameterEstimates_reference_svd_ECSI.rds")

  # Comparer
  expect_equal(actual, expected)
})

# 4️⃣ TESTER PLUSIEURS DATASETS
# ================================================
# Parcourir tous les datasets automatiquement:

datasets <- c("ECSI", "BergamiBagozzi2000", "ITFlex",
              "PoliticalDemocracy", "Russett", "LancelotMiltgenetal2016")

for (dataset_name in datasets) {
  test_that(paste("parameterEstimates fonctionne pour", dataset_name), {
    setup <- setup_test_data_for_reference(dataset_name)

    sem_model <- SemFC$new(
      data = setup$data,
      relation_matrix = setup$relation_matrix,
      estimator = "svd"
    )
    sem_model$fit()

    actual <- sem_model$parameterEstimates(standardized = FALSE)
    expected <- readRDS(paste0("tests/testthat/fixtures/parameterEstimates_reference_svd_",
                                dataset_name, ".rds"))

    expect_equal(actual, expected)
  })
}

# 5️⃣ STRUCTURE COMPLÈTE D'UN TEST FICHIER
# ================================================
# tests/testthat/test_parameterEstimates.R

test_that("parameterEstimates() retourne une table valide - ECSI SVD", {
  setup <- setup_test_data_for_reference("ECSI")
  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    estimator = "svd"
  )
  sem_model$fit()

  actual <- sem_model$parameterEstimates(standardized = FALSE)
  expected <- readRDS("tests/testthat/fixtures/parameterEstimates_reference_svd_ECSI.rds")

  expect_equal(actual, expected)
  expect_is(actual, "data.frame")
  expect_true(nrow(actual) > 0)
})

test_that("parameterEstimates() retourne une table valide - ECSI ML", {
  setup <- setup_test_data_for_reference("ECSI")
  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    estimator = "ml"
  )
  sem_model$fit()

  actual <- sem_model$parameterEstimates(standardized = FALSE)
  expected <- readRDS("tests/testthat/fixtures/parameterEstimates_reference_ml_ECSI.rds")

  expect_equal(actual, expected)
  expect_is(actual, "data.frame")
})

# ============================================================================
# RÉSUMÉ DES 3 MÉTHODES DE TEST
# ============================================================================

# 1️⃣ SNAPSHOT TESTING (Recommandé pour les données complexes)
# ============================================================
# Comparer directement avec une version enregistrée
test_that("resultats stables", {
  setup <- setup_test_data_for_reference("ECSI")
  sem_model <- SemFC$new(data = setup$data, relation_matrix = setup$relation_matrix)
  sem_model$fit()

  result <- sem_model$parameterEstimates()
  expected <- readRDS("tests/testthat/fixtures/parameterEstimates_reference_svd_ECSI.rds")

  expect_equal(result, expected)
})

# 2️⃣ VALIDATION DES PROPRIÉTÉS (Complémentaire)
# =============================================
# Vérifier les propriétés sans stocker la référence
test_that("resultats ont bonnes proprietes", {
  setup <- setup_test_data_for_reference("ECSI")
  sem_model <- SemFC$new(data = setup$data, relation_matrix = setup$relation_matrix)
  sem_model$fit()

  result <- sem_model$parameterEstimates()

  # Vérifications
  expect_is(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(c("estimate", "se", "t_value") %in% colnames(result)))
})

# 3️⃣ TESTS PARAMÉTRÉS MULTI-DATASETS
# ===================================
# Exécuter le même test sur plusieurs datasets
parameterized_test <- function(dataset_name, estimator = "svd") {
  setup <- setup_test_data_for_reference(dataset_name)
  sem_model <- SemFC$new(
    data = setup$data,
    relation_matrix = setup$relation_matrix,
    estimator = estimator
  )
  sem_model$fit()

  result <- sem_model$parameterEstimates()
  ref_file <- paste0("tests/testthat/fixtures/parameterEstimates_reference_",
                     estimator, "_", dataset_name, ".rds")
  expected <- readRDS(ref_file)

  expect_equal(result, expected)
}

# Utiliser la fonction paramétrée:
for (ds in c("ECSI", "BergamiBagozzi2000")) {
  test_that(paste(ds, "- SVD"), {
    parameterized_test(ds, "svd")
  })

  test_that(paste(ds, "- ML"), {
    parameterized_test(ds, "ml")
  })
}

