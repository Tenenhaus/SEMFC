# Guide des Tests pour le Package SEMFC

## Structure des Tests

Ce répertoire contient les tests unitaires pour le package SEMFC, utilisant le framework **testthat**.

### Fichiers de Test

#### 1. `test_parameterEstimates.R`
Tests spécifiques pour la méthode `parameterEstimates()` avec 3 approches :

- **Tests Détaillés** : Vérifications précises de la structure et du contenu
- **Regression Testing** : Comparaison avec une référence RDS stockée

#### 2. `test_semfc_methods.R`
Tests pour toutes les méthodes principales de la classe SemFC :

- `initialize()` - Création d'objets SemFC
- `fit()` - Ajustement du modèle (SVD, ML, avec/sans inference)
- `get_estimate()` - Récupération des estimations
- `check_improper()` - Détection de solutions impropres
- `summary()` - Rapport de résumé
- **Tests d'intégration** - Workflows complets

#### 3. `setup_references.R`
Script pour générer les fichiers de référence (fixtures) utilisés par les tests de regression.

### Dossier `fixtures/`
Contient les données de référence (fichiers RDS) utilisées pour comparer les résultats.

## Comment Exécuter les Tests

### 1. Depuis RStudio

```r
# Exécuter tous les tests
devtools::test()

# Ou avec testthat
testthat::test_dir("tests/testthat")

# Exécuter un fichier de test spécifique
testthat::test_file("tests/testthat/test_parameterEstimates.R")
```



## Configuration Initiale des Références

La première fois, vous devez générer les fichiers de référence RDS :

```r
# Charger le script
source("tests/testthat/setup_references.R")

# Générer toutes les références
create_all_references()

# Ou créer des références spécifiques
create_parameterEstimates_reference_svd()
create_parameterEstimates_reference_ml()
```

Cela créera les fichiers :
- `tests/testthat/fixtures/parameterEstimates_reference_svd.rds`
- `tests/testthat/fixtures/parameterEstimates_reference_ml.rds`



## Mettre à Jour les Références

Si vous changez intentionnellement le comportement de vos fonctions :

```r

# Regénérer les RDS
source("tests/testthat/setup_references.R")
create_all_references()
```

## Bonnes Pratiques pour Écrire des Tests

### 1. Utiliser une fixture commune
```r
setup_test_data <- function() {
  # ... configuration réutilisable ...
}

test_that("mon test", {
  setup <- setup_test_data()
  # ... test ...
})
```

### 2. Tester un seul concept par test
```r
# ❌ Mauvais - teste deux choses
test_that("fit and get_estimate work", { ... })

# ✅ Bon - teste une seule chose
test_that("fit() executes without error", { ... })
test_that("get_estimate() returns correct type", { ... })
```

### 3. Utiliser des noms clairs
```r
# ✅ Bon
test_that("parameterEstimates() returns data.frame with correct structure", { ... })

# ❌ Confus
test_that("it works", { ... })
```

### 4. Tester les cas limites
```r
# Modèle normal
test_that("get_estimate() works with fitted model", { ... })

# Edge case
test_that("get_estimate() fails before model is fitted", { ... })
```

## Interprétation des Résultats

Quand vous exécutez les tests, vous verrez :

```
test_semfc_methods.R:5: ✓ SemFC$new() creates valid object
test_semfc_methods.R:12: ✓ SemFC$new() initializes with different estimators
test_semfc_methods.R:10: ✗ fit() with bootstrap inference fails
  Error: ...
  In test_semfc_methods.R:25
    ...
```

- ✓ = **PASS** : Le test a réussi
- ✗ = **FAIL** : Le test a échoué
- ⚠ = **SKIP** : Le test a été ignoré

## Debugging

Pour déboguer un test qui échoue :

```r
# Activer le mode verbose
options(error = browser)

# Ou utiliser testthat interactivement
library(testthat)
with_reporter(DebugReporter$new(), {
  test_file("tests/testthat/test_parameterEstimates.R")
})
```

## Couverture de Code

Pour mesurer la couverture de code :

```r
# Installer covr
install.packages("covr")

# Générer un rapport
covr::report(covr::package_coverage())
```

## Contact et Support

Pour des questions sur les tests, consultez :
- [Documentation testthat](https://testthat.r-lib.org/)
- [Guide R Packages](https://r-pkgs.org/testing-design.html)

