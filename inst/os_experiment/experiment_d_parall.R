# ============================================================
# Monte Carlo : variation du nombre d'indicateurs par bloc
# ============================================================

rm(list = ls())

devtools::load_all()
source("inst/os_experiment/utils_os.R")
source("inst/os_experiment/population.R")

library(parallel)


# ============================================================
# 1. Configuration
# ============================================================

q_grid <- c(3, 5, 10, 15, 20)
J <- 6L
N <- 2000
n_rep <- 10
n_workers <- 10L

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)


# ============================================================
# 2. Préparer la population pour chaque q
# ============================================================

set.seed(123)

population_states <- lapply(q_grid, function(q_value) {

  # Environnement séparé pour ne pas mélanger les populations
  env <- new.env(parent = globalenv())

  env$q <- q_value
  env$J <- J
  env$N <- N

  sys.source(
    "inst/simulations/data_simulation_mixed.R",
    envir = env
  )

  stopifnot(
    exists("true_param_with_S", envir = env, inherits = FALSE)
  )

  # Sauvegarder les objets effectivement créés par le script
  as.list(env, all.names = TRUE)
})

names(population_states) <- as.character(q_grid)

theta_true_list <- lapply(
  population_states,
  function(x) x$true_param_with_S
)

lavaan_models <- lapply(
  q_grid,
  function(q_value) generate_sem_model_J(J, q_value)
)

names(lavaan_models) <- as.character(q_grid)


# ============================================================
# 3. Stockage des estimations
# ============================================================

estimates_svd <- lapply(theta_true_list, function(theta_true) {

  parameter_names <- names(theta_true)

  if (is.null(parameter_names)) {
    parameter_names <- paste0("theta_", seq_along(theta_true))
  }

  matrix(
    NA_real_,
    nrow = n_rep,
    ncol = length(theta_true),
    dimnames = list(
      replication = seq_len(n_rep),
      parameter = parameter_names
    )
  )
})

estimates_os <- estimates_svd
estimates_rml <- estimates_svd


# ============================================================
# 4. Tâches et graines
# ============================================================

tasks <- expand.grid(
  rep_id = seq_len(n_rep),
  q_index = seq_along(q_grid)
)

tasks$q <- q_grid[tasks$q_index]

set.seed(0)

tasks$seed <- sample.int(
  .Machine$integer.max,
  nrow(tasks),
  replace = FALSE
)


# ============================================================
# 5. Une réplication Monte Carlo
# ============================================================

run_replication_q <- function(task) {

  q <- task$q
  state <- population_states[[task$q_index]]
  d <- length(state$true_param_with_S)

  # Restaurer les objets population dans un environnement local
  sample_env <- list2env(
    state,
    parent = environment()
  )

  sample_env$q <- q
  sample_env$J <- J
  sample_env$N <- N

  set.seed(task$seed)

  sys.source(
    "inst/model/model_mixed.R",
    envir = sample_env
  )

  X_2 <- sample_env$X_2
  Y_2 <- sample_env$Y_2
  C <- sample_env$C
  mode <- sample_env$mode

  # ----------------------------------------------------------
  # Lavaan
  # ----------------------------------------------------------

  fit_lavaan <- safe_estimation_lavaan(
    model = lavaan_models[[task$q_index]],
    data = X_2,
    composites_cov = "fixed"
  )

  # ----------------------------------------------------------
  # SVD
  # ----------------------------------------------------------

  modelsvd <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "svd"
  )

  fit_svd <- safe_estimation(modelsvd, compute_cov = FALSE)

  # ----------------------------------------------------------
  # One-step
  # ----------------------------------------------------------

  modelos <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "one_step"
  )

  fit_os <- safe_estimation(modelos, compute_cov = FALSE)

  # ----------------------------------------------------------
  # RML
  # ----------------------------------------------------------

  modelml <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "ml"
  )

  fit_rml <- safe_estimation(modelml, compute_cov = FALSE)

  # ----------------------------------------------------------
  # Diagnostics dépendant des données
  # ----------------------------------------------------------

  S <- cov(X_2)
  model_spec <- modelsvd$get_model()

  objective_fun <- function(x, S) {
    F1(x, S, model_spec)
  }


  diagnostic_errors <- character()

  safe_diagnostic <- function(label, expr) {
    tryCatch(
      {
        value <- as.numeric(force(expr))

        if (length(value) != 1L || !is.finite(value)) {
          stop("Diagnostic non scalaire ou non fini")
        }

        value
      },
      error = function(e) {
        diagnostic_errors <<- c(
          diagnostic_errors,
          paste0(label, ": ", conditionMessage(e))
        )
        NA_real_
      }
    )
  }



  get_gap <- function(fit, label) {
    if (!isTRUE(fit$success) || !isTRUE(fit_rml$success)) {
      return(NA_real_)
    }

    safe_diagnostic(
      label,
      likelihood_gap(
        theta_os = fit$theta,
        theta_rml = fit_rml$theta,
        S = S,
        objective_fun = objective_fun
      )
    )
  }


  objective_gap <- get_gap(fit_os, "likelihood_gap")
  objective_gap_svd <- get_gap(fit_svd, "likelihood_gap_svd")

  # ----------------------------------------------------------
  # Extraire seulement les résultats utiles
  # ----------------------------------------------------------

  get_theta <- function(fit) {
    if (!isTRUE(fit$success)) {
      return(rep(NA_real_, d))
    }

    stopifnot(length(fit$theta) == d)
    as.numeric(fit$theta)
  }

  error_text <- function(x) {
    if (length(x) == 0L) NA_character_
    else paste(x, collapse = "; ")
  }

  list(
    q_index = task$q_index,
    rep_id = task$rep_id,

    theta_svd = get_theta(fit_svd),
    theta_os = get_theta(fit_os),
    theta_rml = get_theta(fit_rml),

    diagnostics = data.frame(
      J = J,
      q = q,
      p = ncol(X_2),
      d = d,
      n = N,
      replication = task$rep_id,
      seed = task$seed,

      svd_success = isTRUE(fit_svd$success),
      os_success = isTRUE(fit_os$success),
      rml_success = isTRUE(fit_rml$success),
      lavaan_success = isTRUE(fit_lavaan$success),


      likelihood_gap = objective_gap,
      likelihood_gap_svd = objective_gap_svd,

      svd_time = fit_svd$elapsed_time,
      os_time = fit_os$elapsed_time,
      rml_time = fit_rml$elapsed_time,
      lavaan_time = fit_lavaan$elapsed_time,

      svd_error = error_text(fit_svd$error_message),
      os_error = error_text(fit_os$error_message),
      rml_error = error_text(fit_rml$error_message),
      lavaan_error = error_text(fit_lavaan$error_message),
      diagnostic_error = error_text(diagnostic_errors),

      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# 6. Exécution parallèle
# ============================================================

run_parallel_q <- function(tasks, population_states,
                           lavaan_models, N, J,
                           n_workers, project_dir) {

  cl <- makePSOCKcluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(
    cl,
    varlist = c(
      "project_dir",
      "population_states",
      "lavaan_models",
      "N",
      "J"
    ),
    envir = environment()
  )

  clusterEvalQ(cl, {

    setwd(project_dir)

    devtools::load_all(quiet = TRUE)
    library(lavaan)

    source("inst/os_experiment/utils_os.R")
    source("inst/os_experiment/population.R")



    NULL
  })

  task_list <- lapply(
    seq_len(nrow(tasks)),
    function(k) tasks[k, , drop = FALSE]
  )

  parLapplyLB(
    cl,
    task_list,
    fun = run_replication_q
  )
}

total_time <- system.time({

  outputs <- run_parallel_q(
    tasks = tasks,
    population_states = population_states,
    lavaan_models = lavaan_models,
    N = N,
    J = J,
    n_workers = n_workers,
    project_dir = project_dir
  )

})




# ============================================================
# 7. Assembler les résultats
# ============================================================

results <- vector("list", length(outputs))

for (k in seq_along(outputs)) {

  out <- outputs[[k]]

  i <- out$q_index
  r <- out$rep_id

  estimates_svd[[i]][r, ] <- out$theta_svd
  estimates_os[[i]][r, ] <- out$theta_os
  estimates_rml[[i]][r, ] <- out$theta_rml

  results[[k]] <- out$diagnostics
}

results_mc <- do.call(rbind, results)
rownames(results_mc) <- NULL

rm(outputs, results)


# ============================================================
# 8. Vérifier les taux de réussite
# ============================================================

print(
  aggregate(
    cbind(
      svd_success,
      os_success,
      rml_success,
      lavaan_success
    ) ~ q,
    data = results_mc,
    FUN = mean
  )
)


# ============================================================
# 9. Sauvegarder
# ============================================================

# saveRDS(
#   list(
#     results_mc = results_mc,
#
#     estimates_svd = estimates_svd,
#     estimates_os = estimates_os,
#     estimates_rml = estimates_rml,
#
#     theta_true_list = theta_true_list,
#     population_states = population_states,
#     lavaan_models = lavaan_models,
#
#     q_grid = q_grid,
#     J = J,
#     N = N,
#     n_rep = n_rep,
#     tasks = tasks,
#
#     n_workers = n_workers,
#     total_time = total_time,
#     session_info = sessionInfo()
#   ),
#   file = "inst/os_experiment/experiment_q_parallelized.rds"
# )


run_replication_q(tasks[1, , drop = FALSE])