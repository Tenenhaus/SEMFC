# ============================================================
# Monte Carlo : variation du nombre de blocs
# ============================================================

rm(list = ls())

devtools::load_all()
source("inst/os_experiment/utils_os.R")
source("inst/os_experiment/population.R")

library(parallel)


# ============================================================
# 1. Configuration
# ============================================================

J_grid <- c(6, 10, 14, 18, 22)
N <- 2000
q <- 3
n_rep <- 10
n_workers <- 10L

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)


# ============================================================
# 2. Populations : une seule construction par J
# ============================================================

set.seed(123)

populations <- lapply(
  J_grid,
  function(J) make_population_J(J, q)
)

names(populations) <- as.character(J_grid)

theta_true_list <- lapply(
  populations,
  function(population) population$true_param_with_S
)


# ============================================================
# 3. Stockage : la dimension de theta varie avec J
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
  J_index = seq_along(J_grid)
)

tasks$J <- J_grid[tasks$J_index]

set.seed(0)

tasks$seed <- sample.int(
  .Machine$integer.max,
  nrow(tasks),
  replace = FALSE
)


# ============================================================
# 5. Une réplication
# ============================================================

run_replication_J <- function(task) {

  J <- task$J
  population <- populations[[task$J_index]]
  d <- length(population$true_param_with_S)

  set.seed(task$seed)

  sim <- generate_mixed_sample(
    N = N,
    J = J,
    q = q,
    SIGMA = population$SIGMA,
    empirical = FALSE
  )

  X_2 <- sim$X
  Y_2 <- sim$Y
  C <- sim$C
  mode <- sim$mode

  # ----------------------------------------------------------
  # Estimations sur le même échantillon
  # ----------------------------------------------------------

  modelsvd <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "svd"
  )

  fit_svd <- safe_estimation(modelsvd, compute_cov = F)

  modelos <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "one_step"
  )

  fit_os <- safe_estimation(modelos, compute_cov = F)

  modelml <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "ml"
  )

  fit_rml <- safe_estimation(modelml, compute_cov = F)

  # ----------------------------------------------------------
  # Diagnostics dépendant de l'échantillon
  # ----------------------------------------------------------

  S <- cov(X_2)
  model_spec <- modelsvd$get_model()

  objective_fun <- function(x, S) {
    F1(x, S, model_spec)
  }

  constraint_fun <- function(x) {
    heq(x, S, model_spec)
  }

  # Une erreur de diagnostic ne doit pas effacer les estimations
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
  # Extraction des résultats utiles
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
    J_index = task$J_index,
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



      likelihood_gap = objective_gap,
      likelihood_gap_svd = objective_gap_svd,

      svd_time = fit_svd$elapsed_time,
      os_time = fit_os$elapsed_time,
      rml_time = fit_rml$elapsed_time,

      svd_error = error_text(fit_svd$error_message),
      os_error = error_text(fit_os$error_message),
      rml_error = error_text(fit_rml$error_message),
      diagnostic_error = error_text(diagnostic_errors),

      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# 6. Exécution parallèle
# ============================================================

run_parallel_J <- function(tasks, populations, N, q,
                           n_workers, project_dir) {

  cl <- makePSOCKcluster(n_workers)
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(
    cl,
    varlist = c("project_dir", "populations", "N", "q"),
    envir = environment()
  )

  clusterEvalQ(cl, {

    setwd(project_dir)

    devtools::load_all(quiet = TRUE)
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
    fun = run_replication_J
  )
}

total_time <- system.time({

  outputs <- run_parallel_J(
    tasks = tasks,
    populations = populations,
    N = N,
    q = q,
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

  i <- out$J_index
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
    cbind(svd_success, os_success, rml_success) ~ J,
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
#     populations = populations,
#
#     J_grid = J_grid,
#     N = N,
#     q = q,
#     n_rep = n_rep,
#     tasks = tasks,
#
#     n_workers = n_workers,
#     total_time = total_time,
#     session_info = sessionInfo()
#   ),
#   file = "inst/os_experiment/experiment_J_parallelized.rds"
# )