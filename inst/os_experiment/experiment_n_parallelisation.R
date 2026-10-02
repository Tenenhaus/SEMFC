library(parallel)


devtools::load_all()
rm(list = ls())
source('inst/os_experiment/utils_os.R')

set.seed(123)
q <- 3
source('inst/simulations/data_simulation_mixed.R')

# ============================================================
# 1. Simulation parameters
# ============================================================

n_grid <- c(300, 600, 1200, 2400, 4800)
n_rep  <- 10

# True parameter vector.
# It must follow exactly the same ordering as the estimators.
theta_true <- true_param_with_S

d <- length(theta_true)


# ============================================================
# 4. Storage objects
# ============================================================

results <- vector(
  mode = "list",
  length = length(n_grid) * n_rep
)

# Optional storage of all parameter estimates.
# Useful for componentwise bias and coverage analyses.
estimates_svd <- array(
  NA_real_,
  dim = c(length(n_grid), n_rep, d),
  dimnames = list(
    n = as.character(n_grid),
    replication = seq_len(n_rep),
    parameter = names(theta_true)
  )
)

estimates_os <- estimates_svd
estimates_rml <- estimates_svd


Vhat_os <- array(
  NA_real_,
  dim = c(length(n_grid), n_rep, d, d),
  dimnames = list(
    n = as.character(n_grid),
    replication = seq_len(n_rep),
    row_parameter = names(theta_true),
    col_parameter = names(theta_true)
  )
)

Vhat_rml <- Vhat_os


variance_mc_error_os  <- numeric(length(n_grid))
variance_mc_error_rml <- numeric(length(n_grid))

V_mc_os <- array(
  NA_real_,
  dim = c(length(n_grid), d, d),
  dimnames = list(
    n = as.character(n_grid),
    row_parameter = names(theta_true),
    col_parameter = names(theta_true)
  )
)
V_mc_rml <- V_mc_os
norm_V0 <- norm(V_0, type = "F")

# ============================================================
# 1. Configuration
# ============================================================

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)


n_workers <- 10L

# Une ligne = une tâche indépendante
tasks <- expand.grid(
  rep_id = seq_len(n_rep),
  n_index = seq_along(n_grid)
)

tasks$N <- n_grid[tasks$n_index]


set.seed(0)
tasks$seed <- sample.int(
  .Machine$integer.max,
  nrow(tasks),
  replace = FALSE
)


# ============================================================
# 2. Une tâche : simulation + trois estimations
# ============================================================

run_replication <- function(task) {

  N <- task$N
  set.seed(task$seed)

  # local = TRUE est important :
  # les données générées appartiennent à cette réplication
  source("inst/model/model_mixed.R", local = TRUE)

  modelsvd <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "svd"
  )

  fit_svd <- safe_estimation(modelsvd)

  modelos <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "one_step"
  )

  fit_os <- safe_estimation(modelos)

  modelml <- SemFC$new(
    Y_2,
    relation_matrix = C,
    mode = mode,
    estimator = "ml"
  )

  fit_rml <- safe_estimation(modelml)

  # Covariance calculée une seule fois
  S <- cov(X_2)

  objective_fun <- function(x, S) {
    F1(x, S, modelsvd$get_model())
  }

  objective_gap <- if (
    isTRUE(fit_os$success) && isTRUE(fit_rml$success)
  ) {
    likelihood_gap(
      theta_os = fit_os$theta,
      theta_rml = fit_rml$theta,
      S = S,
      objective_fun = objective_fun
    )
  } else {
    NA_real_
  }

  objective_gap_svd <- if (
    isTRUE(fit_svd$success) && isTRUE(fit_rml$success)
  ) {
    likelihood_gap(
      theta_os = fit_svd$theta,
      theta_rml = fit_rml$theta,
      S = S,
      objective_fun = objective_fun
    )
  } else {
    NA_real_
  }

  # Garantir une chaîne même si error_message est NULL
  error_text <- function(x) {
    if (length(x) == 0L) NA_character_
    else paste(x, collapse = "; ")
  }

  # Ne retourner que les informations utiles
  list(
    n_index = task$n_index,
    rep_id = task$rep_id,

    theta_svd = if (isTRUE(fit_svd$success)) {
      as.numeric(fit_svd$theta)
    } else {
      rep(NA_real_, d)
    },

    theta_os = if (isTRUE(fit_os$success)) {
      as.numeric(fit_os$theta)
    } else {
      rep(NA_real_, d)
    },

    theta_rml = if (isTRUE(fit_rml$success)) {
      as.numeric(fit_rml$theta)
    } else {
      rep(NA_real_, d)
    },

    V_os = if (isTRUE(fit_os$success)) {
      as.matrix(fit_os$V)
    } else {
      matrix(NA_real_, d, d)
    },

    V_rml = if (isTRUE(fit_rml$success)) {
      as.matrix(fit_rml$V)
    } else {
      matrix(NA_real_, d, d)
    },

    diagnostics = data.frame(
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

      stringsAsFactors = FALSE
    )
  )
}


# ============================================================
# 3. Initialiser les processus et lancer les tâches
# ============================================================

run_parallel_mc <- function(tasks, n_workers, project_dir, q, d) {

  cl <- makePSOCKcluster(n_workers)

  # Fermer les processus même si une erreur survient
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(
    cl,
    varlist = c("project_dir", "q", "d"),
    envir = environment()
  )

  # Exécuté une seule fois dans chaque processus
  clusterEvalQ(cl, {

    setwd(project_dir)

    devtools::load_all(quiet = TRUE)
    source("inst/os_experiment/utils_os.R")

    # Recréer le même modèle population que dans ton script
    set.seed(123)
    source("inst/simulations/data_simulation_mixed.R")

    NULL
  })

  task_list <- lapply(
    seq_len(nrow(tasks)),
    function(k) tasks[k, , drop = FALSE]
  )

  # Répartir les tâches au fur et à mesure que les processus se libèrent
  parLapplyLB(
    cl,
    task_list,
    fun = run_replication
  )
}

# outputs <- run_parallel_mc(
#   tasks = tasks,
#   n_workers = n_workers,
#   project_dir = project_dir,
#   q = q,
#   d = d
# )


total_time <- system.time({
  outputs <- run_parallel_mc(
    tasks = tasks,
    n_workers = n_workers,
    project_dir = project_dir,
    q = q,
    d = d
  )
})




# ============================================================
# 4. Assembler les résultats dans les tableaux existants
# ============================================================

for (k in seq_along(outputs)) {

  out <- outputs[[k]]

  i <- out$n_index
  r <- out$rep_id

  estimates_svd[i, r, ] <- out$theta_svd
  estimates_os[i, r, ]  <- out$theta_os
  estimates_rml[i, r, ] <- out$theta_rml

  Vhat_os[i, r, , ]  <- out$V_os
  Vhat_rml[i, r, , ] <- out$V_rml

  results[[k]] <- out$diagnostics
}

results_mc <- do.call(rbind, results)
rownames(results_mc) <- NULL

# Libérer la copie intermédiaire des résultats
rm(outputs)

# ============================================================
# 8. Export
# ============================================================

#
saveRDS(
  list(
    results_mc = results_mc,


    estimates_svd = estimates_svd,
    estimates_os = estimates_os,
    estimates_rml = estimates_rml,

    V_0 = V_0,
    Vhat_os = Vhat_os,
    Vhat_rml = Vhat_rml,


    theta_true = theta_true,
    n_grid = n_grid,
    n_rep = n_rep,

    tasks = tasks,
    total_time = total_time,
    session_info = sessionInfo()
  ),
  file = "inst/os_experiment/experiment_n_parallelized.rds"
)