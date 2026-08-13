## ============================================================
## COMPUTE SUMMARY FOR THE PAPER
## ============================================================

files <- c(
  n = "inst/os_experiment/monte_carlo_scaling_n_full.rds",
  J = "inst/os_experiment/monte_carlo_scaling_J_full.rds",
  q = "inst/os_experiment/monte_carlo_scaling_q_full.rds"
)

## ------------------------------------------------------------
## Helper: extract results_mc if the RDS is a list
## ------------------------------------------------------------

load_results <- function(path) {
  x <- readRDS(path)

  if (is.data.frame(x)) {
    return(x)
  }

  if (is.list(x) && "results_mc" %in% names(x)) {
    return(x$results_mc)
  }

  stop(
    "Cannot identify results_mc in: ", path,
    "\nAvailable names: ", paste(names(x), collapse = ", ")
  )
}


results <- lapply(files, load_results)


## ------------------------------------------------------------
## 1. Number of Monte Carlo replications
## ------------------------------------------------------------

cat("\n========================================\n")
cat("MONTE CARLO REPLICATIONS\n")
cat("========================================\n")

replications <- sapply(results, nrow)

print(replications)

cat("\nTotal Monte Carlo replications:",
    sum(replications), "\n")


## ------------------------------------------------------------
## 2. Identify runtime columns
## ------------------------------------------------------------

cat("\n========================================\n")
cat("RUNTIME COLUMNS\n")
cat("========================================\n")

for (exp_name in names(results)) {

  df <- results[[exp_name]]

  time_cols <- grep(
    "time|elapsed",
    names(df),
    ignore.case = TRUE,
    value = TRUE
  )

  cat("\nExperiment:", exp_name, "\n")
  print(time_cols)
}


## ------------------------------------------------------------
## 3. Total estimator runtime
## ------------------------------------------------------------

runtime_summary <- list()

for (exp_name in names(results)) {

  df <- results[[exp_name]]

  time_cols <- grep(
    "time|elapsed",
    names(df),
    ignore.case = TRUE,
    value = TRUE
  )

  ## Keep numeric runtime columns only
  time_cols <- time_cols[
    sapply(df[time_cols], is.numeric)
  ]

  if (length(time_cols) == 0) next

  tmp <- data.frame(
    experiment = exp_name,
    method = time_cols,
    n_fits = sapply(df[time_cols], function(x) sum(is.finite(x))),
    total_seconds = sapply(
      df[time_cols],
      function(x) sum(x[is.finite(x)], na.rm = TRUE)
    ),
    median_seconds = sapply(
      df[time_cols],
      function(x) median(x[is.finite(x)], na.rm = TRUE)
    ),
    mean_seconds = sapply(
      df[time_cols],
      function(x) mean(x[is.finite(x)], na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  runtime_summary[[exp_name]] <- tmp
}

runtime_summary <- do.call(rbind, runtime_summary)

runtime_summary$total_minutes <-
  runtime_summary$total_seconds / 60

runtime_summary$total_hours <-
  runtime_summary$total_seconds / 3600

rownames(runtime_summary) <- NULL

cat("\n========================================\n")
cat("RUNTIME SUMMARY BY EXPERIMENT / METHOD\n")
cat("========================================\n")

print(runtime_summary)


## ------------------------------------------------------------
## 4. Grand totals
## ------------------------------------------------------------

cat("\n========================================\n")
cat("GRAND TOTALS\n")
cat("========================================\n")

total_fits <- sum(runtime_summary$n_fits)

total_seconds <- sum(runtime_summary$total_seconds)

cat("Total estimator fits:", total_fits, "\n")
cat("Total accumulated estimator time:",
    round(total_seconds, 2), "seconds\n")
cat("Total accumulated estimator time:",
    round(total_seconds / 60, 2), "minutes\n")
cat("Total accumulated estimator time:",
    round(total_seconds / 3600, 2), "hours\n")


## ------------------------------------------------------------
## 5. Total by experiment
## ------------------------------------------------------------

time_by_experiment <- aggregate(
  cbind(n_fits, total_seconds) ~ experiment,
  data = runtime_summary,
  FUN = sum
)

time_by_experiment$total_hours <-
  time_by_experiment$total_seconds / 3600

cat("\n========================================\n")
cat("TOTAL BY EXPERIMENT\n")
cat("========================================\n")

print(time_by_experiment)


## ------------------------------------------------------------
## 6. Total by method/column
## ------------------------------------------------------------

time_by_method <- aggregate(
  cbind(n_fits, total_seconds) ~ method,
  data = runtime_summary,
  FUN = sum
)

time_by_method$total_hours <-
  time_by_method$total_seconds / 3600

cat("\n========================================\n")
cat("TOTAL BY METHOD\n")
cat("========================================\n")

print(time_by_method)


## ------------------------------------------------------------
## 7. Available computing resources
## ------------------------------------------------------------

cat("\n========================================\n")
cat("COMPUTING RESOURCES\n")
cat("========================================\n")

cat("R:", R.version.string, "\n")
cat("Logical cores detected:",
    parallel::detectCores(logical = TRUE), "\n")
cat("Physical cores detected:",
    parallel::detectCores(logical = FALSE), "\n")

system(
  'powershell -command "Get-CimInstance Win32_Processor | Select-Object Name,NumberOfCores,NumberOfLogicalProcessors"'
)

system(
  'powershell -command "Get-CimInstance Win32_ComputerSystem | Select-Object @{Name=\'RAM_GB\';Expression={[math]::Round($_.TotalPhysicalMemory/1GB,2)}}"'
)

cat("nloptr:", as.character(packageVersion("nloptr")), "\n")
cat("lavaan:", as.character(packageVersion("lavaan")), "\n")