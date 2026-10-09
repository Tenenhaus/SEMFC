devtools::load_all()


files <- list.files(path = "data", pattern = "\\.rda$", full.names = TRUE)
lapply(files, load, .GlobalEnv)



config_files <- list.files(
  path = "inst/os_experiment/dataset_configs",
  pattern = "\\.R$",
  full.names = TRUE,
  ignore.case = TRUE
)


invisible(lapply(config_files, source))


datasets <- list(
  list(
    name = "BergamiBagozzi2000",
    data = BergamiBagozzi2000,
    config_fun = create_bergami_config
  ),
  list(
    name = "ECSI",
    data = ECSI,
    config_fun = create_ecsi_full_config
  ),
  list(
    name = "ITFlex",
    data = ITFlex,
    config_fun = create_itflex_config
  ),
  list(
    name = "LancelotMiltgenetal2016",
    data = LancelotMiltgenetal2016,
    config_fun = create_lancelot_config
  ),
  list(
    name = "PoliticalDemocracy",
    data = PoliticalDemocracy,
    config_fun = create_politicaldemocracy_config
  ),
  list(
    name = "Russett",
    data = Russett,
    config_fun = create_russett_config
  )
)

comparison_real_data  <- function(config, dataset_name) {

  results <- lapply(c("svd", "one_step", "ml"), function(estimator) {

    model <- SemFC$new(
      config$data,
      relation_matrix = config$relation_matrix,
      mode = config$mode,
      estimator = estimator
    )

    fit <- safe_estimation(model)
    fit$model$summary(standardized = T)

    data.frame(
      dataset = dataset_name,
      estimator = switch(estimator,
                         svd = "SVD",
                         one_step = "OS",
                         ml = "RML"),
      f = F1(
        fit$model$get_estimate("theta"),
        fit$model$get_data()$cov_S,
        fit$model$get_model()
      ),
      time = fit$elapsed_time
    )
  })

  fit <- safe_estimation_lavaan(
    model = config$sem,
    data = data.frame(Reduce("cbind", config$data)),
    composites_cov = "free"
  )
  show(summary(fit$fit, standardized = TRUE))

  results[[4]] <- data.frame(
    dataset = dataset_name,
    estimator = "lavaan",
    f = 2 * lavaan::lavInspect(fit$fit, "optim")$fx,
    time = fit$elapsed_time
  )

  do.call(rbind, results)
}


set.seed(0)
results_all <- do.call(rbind, lapply(datasets, function(x) {
  message("Dataset: ", x$name)

  config <- x$config_fun(x$data)
  comparison_real_data(config, x$name)
}))