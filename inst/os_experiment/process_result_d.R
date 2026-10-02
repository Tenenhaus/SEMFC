table_rmse_q <- function(estimates, theta_true_list, q_grid,
                         level = 0.95) {

  tab <- table_rmse_J(
    estimates = estimates,
    theta_true_list = theta_true_list,
    J_grid = q_grid,
    level = level
  )

  names(tab)[names(tab) == "J"] <- "q"
  tab
}



table_time_q <- function(results_mc, column) {


  dat <- results_mc
  dat$J <- dat$q

  tab <- table_time_J(dat, column)
  names(tab)[names(tab) == "J"] <- "q"

  tab
}




rmse_os <- table_rmse_q(estimates_os, theta_true_list, q_grid)
rmse_rml <- table_rmse_q(estimates_rml, theta_true_list, q_grid)
rmse_svd <- table_rmse_q(estimates_svd, theta_true_list, q_grid)


time_os <- table_time_q(results_mc, "os_time")
time_rml <- table_time_q(results_mc, "rml_time")
time_svd <- table_time_q(results_mc, "svd_time")
time_lavaan <- table_time_q(results_mc, "lavaan_time")
