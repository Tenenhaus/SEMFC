# ============================================================
# RMSE par J : estimation, MCSE et intervalle Monte Carlo
# ============================================================

table_rmse_J <- function(estimates, theta_true_list, J_grid,
                         level = 0.95) {

  stopifnot(
    length(estimates) == length(J_grid),
    length(theta_true_list) == length(J_grid)
  )

  tables <- lapply(seq_along(J_grid), function(i) {

    theta <- as.matrix(estimates[[i]])
    theta_true <- as.numeric(theta_true_list[[i]])

    stopifnot(ncol(theta) == length(theta_true))

    errors <- sweep(theta, 2, theta_true, "-")

    # MSE par réplication, moyenne sur les paramètres
    mse <- rowMeans(errors^2)

    # Exclure les vecteurs incomplets ou non finis
    mse[rowSums(!is.finite(errors)) > 0] <- NA_real_

    tab <- summarise_mc(
      x = matrix(mse, nrow = 1),
      n_grid = J_grid[i],
      level = level,
      root = TRUE
    )

    names(tab)[names(tab) == "n"] <- "J"
    tab$d <- ncol(theta)

    tab[, c(
      "J", "d", "n_valid", "n_invalid",
      "estimate", "mcse", "lower", "upper"
    )]
  })

  do.call(rbind, tables)
}


# ============================================================
# Temps par J : médiane et quartiles
# ============================================================

table_time_J <- function(results_mc, column) {

  stopifnot(all(c("J", column) %in% names(results_mc)))

  tables <- lapply(sort(unique(results_mc$J)), function(J_value) {

    times <- results_mc[[column]][results_mc$J == J_value]

    valid <- is.finite(times) & times >= 0
    x <- times[valid]

    quartiles <- if (length(x) > 0L) {
      quantile(
        x,
        probs = c(0.25, 0.50, 0.75),
        names = FALSE
      )
    } else {
      rep(NA_real_, 3)
    }

    data.frame(
      J = J_value,
      n_valid = length(x),
      n_invalid = sum(!valid),
      median = quartiles[2],
      q25 = quartiles[1],
      q75 = quartiles[3]
    )
  })

  do.call(rbind, tables)
}