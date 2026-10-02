# ============================================================
# 1. Fonction commune : résumer une métrique par taille n
# x : matrice [n, réplication]
#
# root = TRUE : transforme une moyenne de MSE en RMSE
#               avec MCSE obtenue par méthode delta.
# ============================================================

summarise_mc <- function(x, n_grid, level = 0.95,
                         root = FALSE, nonnegative = FALSE) {

  stopifnot(is.matrix(x), nrow(x) == length(n_grid))

  z <- qnorm((1 + level) / 2)

  do.call(rbind, lapply(seq_along(n_grid), function(i) {

    values <- x[i, ]
    values <- values[is.finite(values)]
    m <- length(values)

    avg <- if (m > 0) mean(values) else NA_real_
    se_avg <- if (m > 1) sd(values) / sqrt(m) else NA_real_

    estimate <- if (root) sqrt(avg) else avg

    mcse <- if (root) {
      if (is.finite(avg) && avg > 0) {
        se_avg / (2 * sqrt(avg))
      } else {
        NA_real_
      }
    } else {
      se_avg
    }

    lower <- estimate - z * mcse
    upper <- estimate + z * mcse

    if (root || nonnegative) {
      lower <- pmax(0, lower)
    }

    data.frame(
      n = n_grid[i],
      n_valid = m,
      n_invalid = ncol(x) - m,
      estimate = estimate,
      mcse = mcse,
      lower = lower,
      upper = upper
    )
  }))
}


# ============================================================
# 2. RMSE global
# sqrt(moyenne sur réplications et paramètres des erreurs²)
#
# Une réplication est exclue si son vecteur est incomplet.
# ============================================================

table_rmse <- function(estimates, theta_true, n_grid,
                       level = 0.95) {

  dims <- dim(estimates)
  stopifnot(length(dims) == 3L,
            dims[1] == length(n_grid),
            dims[3] == length(theta_true))

  errors <- matrix(
    sweep(estimates, 3, as.numeric(theta_true), "-"),
    nrow = dims[1] * dims[2],
    ncol = dims[3]
  )

  mse <- rowMeans(errors^2)
  mse[rowSums(!is.finite(errors)) > 0] <- NA_real_

  summarise_mc(
    matrix(mse, nrow = dims[1], ncol = dims[2]),
    n_grid,
    level = level,
    root = TRUE
  )
}


# ============================================================
# 3. Distance moyenne entre deux estimateurs
# scaled = TRUE : sqrt(n) * ||theta_a - theta_b||_2
# ============================================================

table_distance <- function(estimates_a, estimates_b, n_grid,
                           scaled = TRUE, level = 0.95) {

  stopifnot(identical(dim(estimates_a), dim(estimates_b)))

  dims <- dim(estimates_a)

  differences <- matrix(
    estimates_a - estimates_b,
    nrow = dims[1] * dims[2],
    ncol = dims[3]
  )

  distances <- sqrt(rowSums(differences^2))
  distances[rowSums(!is.finite(differences)) > 0] <- NA_real_

  distances <- matrix(
    distances, nrow = dims[1], ncol = dims[2]
  )

  if (scaled) {
    distances <- sweep(distances, 1, sqrt(n_grid), "*")
  }

  summarise_mc(
    distances, n_grid,
    level = level, nonnegative = TRUE
  )
}


# ============================================================
# 4. Erreur relative de covariance plug-in
# ||Vhat - V0||_F / ||V0||_F
# Vhat : [n, réplication, paramètre, paramètre]
# ============================================================

table_covariance_plugin <- function(Vhat, V0, n_grid,
                                    level = 0.95) {

  dims <- dim(Vhat)

  stopifnot(
    length(dims) == 4L,
    dims[1] == length(n_grid),
    all(dim(V0) == dims[3:4])
  )

  norm_V0 <- norm(as.matrix(V0), type = "F")
  stopifnot(is.finite(norm_V0), norm_V0 > 0)

  differences <- sweep(
    matrix(Vhat, nrow = dims[1] * dims[2]),
    2, as.vector(as.matrix(V0)), "-"
  )

  errors <- sqrt(rowSums(differences^2)) / norm_V0
  errors[rowSums(!is.finite(differences)) > 0] <- NA_real_

  summarise_mc(
    matrix(errors, nrow = dims[1], ncol = dims[2]),
    n_grid,
    level = level, nonnegative = TRUE
  )
}




# ============================================================
# 6. Moyenne d'une colonne de results_mc, avec intervalle MC
# Une ligne doit correspondre à une réplication.
# ============================================================

table_scalar <- function(results_mc, column,
                         level = 0.95, nonnegative = FALSE) {

  stopifnot(all(c("n", column) %in% names(results_mc)))

  ns <- sort(unique(results_mc$n))

  do.call(rbind, lapply(ns, function(n_value) {

    values <- results_mc[[column]][results_mc$n == n_value]

    summarise_mc(
      matrix(values, nrow = 1),
      n_grid = n_value,
      level = level,
      nonnegative = nonnegative
    )
  }))
}



table_covariance_mc <- function(estimates, V0, n_grid) {

  dims <- dim(estimates)
  stopifnot(length(dims) == 3L)

  d <- dims[3]
  V0 <- as.matrix(V0)

  stopifnot(
    dims[1] == length(n_grid),
    all(dim(V0) == c(d, d)),
    all(is.finite(V0))
  )

  norm_V0 <- norm(V0, type = "F")
  stopifnot(norm_V0 > 0)

  V_mc <- array(
    NA_real_,
    dim = c(length(n_grid), d, d),
    dimnames = list(
      n = as.character(n_grid),
      row_parameter = dimnames(estimates)[[3]],
      col_parameter = dimnames(estimates)[[3]]
    )
  )

  result <- data.frame(
    n = n_grid,
    n_valid = integer(length(n_grid)),
    n_invalid = integer(length(n_grid)),
    variance_mc_error = NA_real_
  )

  for (i in seq_along(n_grid)) {

    # Réplications × paramètres
    theta <- matrix(
      estimates[i, , , drop = FALSE],
      nrow = dims[2],
      ncol = d
    )

    # Exclure les vecteurs incomplets ou non finis
    valid <- rowSums(!is.finite(theta)) == 0
    theta <- theta[valid, , drop = FALSE]

    result$n_valid[i] <- nrow(theta)
    result$n_invalid[i] <- sum(!valid)

    if (nrow(theta) < 2L) next

    V_mc[i, , ] <- n_grid[i] * cov(theta)

    result$variance_mc_error[i] <-
      norm(V_mc[i, , ] - V0, type = "F") / norm_V0
  }

  list(table = result, V_mc = V_mc)
}


qqplot_parameter <- function(estimates, theta_true, V0, n_grid,
                             j = 1L, ncol = 3L) {

  dims <- dim(estimates)
  V0 <- as.matrix(V0)

  stopifnot(
    length(dims) == 3L,
    dims[1] == length(n_grid),
    length(theta_true) == dims[3],
    all(dim(V0) == c(dims[3], dims[3])),
    length(j) == 1L,
    j >= 1, j <= dims[3], j == as.integer(j),
    is.finite(V0[j, j]), V0[j, j] > 0
  )

  parameter_names <- dimnames(estimates)[[3]]
  parameter_label <- if (is.null(parameter_names)) {
    paste0("theta_", j)
  } else {
    parameter_names[j]
  }

  # Matrice n × réplications
  theta_j <- matrix(
    estimates[, , j, drop = FALSE],
    nrow = dims[1],
    ncol = dims[2]
  )

  Z <- sweep(
    theta_j - theta_true[j],
    1,
    sqrt(n_grid / V0[j, j]),
    "*"
  )

  Z[!is.finite(Z)] <- NA_real_

  # Configuration graphique, restaurée à la sortie
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  ncol <- min(ncol, length(n_grid))

  par(
    mfrow = c(ceiling(length(n_grid) / ncol), ncol),
    mar = c(3.2, 3.2, 2.2, 0.8),
    mgp = c(1.9, 0.6, 0)
  )

  summary <- data.frame(
    n = n_grid,
    n_valid = integer(length(n_grid)),
    mean = NA_real_,
    sd = NA_real_
  )

  for (i in seq_along(n_grid)) {

    z <- Z[i, ]
    z <- z[is.finite(z)]

    summary$n_valid[i] <- length(z)
    if (length(z) > 0) summary$mean[i] <- mean(z)
    if (length(z) > 1) summary$sd[i] <- sd(z)

    if (length(z) < 2L) {
      plot.new()
      title(main = paste0("n = ", n_grid[i]))
      text(0.5, 0.5, "don't have enough valid replications")
      next
    }

    qqnorm(
      z,
      main = paste0(parameter_label, " - n = ", n_grid[i]),
      xlab = "N(0,1) quantiles",
      ylab = "Observed quantiles",
      pch = 16,
      cex = 0.6
    )

    abline(a = 0, b = 1, col = "red", lwd = 2)
  }

  invisible(summary)
}


table_time <- function(results_mc, column) {

  stopifnot(all(c("n", column) %in% names(results_mc)))

  do.call(rbind, lapply(sort(unique(results_mc$n)), function(n_value) {

    times <- results_mc[[column]][results_mc$n == n_value]
    valid <- is.finite(times) & times >= 0
    x <- times[valid]

    quartiles <- if (length(x) > 0L) {
      quantile(x, probs = c(0.25, 0.50, 0.75), names = FALSE)
    } else {
      rep(NA_real_, 3)
    }

    data.frame(
      n = n_value,
      n_valid = length(x),
      n_invalid = sum(!valid),
      median = quartiles[2],
      q25 = quartiles[1],
      q75 = quartiles[3]
    )
  }))
}

table_coverage <- function(estimates, Vhat, theta_true, n_grid,
                           ci_level = 0.95, mc_level = 0.95) {

  dims <- dim(estimates)
  d <- dims[3]

  stopifnot(
    length(dims) == 3L,
    dims[1] == length(n_grid),
    length(theta_true) == d,
    all(dim(Vhat) == c(dims[1:2], d, d))
  )

  parameter_names <- dimnames(estimates)[[3]]
  if (is.null(parameter_names)) {
    parameter_names <- paste0("theta_", seq_len(d))
  }

  # Diagonales des covariances : [n, réplication, paramètre]
  v_diag <- array(
    matrix(Vhat, nrow = dims[1] * dims[2])[
      , seq.int(1L, d * d, by = d + 1L),
      drop = FALSE
    ],
    dim = dims
  )

  v_diag[!is.finite(v_diag) | v_diag <= 0] <- NA_real_
  se <- sqrt(sweep(v_diag, 1, n_grid, "/"))

  errors <- sweep(estimates, 3, as.numeric(theta_true), "-")
  errors[!is.finite(errors)] <- NA_real_

  covered <- abs(errors) <= qnorm((1 + ci_level) / 2) * se

  m <- apply(!is.na(covered), c(1, 3), sum)
  k <- apply(covered, c(1, 3), sum, na.rm = TRUE)

  coverage <- k / m
  coverage[m == 0] <- NA_real_

  mcse <- sqrt(coverage * (1 - coverage) / m)

  # Intervalle de Wilson
  z <- qnorm((1 + mc_level) / 2)
  denominator <- 1 + z^2 / m
  center <- (coverage + z^2 / (2 * m)) / denominator
  half_width <- z / denominator *
    sqrt(coverage * (1 - coverage) / m + z^2 / (4 * m^2))

  data.frame(
    n = rep(n_grid, times = d),
    parameter = rep(parameter_names, each = length(n_grid)),
    n_valid = as.vector(m),
    n_invalid = as.vector(dims[2] - m),
    estimate = as.vector(coverage),
    mcse = as.vector(mcse),
    lower = as.vector(center - half_width),
    upper = as.vector(center + half_width)
  )
}


summary_coverage <- function(coverage_os, coverage_rml,
                             digits = 1) {

  required <- c("n", "parameter", "estimate")

  stopifnot(
    all(required %in% names(coverage_os)),
    all(required %in% names(coverage_rml))
  )

  # Comparaison sur les mêmes paramètres pour chaque n
  common <- merge(
    coverage_os[, required],
    coverage_rml[, required],
    by = c("n", "parameter"),
    suffixes = c("_os", "_rml")
  )

  ns <- sort(unique(c(coverage_os$n, coverage_rml$n)))

  fmt <- function(x) {
    sprintf(paste0("%.", digits, "f"), 100 * x)
  }

  do.call(rbind, lapply(ns, function(n_value) {

    x <- common[
      common$n == n_value &
        is.finite(common$estimate_os) &
        is.finite(common$estimate_rml),
      ,
      drop = FALSE
    ]

    summarize <- function(values) {
      if (length(values) == 0L) {
        return(c(median = NA_character_, range = NA_character_))
      }

      c(
        median = fmt(median(values)),
        range = paste0(fmt(min(values)), "–", fmt(max(values)))
      )
    }

    os <- summarize(x$estimate_os)
    rml <- summarize(x$estimate_rml)

    data.frame(
      n = n_value,
      n_parameters = nrow(x),
      OS_median = unname(os["median"]),
      RML_median = unname(rml["median"]),
      OS_range = unname(os["range"]),
      RML_range = unname(rml["range"])
    )
  }))
}