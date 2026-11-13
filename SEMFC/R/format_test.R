#' Extract Non-Zero Matrix Elements with Row and Column Names
#'
#' Creates a grid of row-column pairs for non-zero elements in a matrix.
#'
#' @param matrix A matrix with named rows and columns.
#'
#' @return Data frame with columns: lhs (row names), rhs (column names).
#'
#' @keywords internal
extract_matrix_pairs <- function(matrix, non_zero = TRUE, rang = NA_integer_) {
    row_labels <- rownames(matrix)
    col_labels <- colnames(matrix)
    if (is.null(row_labels)) {
        row_labels <- paste0("row", seq_len(nrow(matrix)))
    }
    if (is.null(col_labels)) {
        col_labels <- paste0("col", seq_len(ncol(matrix)))
    }

    indices <- if (non_zero) {
        which(matrix != 0, arr.ind = TRUE)
    } else {
        expand.grid(row = seq_len(nrow(matrix)), col = seq_len(ncol(matrix)))
    }

    if (nrow(indices) == 0) {
        return(data.frame(
            lhs = character(0),
            rhs = character(0),
            rang = integer(0)
        ))
    }

    data.frame(
        lhs = row_labels[indices[, 1]],
        rhs = col_labels[indices[, 2]],
        rang = rang
    )
}


#' Create Grid Based on Type
#'
#' Creates a grid of lhs-rhs pairs based on the estimation type.
#'
#' @param type Character string indicating the type: "lambda", "gamma", "beta", "residual", "total", "indirect", or "omega".
#' @param fit_component The relevant component from the fit object.
#'
#' @return Data frame with columns: lhs, rhs, op.
#'
#' @keywords internal
create_grid <- function(type, fit_component) {
    switch(type,
        "lambda" = {
            components <- if (is.list(fit_component)) fit_component else list(fit_component)
            component_names <- names(components)
            if (is.null(component_names)) {
                component_names <- paste0("block", seq_along(components))
            }
            grid <- do.call(rbind, lapply(seq_along(components), function(index) {
                mat <- components[[index]]
                row_labels <- rownames(mat)
                if (is.null(row_labels)) {
                    row_labels <- paste0("indicator", seq_len(nrow(mat)))
                }
                col_labels <- colnames(mat)
                if (is.null(col_labels)) {
                    col_labels <- paste0(component_names[index], "_", seq_len(ncol(mat)))
                }
                data.frame(
                    lhs = rep(col_labels, each = nrow(mat)),
                    rhs = rep(row_labels, times = ncol(mat))
                )
            }))
            grid$op <- "=~"
            grid
        },
        "regression" = {
            components <- if (is.list(fit_component)) fit_component else list(fit_component)
            grid <- do.call(rbind, lapply(seq_along(components), function(rang) {
                extract_matrix_pairs(components[[rang]], rang = rang)
            }))
            grid$op <- "~"
            grid
        },
        "residual" = {
            grid_temp <- do.call(rbind, lapply(names(fit_component), function(lv_name) {
                expand.grid(lhs = lv_name, rhs = gsub("^\\.", "", names(fit_component[[lv_name]])))
            }))
            grid <- data.frame(lhs = grid_temp$rhs, rhs = grid_temp$rhs, op = "~~")
            grid
        },
        "disturbance" = {
            components <- if (is.list(fit_component)) fit_component else list(fit_component)
            grid <- do.call(rbind, lapply(seq_along(components), function(rang) {
                mat <- components[[rang]]
                row_labels <- rownames(mat)
                if (is.null(row_labels)) {
                    row_labels <- paste0("variable", seq_len(nrow(mat)))
                }
                col_labels <- colnames(mat)
                if (is.null(col_labels)) {
                    col_labels <- paste0("variable", seq_len(ncol(mat)))
                }
                data.frame(
                    lhs = rep(row_labels, times = ncol(mat)),
                    rhs = rep(col_labels, each = nrow(mat)),
                    op = "~~",
                    rang = rang
                )
            }))
            grid
        },
        "effect" = {
            grid <- expand.grid(lhs = rownames(fit_component), rhs = colnames(fit_component))
            grid$op <- "~"
            grid
        },
        "omega" = {
            grid <- do.call(rbind, lapply(names(fit_component), function(lv_name) {
                expand.grid(lhs = lv_name, rhs = rownames(fit_component[[lv_name]]))
            }))
            grid$op <- "<~"
            grid
        }
    )
}


#' Format Parameter Estimates Table
#'
#' Creates a standardized data frame with parameter estimates, standard errors,
#' z-scores, confidence intervals and p-values based on the estimation type.
#'
#' @param type Character string indicating the type: "lambda", "regression", "residual", "effect", or "omega".
#' @param fit_component The relevant component from the fit object (e.g., fit$lambda, fit$gamma, fit$beta).
#' @param fit_std_component Standardized estimates component (optional, defaults to NULL).
#' @param se Numeric vector of standard errors (optional, defaults to NA).
#' @param alpha Significance level for confidence intervals (default: 0.05).
#'
#' @return Data frame with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all.
#'   Returns an empty data frame with the same structure if length(est) = 0.
#' @importFrom stats qnorm
#' @keywords internal
format_estimates_table <- function(type, fit_component, fit_std_component = NULL, se = NA, alpha = 0.05, digits = 3) {
    components <- if (is.list(fit_component)) fit_component else list(fit_component)
    std_components <- if (is.null(fit_std_component)) {
        NULL
    } else if (is.list(fit_std_component)) {
        fit_std_component
    } else {
        list(fit_std_component)
    }
    se_components <- if (is.list(se)) se else list(se)

    flatten <- function(values, type, non_zero = FALSE, reference = values) {
        unlist(lapply(seq_along(values), function(index) {
            mat <- values[[index]]
            if (non_zero) {
                mat[reference[[index]] != 0]
            } else {
                as.vector(mat)
            }
        }), use.names = FALSE)
    }

    non_zero <- type == "regression"
    est <- flatten(components, type, non_zero, components)


    # Return empty data frame if no estimates
    if (length(est) == 0) {
        return(data.frame(
            lhs      = character(0),
            op       = character(0),
            rhs      = character(0),
            rang     = integer(0),
            est      = numeric(0),
            se       = numeric(0),
            z        = numeric(0),
            ci.lower = numeric(0),
            ci.upper = numeric(0),
            pvalue   = numeric(0),
            std.all  = numeric(0)
        ))
    }

    grid <- create_grid(type, fit_component)


    std.all <- if (!is.null(std_components)) {
        flatten(std_components, type, non_zero, components)
    } else {
        rep(NA_real_, length(est))
    }

    se <- if (is.null(se) || (length(se) == 1 && is.na(se))) {
        rep(NA_real_, length(est))
    } else {
        flatten(se_components, type, non_zero, components)
    }

    z_critical <- qnorm(1 - alpha / 2)
    if (all(is.na(se))) {
        se <- rep(NA, length(est))
        z_scores <- rep(NA, length(est))
        ci_lower <- rep(NA, length(est))
        ci_upper <- rep(NA, length(est))
        pvalue <- rep(NA, length(est))
    } else {
        z_scores <- est / se
        ci_lower <- est - z_critical * se
        ci_upper <- est + z_critical * se
        pvalue <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
    }

    data.frame(
        lhs = grid$lhs,
        op = grid$op,
        rhs = grid$rhs,
        rang = if ("rang" %in% names(grid)) grid$rang else NA_integer_,
        est = round(est, digits),
        se = round(se, digits),
        z = round(z_scores, digits),
        ci.lower = round(ci_lower, digits),
        ci.upper = round(ci_upper, digits),
        pvalue = round(pvalue, digits),
        std.all = round(std.all, digits)
    )
}

#' Format Model Estimates
#'
#' Processes a fitted model object to extract and format parameter estimates
#' into standardized data frames. Creates tables for loadings, path coefficients,
#' regression coefficients, residual variances, total effects, indirect effects,
#' and composite weights.
#'
#' @param fit A fitted model object containing the following components:
#'   - `lambda`: A named list of loadings for each latent variable.
#'   - `std_lambda`: Standardized loadings (optional).
#'   - `gamma`: A matrix of path coefficients between latent variables.
#'   - `beta`: A matrix of regression coefficients.
#'   - `residual_variance`: A named list of residual variances.
#'   - `effect$total_effect`: A matrix of total effects.
#'   - `effect$indirect_effect`: A matrix of indirect effects.
#'   - `omega`: A list of matrices containing weights of indicators for composite blocks.
#' @param se_list A list containing standard errors for each component (optional, defaults to empty list).
#'   Expected elements: `sd_lambda`, `sd_gamma`, `sd_beta`, `sd_residual_variance`,
#'   `sd_total_effects`, `sd_indirect_effects`, `sd_omega`.
#'
#' @return A list of data frames with columns: lhs, op, rhs, est, se, z, ci.lower, ci.upper, pvalue, std.all:
#'   - `lambda`: Loadings (op = "=~").
#'   - `gamma`: Path coefficients between latent variables (op = "~").
#'   - `beta`: Regression coefficients (op = "~").
#'   - `residual_variance`: Residual variances (op = "~~").
#'   - `total_effects`: Total effects (op = "~").
#'   - `indirect_effects`: Indirect effects (op = "~").
#'   - `omega`: Composite weights (op = "<~").
#'
#' @keywords internal

formatting_estimate <- function(fit, se_list = list(), digits = 3) {
    table_lambda <- format_estimates_table(
        type = "lambda",
        fit_component = fit$li_Lambda,
        fit_std_component = fit$li_Lambda_star,
        se = se_list$li_sd_Lambda,
        digits = digits
    )

    # Gamma
    table_gamma <- format_estimates_table(
        type = "regression",
        fit_component = fit$li_gamma,
        se = se_list$li_sd_gamma,
        digits = digits
    )

    # Beta
    table_beta <- format_estimates_table(
        type = "regression",
        fit_component = fit$li_beta,
        se = se_list$li_sd_beta,
        digits = digits
    )

    # Psi
    result <- list(
        lambda = table_lambda,
        gamma = table_gamma,
        beta = table_beta
    )

    if (!is.null(fit$li_PSI)) {
        result$psi <- format_estimates_table(
            type = "disturbance",
            fit_component = fit$li_PSI,
            se = se_list$li_sd_psi,
            digits = digits
        )
    }

    result
}
