# Deux baselines de vraisemblance penalisee, independantes du code SPC / SEMFC.

# Reconstruit la covariance impliquee par les parametres finaux de cv_regsem,
# en utilisant les matrices RAM (A, S, F) du modele lavaan de depart.
regsem_covariance <- function(lav, final_pars, observed_names) {
  ram <- regsem::extractMatrices(lav)
  A <- ram$A_est
  S <- ram$S_est
  a_free <- ram$A > 0
  s_free <- ram$S > 0

  A[a_free] <- final_pars[ram$A[a_free]]
  S[s_free] <- final_pars[ram$S[s_free]]

  F <- ram$Fmat[observed_names, , drop = FALSE]
  inv <- solve(diag(nrow(A)) - A)
  sigma <- F %*% inv %*% S %*% t(inv) %*% t(F)
  dimnames(sigma) <- list(observed_names, observed_names)

  list(sigma = sigma, A = A)
}

# METHODE : regsem, LASSO + choix de penalisation par BIC (cv_regsem).
fit_regsem_baseline <- function(X, models, n_lambda = 23, jump = 0.05) {
  lav <- lavaan::sem(models$sem.model.regsem, data = X, do.fit = FALSE)
  out <- regsem::cv_regsem(
    lav, type = "lasso", pars_pen = models$pen_labels,
    n.lambda = n_lambda, jump = jump, metric = "BIC"
  )

  covariance <- regsem_covariance(lav, out$final_pars, colnames(X))
  loadings <- lapply(seq_along(models$indicators), function(j) {
    drop(covariance$A[models$indicators[[j]], paste0("eta", j)])
  })
  names(loadings) <- paste0("LV", seq_along(loadings))
  list(lambda = loadings, sigma = covariance$sigma, raw_fit = out)
}

# METHODE : lslx, MCP + choix de penalisation par BIC.
fit_lslx_baseline <- function(X, models, penalty_method = "mcp",
                              lambda_grid = exp(seq(log(0.001), log(1), length.out = 10)),
                              delta_grid = c(1, 2, 3, 5, 10, Inf)) {
  out <- lslx::lslx$new(model = models$sem.model.lslx, data = X)
  if (penalty_method == "mcp") {
    out$fit(penalty_method = "mcp", lambda_grid = lambda_grid,
            delta_grid = delta_grid)
  } else {
    out$fit(penalty_method = "lasso", lambda_grid = lambda_grid)
  }

  G <- out$extract_coefficient_matrix(selector = "bic", block = "y<-f")$g
  sigma <- out$extract_implied_cov(selector = "bic")$g
  sigma <- sigma[colnames(X), colnames(X), drop = FALSE]
  loadings <- lapply(seq_along(models$indicators), function(j) {
    drop(G[models$indicators[[j]], paste0("eta", j)])
  })
  names(loadings) <- paste0("LV", seq_along(loadings))
  list(lambda = loadings, sigma = sigma, raw_fit = out)
}
