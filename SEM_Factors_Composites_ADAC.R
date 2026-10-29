##################################################################
# Reproducible results for Tenenhaus, Tenenhaus & Dijkstra paper #
#              submitted to ADAC in July 15th, 2024              #
# Work alex
##################################################################

#############################################
# Remove all objects from the R environment #
#############################################
rm(list = ls())

######################################
# Load useful packages and functions #
######################################

# packages
library(readxl)
library(roxygen2)
library(devtools)
library(Matrix)
library(knitr)
library(pheatmap)
library(mvtnorm)
library(MASS)
library(lavaan)
library(Rsolnp)

setwd("~/Documents/SEMFC")
load_all("SEMFC")


# functions

# source("functions/data_simulation.R")
source("functions/F1.R")
source("functions/h_theta.R")
source("functions/s_implied.R")
source("functions/h_constraints.R")
source("functions/reconstruct_blocs.R")
source("functions/calculate_Sigma.R")

#############################################
########## MONTE-CARLO SIMULATION ###########
#############################################
set.seed(4)
n_simu <- 100
N <- 100000
R_try <- 2
R_true <- 2
n_signif <- 3
do_ML <- FALSE


ortho_generate <- FALSE
ortho_solve_SVD <- FALSE
ortho_solve_ML <- FALSE
diag_generate <- TRUE

show_svd <- FALSE
show_ml <- FALSE
show_TRUE <- FALSE
show_norme <- FALSE
show_svd <- TRUE
show_ml <- TRUE
show_TRUE <- TRUE

do_full_reorder <- TRUE
ancien_P <- FALSE

triche <- FALSE
no_svd_refining <- FALSE
mauvais_depart <- FALSE


source("functions/data_simulation_dim_n.R")


sol_svd <- matrix(0, 61, n_simu)
sol_ml <- matrix(0, 61, n_simu)


lili_lambda_svd <- list()
lili_beta_svd <- list()
lili_gamma_svd <- list()
lili_psi_svd <- list()
lili_beta_svd <- list()
lili_gamma_svd <- list()
li_P_svd <- list()
lili_psi_svd <- list()
error_SVD <- rep(NA, n_simu)
ecart_true_empirique <- rep(NA, n_simu)
ecart_empirique_svd <- rep(NA, n_simu)
ecarts_sigma_per_bloc <- matrix(NA, length(is_endogenes), n_simu)
lili_theta_svd <- list()
lili_phi_endo_svd <- list()
lili_phi_exo_svd <- list()
lili_S <- list()

lili_lambda_ml <- list()
lili_phi_exo_ml <- list()
lili_phi_endo_ml <- list()
lili_theta_ml <- list()
lili_gamma_ml <- list()
lili_beta_ml <- list()
lili_is_dag_ml <- list()
lili_correl_same_diag_ml <- list()

# Maximum likelihood
lambda_hat_ml <- matrix(0, 18, n_simu) # les correlations entre facteurs et y
omega_hat_ml <- matrix(0, 12, n_simu) # les w de poids pour les blocs formatifs
rho_hat_ml <- matrix(0, 15, n_simu) # les correlations entre facteurs latents
beta_hat_ml <- matrix(0, 2, n_simu) # les  elements libres de B
gamma_hat_ml <- matrix(0, 4, n_simu) # Les elements libres de Gamma
r2_hat_ml <- matrix(0, 2, n_simu) # Le R^2 de chaque facteur endogene
psi_hat_ml <- matrix(0, 3, n_simu) # les elements libres de Psi
sigma_hat <- matrix(0, 12, n_simu) # Plusieurs ecarts quadratiques definis plus loin entre valeurs trouvee de Sigma sur blocs formatifs
var_hat_ml <- matrix(0, 6, n_simu) # Les elements libres de Theta
std_err_ml <- matrix(NA, 61, n_simu) # Estimation des standard errors des parametres ML
f_ml <- rep(0, n_simu) # goodness of fit?
param_ml <- matrix(0, 61, n_simu)
lrt_pval <- rep(NA, n_simu)
vcov_pval <- rep(NA, n_simu)
error_ML <- rep(NA, n_simu)


# Goodness of fit
gof <- matrix(NA, 2, n_simu)

J <- length(is_endogenes)

for (b in seq_len(n_simu)) {
  Sigma_SVD <- NULL
  fit.svd <- NULL
  if (b %% 10 == 0) print(b)
  # try(
  #   {

  # print(dim(SIGMA_TRUE))
  # print(sum(vec_indicators_per_bloc))
  X <- mvrnorm(N, rep(0, sum(vec_indicators_per_bloc)), SIGMA_TRUE, empirical = FALSE)
  S <- cov(X)
  if (triche) {
    S <- SIGMA_TRUE
  }
  ecart_true_empirique[b] <- norm(as.matrix(S) - as.matrix(SIGMA_TRUE), "F")

  indicator_blocks <- split(
    seq_len(ncol(X)),
    rep(seq_along(vec_indicators_per_bloc), vec_indicators_per_bloc)
  )
  Y <- setNames(
    lapply(indicator_blocks, function(indices) X[, indices, drop = FALSE]),
    paste0("X", seq_along(indicator_blocks))
  )


  li_C <- lapply(1:R_try, function(r) {
    if (r <= R_true) {
      C <- li_C_TRUE[[r]]
    } else {
      C <- li_C_TRUE[[R_true]]
    }
  })

  # fit.svd <- svdSEM_multi(Y, li_C,
  #   scale = FALSE,
  #   mode = rep("reflective", n_blocs),
  #   bias = FALSE
  # )

  tryCatch(
    {
      # if (b <= 7) {
      #   stop("stop svd")
      # }
      fit.svd <- svdSEM_multi(Y, li_C,
        scale = FALSE,
        mode = rep("reflective", n_blocs),
        bias = FALSE
      )

      lili_lambda_svd[[b]] <- fit.svd$li_Lambda
      lili_beta_svd[[b]] <- fit.svd$li_beta
      lili_gamma_svd[[b]] <- fit.svd$li_gamma
      lili_psi_svd[[b]] <- fit.svd$li_psi
      lili_theta_svd[[b]] <- fit.svd$li_theta
      li_P_svd[[b]] <- fit.svd$P_IMPLIED
      lili_beta_svd[[b]] <- fit.svd$li_beta
      lili_gamma_svd[[b]] <- fit.svd$li_gamma
      lili_psi_svd[[b]] <- fit.svd$li_PSI
      lili_phi_endo_svd[[b]] <- fit.svd$li_phi_endo
      lili_phi_exo_svd[[b]] <- fit.svd$li_phi_exo
      lili_S[[b]] <- fit.svd$li_S

      for (j in 1:J) {
        Sigma_svd_ii <- fit.svd$SIGMA_IMPLIED[sum(vec_indicators_per_bloc[1:(j - 1)]) + 1:vec_indicators_per_bloc[j], sum(vec_indicators_per_bloc[1:(j - 1)]) + 1:vec_indicators_per_bloc[j]]
        Sigma_true_ii <- SIGMA_TRUE[sum(vec_indicators_per_bloc[1:(j - 1)]) + 1:vec_indicators_per_bloc[j], sum(vec_indicators_per_bloc[1:(j - 1)]) + 1:vec_indicators_per_bloc[j]]
        ecart_ii <- norm(as.matrix(Sigma_svd_ii) - as.matrix(Sigma_true_ii), "F")
        ecarts_sigma_per_bloc[j, b] <- ecart_ii
      }
      Sigma_SVD <- fit.svd$SIGMA_IMPLIED

      ## TESTING
      li_Lambda <- fit.svd$li_Lambda
      mat_Lambda <- Matrix::bdiag(li_Lambda)
      P_tilde_approx <- t(mat_Lambda) %*% SIGMA_TRUE %*% mat_Lambda
      # print(round(as.matrix(P_tilde_approx), 2))
      P_tilde_obtained <- t(mat_Lambda) %*% Sigma_SVD %*% mat_Lambda
      # print(round(as.matrix(P_tilde_obtained), 2))
      # print(eigen(as.matrix(P_tilde_obtained), symmetric = TRUE)$values)
      # print(round(as.matrix(abs(P_tilde_approx - P_tilde_obtained)), 2))


      # relat_P <- round(as.matrix(abs(fit.svd$P_tilde - P_TRUE) / abs(P_TRUE)), 2)
      # mask <- as.matrix(P_TRUE) < 1e-2
      # relat_P[mask] <- 0
      # print(paste("Relative error:", norm(as.matrix(relat_P), "F")))

      # print("relative error")
      # relat <- round(as.matrix(abs(Sigma_SVD - SIGMA_TRUE) / abs(SIGMA_TRUE)), 2)
      # mask <- as.matrix(SIGMA_TRUE) < 1e-2
      # relat[mask] <- 0
      # print(norm(as.matrix(relat), "F"))
      ### END TESTING


      erreur_abs <- norm(as.matrix(Sigma_SVD) - as.matrix(SIGMA_TRUE), "F")

      error_SVD[b] <- erreur_abs
      ecart_empirique_svd[b] <- norm(as.matrix(Sigma_SVD) - as.matrix(S), "F")
    },
    error = function(e) {
      message("Erreur dans svd à l'itération ", b)
      warning(paste("Erreur dans svd à l'itération", b, ":", conditionMessage(e)))
    }
  )


  if (do_ML) {
    tryCatch(
      {
        if (is.null(Sigma_SVD)) {
          # warning(paste("SVD failed at iteration", b, "- ML estim from scratch"))
          init_base <- TRUE
        } else {
          min_eigenvalue <- min(eigen(Sigma_SVD, symmetric = TRUE, only.values = TRUE)$values)
          if (min_eigenvalue < 1e-8) {
            warning(paste("SVD implied covariance matrix is not positive definite at iteration", b, "value:", min_eigenvalue, "< 1e-8 - ML estim from scratch"))
            init_base <- TRUE
          } else {
            init_base <- FALSE
          }
        }
        if (!is.null(fit.svd)) {
          big_mat_Lambda <- do.call(cbind, fit.svd$li_Lambda)
          big_vec_Lambda <- as.vector(big_mat_Lambda) # colonne par colonne : ordre lexico j puis r
        }

        if (mauvais_depart | init_base) {
          li_vec_lambda <- lapply(1:J, function(j) {
            Lambda_j <- matrix(0, nrow = vec_indicators_per_bloc[j], ncol = R_try)
            for (r in 1:R_try) {
              delta <- rep(0, vec_indicators_per_bloc[j])
              delta[r] <- 1
              Lambda_j[, r] <- delta
            }
            vec_lambda_j <- as.vector(Lambda_j)
            return(vec_lambda_j)
          })
          big_vec_Lambda <- do.call(c, li_vec_lambda)
        }

        ####
        if (!is.null(fit.svd)) {
          li_vec_phi_exo <- lapply(1:R_try, function(r) {
            phi_exo <- fit.svd$li_phi_exo[[r]]
            upper_tri <- phi_exo[upper.tri(phi_exo, diag = FALSE)]
            return(upper_tri)
          })
        }
        n_exogenes <- sum(!is_endogenes)
        n_endogenes <- sum(is_endogenes)


        if (mauvais_depart | init_base) {
          li_vec_phi_exo <- lapply(1:R_try, function(r) {
            phi_exo <- matrix(0, nrow = n_exogenes, ncol = n_exogenes)
            upper_tri <- phi_exo[upper.tri(phi_exo, diag = FALSE)]
            upper_tri[] <- 0
            return(upper_tri)
          })
        }

        big_vec_phi_exo <- do.call(c, li_vec_phi_exo)

        if (!is.null(fit.svd)) {
          li_vec_phi_endo <- lapply(1:R_try, function(r) {
            phi_endo <- fit.svd$li_phi_endo[[r]]
            upper_tri <- phi_endo[upper.tri(phi_endo, diag = FALSE)]
            gr <- fit.svd$li_gr[[r]]
            if (!igraph::is_dag(gr)) {
              return(upper_tri)
            } else {
              return(c())
            }
            return(upper_tri)
          })
        }

        li_gr <- lapply(li_C, function(C) igraph::graph_from_adjacency_matrix(C))

        if (mauvais_depart | init_base) {
          li_vec_phi_endo <- lapply(1:R_try, function(r) {
            phi_endo <- matrix(0, nrow = n_endogenes, ncol = n_endogenes)
            upper_tri <- phi_endo[upper.tri(phi_endo, diag = FALSE)]
            gr <- li_gr[[r]]
            if (!igraph::is_dag(gr)) {
              upper_tri[] <- 0
              return(upper_tri)
            } else {
              return(c())
            }
          })
        }

        big_vec_phi_endo <- do.call(c, li_vec_phi_endo)

        if (!is.null(fit.svd)) {
          theta <- do.call(c, fit.svd$li_theta)
        }
        if (mauvais_depart | init_base) {
          theta <- rep(0.1, sum(vec_indicators_per_bloc))
        }


        li_non_null_beta_gamma <- lapply(1:R_try, function(r) {
          # non_null_beta <- matrix(0, nrow = nrow(li_beta[[r]]), ncol = ncol(li_beta[[r]]))
          # non_null_gamma <- matrix(0, nrow = nrow(li_gamma[[r]]), ncol = ncol(li_gamma[[r]]))
          C <- li_C[[r]]
          exo_endo <- ind_exo_endo(C)
          ind_exo <- exo_endo$ind_exo
          ind_endo <- exo_endo$ind_endo
          C_endo <- C[ind_endo, ind_endo]
          C_exo_to_endo <- C[ind_exo, ind_endo]
          non_null_gamma <- t(C_exo_to_endo)
          non_null_beta <- t(C_endo)
          return(list(non_null_gamma = non_null_gamma, non_null_beta = non_null_beta))
        })

        if (!is.null(fit.svd)) {
          li_beta <- fit.svd$li_beta
          li_vec_beta <- lapply(1:R_try, function(r) {
            non_null_beta <- li_non_null_beta_gamma[[r]]$non_null_beta
            beta_r <- li_beta[[r]]
            beta_r[non_null_beta == 0] <- NA
            vec_beta_r <- as.vector(beta_r)
            vec_beta_r <- vec_beta_r[!is.na(vec_beta_r)]
            return(vec_beta_r)
          })

          vec_beta <- do.call(c, li_vec_beta)
        }


        if (mauvais_depart | init_base) {
          vec_beta <- rep(0, sum(sapply(li_non_null_beta_gamma, function(x) sum(x$non_null_beta))))
        }

        if (!is.null(fit.svd)) {
          li_gamma <- fit.svd$li_gamma
          li_vec_gamma <- lapply(1:R_try, function(r) {
            non_null_gamma <- li_non_null_beta_gamma[[r]]$non_null_gamma
            gamma_r <- li_gamma[[r]]
            gamma_r[non_null_gamma == 0] <- NA
            vec_gamma_r <- as.vector(gamma_r)
            vec_gamma_r <- vec_gamma_r[!is.na(vec_gamma_r)]
            return(vec_gamma_r)
          })
          vec_gamma <- do.call(c, li_vec_gamma)
        }
        if (mauvais_depart | init_base) {
          vec_gamma <- rep(0, sum(sapply(li_non_null_beta_gamma, function(x) sum(x$non_null_gamma))))
        }


        if (!is.null(fit.svd)) {
          P_IMPLIED <- fit.svd$P_IMPLIED
          li_correl_same_diag <- lapply(1:J, function(j) {
            index_low <- R_try * (j - 1) + 1
            index_high <- R_try * j
            diag_bloc <- P_IMPLIED[index_low:index_high, index_low:index_high]
            diag_vec <- diag_bloc[upper.tri(diag_bloc)]
            return(diag_vec)
          })
        }

        if (mauvais_depart | init_base) {
          li_correl_same_diag <- lapply(1:J, function(j) {
            index_low <- R_try * (j - 1) + 1
            index_high <- R_try * j
            diag_bloc <- diag(1, nrow = R_try, ncol = R_try)
            diag_vec <- diag_bloc[upper.tri(diag_bloc)]
            diag_vec[] <- 0
            return(diag_vec)
          })
        }

        vec_correl_same_diag <- do.call(c, li_correl_same_diag)


        init_ml_with_S <- c(big_vec_Lambda, big_vec_phi_exo, big_vec_phi_endo, theta, vec_gamma, vec_beta, vec_correl_same_diag)

        vrais_at_opt_SVD <- F1(init_ml_with_S, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
        print(paste("vrais_at_opt SVD", vrais_at_opt_SVD))

        # stop("F1")

        if (ortho_solve_ML) {
          n_zeros_eq <- (J * R_try * (R_try - 1)) %/% 2
        } else {
          n_zeros_eq <- 0
        }

        which_theta <- (length(big_vec_Lambda) + length(big_vec_phi_exo) + length(big_vec_phi_endo) + 1):(length(big_vec_Lambda) + length(big_vec_phi_exo) + length(big_vec_phi_endo) + length(theta))

        ineq_LB <- rep(-Inf, length(init_ml_with_S))
        ineq_LB[which_theta] <- 0


        setwd("~/Documents/SEMFC/SEMFC/easy_solver")
        source("easy_solver.R")

        # if (b <= 7) {
        #   stop("stop ml")
        # }

        fit.ml <- easy_solver(
          pars = init_ml_with_S,
          fun = F1,
          eqfun = heq1,
          eqB = rep(0, n_zeros_eq),
          ineqLB = ineq_LB,
          S = S,
          J = J,
          R_try = R_try,
          vec_indicators_per_bloc = vec_indicators_per_bloc,
          li_gr = li_gr,
          li_non_null_beta_gamma = li_non_null_beta_gamma,
          maxit = 1500,
          kkt_tol = 1e-5,
          stepsize = 1e-3,
          max_step = 10,
          f_min = -1e-10,
          trace = 10
        )

        convergence <- fit.ml$convergence
        if (convergence != 0) {
          warning(paste("Convergence issue at iteration", b, "with convergence code:", convergence))
        }
        opt_param <- fit.ml$pars
        opt_vrais <- fit.ml$values[length(fit.ml$values)]

        reconstructed_params <- reconstruct_blocs(opt_param, S, J, R_try, vec_indicators_per_bloc, li_gr, li_non_null_beta_gamma)
        li_Lambda <- reconstructed_params$li_Lambda
        li_phi_exo <- reconstructed_params$li_phi_exo
        li_phi_endo <- reconstructed_params$li_phi_endo
        vec_theta <- reconstructed_params$vec_theta
        li_gamma <- reconstructed_params$li_gamma
        li_beta <- reconstructed_params$li_beta
        is_dag <- reconstructed_params$is_dag
        li_correl_same_diag <- reconstructed_params$li_correl_same_diag

        lili_lambda_ml[[b]] <- li_Lambda
        lili_phi_exo_ml[[b]] <- li_phi_exo
        lili_phi_endo_ml[[b]] <- li_phi_endo
        lili_theta_ml[[b]] <- vec_theta
        lili_gamma_ml[[b]] <- li_gamma
        lili_beta_ml[[b]] <- li_beta
        lili_is_dag_ml[[b]] <- is_dag
        lili_correl_same_diag_ml[[b]] <- li_correl_same_diag
        Sigma_ML <- calculate_Sigma(li_Lambda, li_phi_exo, li_phi_endo, vec_theta, li_gamma, li_beta, is_dag, li_correl_same_diag, li_C, is_endogenes)

        error_ML[b] <- norm(as.matrix(Sigma_ML) - as.matrix(SIGMA_TRUE), "F")

        print(paste("opt_vrais", opt_vrais))
      },
      error = function(e) {
        message("Erreur dans ML à l'itération ", b, ": ", conditionMessage(e))
        warning(paste("Erreur dans ML à l'itération", b, ":", conditionMessage(e)))
      }
    )
  }
}

### On calcule les stats SVD

li_P_svd_non_null <- li_P_svd[!sapply(li_P_svd, is.null)]
P_tridim <- simplify2array(li_P_svd_non_null)
P_sd <- apply(P_tridim, c(1, 2), sd, na.rm = TRUE)
P_mean <- apply(P_tridim, c(1, 2), mean, na.rm = TRUE)

lili_Lambda_aggregate <- lapply(1:J, function(j) {
  li_Lambda_j <- lapply(lili_lambda_svd, function(x) x[[j]])
  li_Lambda_j_non_null <- li_Lambda_j[!sapply(li_Lambda_j, is.null)]
  li_Lambda_j_tridim <- simplify2array(li_Lambda_j_non_null)
  Lambda_j_sd <- round(apply(li_Lambda_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  Lambda_j_mean <- round(apply(li_Lambda_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(Lambda_j_sd = Lambda_j_sd, Lambda_j_mean = Lambda_j_mean))
})

li_lambda_norm_TRUE <- lapply(li_lambda_TRUE, function(mat) {
  vec_norm <- apply(mat, 2, function(x) sqrt(sum(x^2)))
  return(vec_norm)
})

lili_norm_lambda_aggregate <- lapply(lili_Lambda_aggregate, function(mat) {
  vec_norm <- apply(mat$Lambda_j_mean, 2, function(x) sqrt(sum(x^2)))
  return(vec_norm)
})

lili_S_aggregate <- lapply(1:J, function(i) {
  lapply(1:J, function(j) {
    mats <- lapply(lili_S, function(x) x[[i]][[j]])
    mats <- mats[!sapply(mats, is.null)]
    arr <- simplify2array(mats)
    list(
      mean = apply(arr, c(1, 2), mean, na.rm = TRUE),
      sd   = apply(arr, c(1, 2), sd, na.rm = TRUE)
    )
  })
})

lili_gamma_aggregate <- lapply(1:R_try, function(r) {
  li_gamma_j <- lapply(lili_gamma_svd, function(x) x[[r]])
  li_gamma_j_non_null <- li_gamma_j[!sapply(li_gamma_j, is.null)]
  li_gamma_j_tridim <- simplify2array(li_gamma_j_non_null)
  gamma_j_sd <- round(apply(li_gamma_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  gamma_j_mean <- round(apply(li_gamma_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(gamma_j_sd = gamma_j_sd, gamma_j_mean = gamma_j_mean))
})

lili_beta_aggregate <- lapply(1:R_try, function(r) {
  li_beta_j <- lapply(lili_beta_svd, function(x) x[[r]])
  li_beta_j_non_null <- li_beta_j[!sapply(li_beta_j, is.null)]
  li_beta_j_tridim <- simplify2array(li_beta_j_non_null)
  beta_j_sd <- round(apply(li_beta_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  beta_j_mean <- round(apply(li_beta_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(beta_j_sd = beta_j_sd, beta_j_mean = beta_j_mean))
})

lili_psi_aggregate <- lapply(1:R_try, function(r) {
  li_psi_j <- lapply(lili_psi_svd, function(x) x[[r]])
  li_psi_j_non_null <- li_psi_j[!sapply(li_psi_j, is.null)]
  li_psi_j_tridim <- simplify2array(li_psi_j_non_null)
  psi_j_sd <- round(apply(li_psi_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  psi_j_mean <- round(apply(li_psi_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(psi_j_sd = psi_j_sd, psi_j_mean = psi_j_mean))
})

lili_endo_aggregate <- lapply(1:R_try, function(r) {
  li_endo_j <- lapply(lili_phi_endo_svd, function(x) x[[r]])
  li_endo_j_non_null <- li_endo_j[!sapply(li_endo_j, is.null)]
  li_endo_j_tridim <- simplify2array(li_endo_j_non_null)
  endo_j_sd <- round(apply(li_endo_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  endo_j_mean <- round(apply(li_endo_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(endo_j_sd = endo_j_sd, endo_j_mean = endo_j_mean))
})

lili_exo_aggregate <- lapply(1:R_try, function(r) {
  li_exo_j <- lapply(lili_phi_exo_svd, function(x) x[[r]])
  li_exo_j_non_null <- li_exo_j[!sapply(li_exo_j, is.null)]
  li_exo_j_tridim <- simplify2array(li_exo_j_non_null)
  exo_j_sd <- round(apply(li_exo_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
  exo_j_mean <- round(apply(li_exo_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
  return(list(exo_j_sd = exo_j_sd, exo_j_mean = exo_j_mean))
})


# li_Lambda_star_2 <- lapply(lili_lambda_svd, function(x) x[[2]] / sum(x[[2]]^2)^0.5)
# li_lambda_star_2_non_null <- li_Lambda_star_2[!sapply(li_Lambda_star_2, is.null)]
# li_lambda_star_2_tridim <- simplify2array(li_lambda_star_2_non_null)
# print(li_lambda_star_2_tridim)
# print(li_lambda_star_2_tridim[, , 1])
# Lambda_star_2_sd <- apply(li_lambda_star_2_tridim, c(1, 2), sd, na.rm = TRUE)
# Lambda_star_2_mean <- apply(li_lambda_star_2_tridim, c(1, 2), mean, na.rm = TRUE)

li_beta_1 <- lapply(lili_beta_svd, function(x) x[[1]])
li_beta_1_non_null <- li_beta_1[!sapply(li_beta_1, is.null)]
li_beta_1_tridim <- simplify2array(li_beta_1_non_null)
beta_1_sd <- apply(li_beta_1_tridim, c(1, 2), sd, na.rm = TRUE)
beta_1_mean <- apply(li_beta_1_tridim, c(1, 2), mean, na.rm = TRUE)

li_gamma_1 <- lapply(lili_gamma_svd, function(x) x[[1]])
li_gamma_1_non_null <- li_gamma_1[!sapply(li_gamma_1, is.null)]
li_gamma_1_tridim <- simplify2array(li_gamma_1_non_null)
gamma_1_sd <- apply(li_gamma_1_tridim, c(1, 2), sd, na.rm = TRUE)
gamma_1_mean <- apply(li_gamma_1_tridim, c(1, 2), mean, na.rm = TRUE)

li_theta_svd <- lapply(1:J, function(j) {
  li_theta_j <- lapply(lili_theta_svd, function(x) x[[j]])
  li_theta_j_non_null <- li_theta_j[!sapply(li_theta_j, is.null)]
  li_theta_j_bidim <- simplify2array(li_theta_j_non_null)
  theta_j_sd <- round(apply(li_theta_j_bidim, c(1), sd, na.rm = TRUE), n_signif)
  theta_j_mean <- round(apply(li_theta_j_bidim, c(1), mean, na.rm = TRUE), n_signif)
  return(list(theta_j_sd = theta_j_sd, theta_j_mean = theta_j_mean))
})


### On calcule les stats ML

if (do_ML) {
  lili_Lambda_ml_aggregate <- lapply(1:J, function(j) {
    li_Lambda_j <- lapply(lili_lambda_ml, function(x) x[[j]])
    li_Lambda_j_non_null <- li_Lambda_j[!sapply(li_Lambda_j, is.null)]
    li_Lambda_j_tridim <- simplify2array(li_Lambda_j_non_null)
    Lambda_j_sd <- round(apply(li_Lambda_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    Lambda_j_mean <- round(apply(li_Lambda_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(Lambda_j_sd = Lambda_j_sd, Lambda_j_mean = Lambda_j_mean))
  })

  lili_norm_lambda_ml_aggregate <- lapply(lili_Lambda_ml_aggregate, function(mat) {
    vec_norm <- apply(mat$Lambda_j_mean, 2, function(x) sqrt(sum(x^2)))
    return(vec_norm)
  })

  lili_gamma_ml_aggregate <- lapply(1:R_try, function(r) {
    li_gamma_j <- lapply(lili_gamma_ml, function(x) x[[r]])
    li_gamma_j_non_null <- li_gamma_j[!sapply(li_gamma_j, is.null)]
    li_gamma_j_tridim <- simplify2array(li_gamma_j_non_null)
    gamma_j_sd <- round(apply(li_gamma_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    gamma_j_mean <- round(apply(li_gamma_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(gamma_j_sd = gamma_j_sd, gamma_j_mean = gamma_j_mean))
  })

  lili_beta_ml_aggregate <- lapply(1:R_try, function(r) {
    li_beta_j <- lapply(lili_beta_ml, function(x) x[[r]])
    li_beta_j_non_null <- li_beta_j[!sapply(li_beta_j, is.null)]
    li_beta_j_tridim <- simplify2array(li_beta_j_non_null)
    beta_j_sd <- round(apply(li_beta_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    beta_j_mean <- round(apply(li_beta_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(beta_j_sd = beta_j_sd, beta_j_mean = beta_j_mean))
  })

  lili_theta_ml_non_null <- lili_theta_ml[!sapply(lili_theta_ml, is.null)]
  theta_ml_tridim <- simplify2array(lili_theta_ml_non_null)
  theta_ml_sd <- round(apply(theta_ml_tridim, 1, sd, na.rm = TRUE), n_signif)
  theta_ml_mean <- round(apply(theta_ml_tridim, 1, mean, na.rm = TRUE), n_signif)

  lili_endo_ml_aggregate <- lapply(1:R_try, function(r) {
    li_endo_j <- lapply(lili_phi_endo_ml, function(x) x[[r]])
    li_endo_j_non_null <- li_endo_j[!sapply(li_endo_j, is.null)]
    if (length(li_endo_j_non_null) == 0) {
      return(list(endo_j_sd = NULL, endo_j_mean = NULL))
    }
    li_endo_j_tridim <- simplify2array(li_endo_j_non_null)
    endo_j_sd <- round(apply(li_endo_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    endo_j_mean <- round(apply(li_endo_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(endo_j_sd = endo_j_sd, endo_j_mean = endo_j_mean))
  })

  lili_exo_ml_aggregate <- lapply(1:R_try, function(r) {
    li_exo_j <- lapply(lili_phi_exo_ml, function(x) x[[r]])
    li_exo_j_non_null <- li_exo_j[!sapply(li_exo_j, is.null)]
    li_exo_j_tridim <- simplify2array(li_exo_j_non_null)
    exo_j_sd <- round(apply(li_exo_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    exo_j_mean <- round(apply(li_exo_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(exo_j_sd = exo_j_sd, exo_j_mean = exo_j_mean))
  })

  lili_is_dag_ml_non_null <- lili_is_dag_ml[!sapply(lili_is_dag_ml, is.null)]
  is_dag_ml_tridim <- simplify2array(lili_is_dag_ml_non_null)
  is_dag_ml_mean <- apply(is_dag_ml_tridim, 1, mean, na.rm = TRUE)

  lili_correl_same_diag_ml_aggregate <- lapply(1:J, function(j) {
    li_correl_j <- lapply(lili_correl_same_diag_ml, function(x) x[[j]])
    li_correl_j_non_null <- li_correl_j[!sapply(li_correl_j, is.null)]
    li_correl_j_tridim <- simplify2array(li_correl_j_non_null)
    correl_j_sd <- round(apply(li_correl_j_tridim, c(1, 2), sd, na.rm = TRUE), n_signif)
    correl_j_mean <- round(apply(li_correl_j_tridim, c(1, 2), mean, na.rm = TRUE), n_signif)
    return(list(correl_j_sd = correl_j_sd, correl_j_mean = correl_j_mean))
  })
}

print(paste("Nombre d'itérations SVD échouées :", sum(sapply(lili_lambda_svd, is.null))))
print(paste("erreur SVD moyenne L2 sur Sigma:", mean(error_SVD, na.rm = TRUE)))

if (do_ML) {
  print(paste("Nombre d'itérations ML échouées :", sum(sapply(lili_lambda_ml, is.null))))
  print(paste("erreur ML moyenne L2 sur Sigma:", mean(error_ML, na.rm = TRUE)))
}


if (show_svd) {
  print("SVD Lambda")
  print(lili_Lambda_aggregate)

  print("SVD normes Lambda")
  print(lili_norm_lambda_aggregate)

  print("SVD Phi endo")
  print(lili_endo_aggregate)

  print("SVD Phi exo")
  print(lili_exo_aggregate)

  print("SVD Gamma")
  print(lili_gamma_aggregate)

  print("SVD Beta")
  print(lili_beta_aggregate)

  print("SVD Theta")
  print(li_theta_svd)

  print("SVD Psi")
  print(lili_psi_aggregate)
}

if (show_ml & do_ML) {
  print("ML Lambda")
  print(lili_Lambda_ml_aggregate)

  print("ML normes Lambda")
  print(lili_norm_lambda_ml_aggregate)

  print("ML Phi endo")
  print(lili_endo_ml_aggregate)

  print("ML Phi exo")
  print(lili_exo_ml_aggregate)

  print("ML Gamma")
  print(lili_gamma_ml_aggregate)

  print("ML Beta")
  print(lili_beta_ml_aggregate)

  print("ML Theta")
  print(list(theta_ml_sd = theta_ml_sd, theta_ml_mean = theta_ml_mean))

  print("ML is_dag")
  print(is_dag_ml_mean)

  print("ML correlations within blocks")
  print(lili_correl_same_diag_ml_aggregate)


  ### UNIQUEMENT si TRICHE = TRUE pour la suite

  # li_Lambda <- lapply(1:J, function(j) {
  #   li_Lambda_j <- lili_Lambda_ml_aggregate[[j]]$Lambda_j_mean
  # })
  # mat_Lambda <- as.matrix(Matrix::bdiag(li_Lambda))
  # P_estimate <- t(mat_Lambda) %*% SIGMA_TRUE %*% mat_Lambda
  # P_all_approx <- t(mat_Lambda) %*% Sigma_ML %*% mat_Lambda
  # print(round(as.matrix(P_estimate), 2))
  # print(round(as.matrix(P_all_approx), 2))
}

if (show_TRUE) {
  round_list <- function(x) lapply(x, round, digits = n_signif)

  print("TRUE Lambda")
  print(round_list(li_lambda_TRUE))

  print("TRUE normes Lambda")
  print(round_list(li_lambda_norm_TRUE))

  print("TRUE Phi endo")
  print(round_list(li_R22_TRUE))

  print("TRUE Phi exo")
  print(round_list(li_PHI_TRUE))

  print("TRUE Gamma")
  print(round_list(li_GAMMA_TRUE))

  print("TRUE Beta")
  print(round_list(li_BETA_TRUE))

  print("TRUE Theta")
  print(round(Theta_TRUE, digits = n_signif))

  print("TRUE Psi")
  print(round_list(li_PSI_TRUE))

  print("TRUE R2")
  print(round_list(li_R2_TRUE))

  print("TRUE is_dag")
  print(vec_dag)

  print("TRUE correlations within blocks")
  print(round_list(lapply(seq_along(is_endogenes), function(j) {
    indices <- ((j - 1) * R_true + 1):(j * R_true)
    P_TRUE[indices, indices]
  })))

  print("TRUE latent correlation matrix")
  print(round(P_TRUE, digits = n_signif))
}
