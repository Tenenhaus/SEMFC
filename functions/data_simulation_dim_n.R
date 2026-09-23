######################################
####       data simulation        ####
#### true value of the parameters ####
######################################
library(Matrix)

R_true <- 2

is_endogenes <- c(rep(FALSE, 4), rep(TRUE, 2)) # 4 blocs exogenes et 2 endogenes
n_blocs <- length(is_endogenes)
vec_indicators_per_bloc <- rep(4, n_blocs) # 7 indicateurs par bloc




n_endogenes <- sum(is_endogenes)
n_exogenes <- sum(!is_endogenes)
do_composite <- FALSE
index_var <- rep(0, length(is_endogenes))
count_endo <- 0
count_exo <- 0
for (a in 1:length(is_endogenes)) {
      if (is_endogenes[a]) {
            count_endo <- count_endo + 1
            index_var[a] <- count_endo
      } else {
            count_exo <- count_exo + 1
            index_var[a] <- count_exo
      }
}


index_of_endo <- which(is_endogenes)
index_of_exo <- which(!is_endogenes)


li_BETA_TRUE <- lapply(1:R_true, function(r) {
      if (r == 1) {
            vec_weights <- c(0, 0.3, 0.2, 0)
      } else if (r == 2) {
            vec_weights <- c(0, 0.1, 0.3, 0)
      } else {
            vec_weights <- c(0, 0.5, -0.1, 0)
      }
      Beta_loc <- matrix(vec_weights, n_endogenes, n_endogenes, byrow = TRUE)
      diag(Beta_loc) <- 0
      # print(Beta_loc)
      return(Beta_loc)
})


li_GAMMA_TRUE <- lapply(1:R_true, function(r) {
      if (r == 1) {
            vec_weights <- c(0.2, -0.2, 0, 0, 0, 0, 1, -1)
      } else if (r == 2) {
            vec_weights <- c(0, 0, -0.5, 0.5, -1, -1, 0, 0)
      } else {
            vec_weights <- c(0, 1, 0, 0, 0, 0, 1, 0)
      }
      # print(vec_weights)
      Gamma_loc <- matrix(vec_weights, nrow = n_endogenes, ncol = n_exogenes, byrow = TRUE)
      return(Gamma_loc)
})


li_C_TRUE <- lapply(1:R_true, function(r) {
      Beta_loc <- li_BETA_TRUE[[r]]
      Gamma_loc <- li_GAMMA_TRUE[[r]]
      mask_Beta <- (Beta_loc != 0)
      mask_Gamma <- (Gamma_loc != 0)
      C_loc <- matrix(0, nrow = n_endogenes + n_exogenes, ncol = n_exogenes + n_endogenes)
      for (i in 1:nrow(C_loc)) {
            for (j in 1:ncol(C_loc)) {
                  if (is_endogenes[i] && is_endogenes[j]) {
                        C_loc[i, j] <- mask_Beta[index_var[j], index_var[i]]
                        # print(index_var[j])
                        # print(index_var[i])
                  } else if (!is_endogenes[i] && is_endogenes[j]) {
                        C_loc[i, j] <- mask_Gamma[index_var[j], index_var[i]]
                  } else {
                        C_loc[i, j] <- 0
                  }
            }
      }
      # print(C_loc)
      return(C_loc)
})




li_graphs_TRUE <- lapply(li_C_TRUE, function(C) igraph::graph_from_adjacency_matrix(C))
vec_dag <- sapply(li_graphs_TRUE, igraph::is_dag)
if (any(vec_dag)) {
      print("attention au moins un dag est créé sans recursive structure. For ranks:")
      print(which(vec_dag))
}

# li_R22_TRUE <- lapply(1:R_true, function(r) {
#       if (r == 1) {
#             vec_values <- c(0, 2, 2, 0) * 0.1 / sqrt((r + 1) %/% 2)
#       } else {
#             vec_values <- c(0, 3, 3, 0) * 0.1 / sqrt((r + 1) %/% 2)
#       }
#       # vec_values[1:(length(vec_values) %/% 2)] <- 0.1
#       Phi_endo_r <- matrix(vec_values, n_endogenes, n_endogenes, byrow = TRUE)
#       diag(Phi_endo_r) <- 1
#       Phi_endo_r <- (Phi_endo_r + t(Phi_endo_r)) / 2 # symmetrize
#       return(Phi_endo_r) # cov endo
# })


li_PHI_TRUE <- lapply(1:R_true, function(r) {
      if (r == 1) {
            vec_values <- c(
                  1, 0.2, 0.1, -0.2,
                  0.2, 1, -0.3, 0.2,
                  0.1, -0.3, 1, -0.2,
                  -0.2, 0.2, -0.2, 1
            )
      } else {
            vec_values <- c(
                  1, 0.1, -0.2, 0.3,
                  0.1, 1, 0.3, -0.1,
                  -0.2, 0.3, 1, 0.15,
                  0.3, -0.1, 0.15, 1
            ) + runif(16, min = -0.1, max = 0.1)
      }
      # vec_values[1:(length(vec_values) %/% 2)] <- 0.3
      Phi_exo_r <- matrix(vec_values, n_exogenes, n_exogenes, byrow = TRUE) # Phi_exo
      diag(Phi_exo_r) <- 1
      Phi_exo_r <- (Phi_exo_r + t(Phi_exo_r)) / 2 # symmetrize
      return(Phi_exo_r)
}) # cov exo


# DEBUT AJOUT

# --- Génération "vers l'avant" : B, Gamma, Phi_exo, Psi sont les paramètres
# imposés. Phi_endo (li_R22_TRUE) est calculé, puis tout est standardisé
# pour que les variables latentes endogènes aient variance 1.

li_PSI0_TRUE <- lapply(1:R_true, function(r) {
      # forme raisonnable : disturbances moderement correlees (modele non recursif,
      # la boucle B12/B21 autorise une correlation hors-diagonale sur Psi)
      gr <- li_graphs_TRUE[[r]]
      if (!igraph::is_dag(gr)) {
            Psi0 <- matrix(0.15, n_endogenes, n_endogenes) / sqrt(r)
            diag(Psi0) <- 0.2 / sqrt(r)
      } else {
            Psi0 <- diag(0.2, n_endogenes) / sqrt(r)
      }
      return(Psi0)
})

li_R22_TRUE <- vector("list", R_true) # redevient Phi_endo, mais calcule (pas fixe)
li_PSI_TRUE <- vector("list", R_true)

for (r in 1:R_true) {
      B <- li_BETA_TRUE[[r]]
      G <- li_GAMMA_TRUE[[r]]
      Pe <- li_PHI_TRUE[[r]]
      Psi0 <- li_PSI0_TRUE[[r]]

      IB <- solve(diag(n_endogenes) - B)
      Phi_endo_raw <- IB %*% (G %*% Pe %*% t(G) + Psi0) %*% t(IB)
      Phi_endo_raw <- (Phi_endo_raw + t(Phi_endo_raw)) / 2 # nettoie les arrondis

      d <- sqrt(diag(Phi_endo_raw))
      D <- diag(d, n_endogenes)
      Dinv <- diag(1 / d, n_endogenes)

      # B, Gamma, Psi sont remplaces par leur version standardisee
      li_BETA_TRUE[[r]] <- Dinv %*% B %*% D
      li_GAMMA_TRUE[[r]] <- Dinv %*% G
      li_PSI_TRUE[[r]] <- Dinv %*% Psi0 %*% Dinv

      Phi_endo_std <- Dinv %*% Phi_endo_raw %*% Dinv
      diag(Phi_endo_std) <- 1 # doit deja valoir ~1 ; ceci nettoie les residus numeriques
      li_R22_TRUE[[r]] <- Phi_endo_std
}

# verification: diag(Psi) doit etre dans (0,1), ni proche de 0 ni de 1


# Fin AJOUT


# li_PSI_TRUE <- lapply(1:R_true, function(r) {
#       Psi_loc <- (diag(1, ncol = n_endogenes, nrow = n_endogenes) - li_BETA_TRUE[[r]]) %*% li_R22_TRUE[[r]] %*% t(diag(1, ncol = n_endogenes, nrow = n_endogenes) - li_BETA_TRUE[[r]]) - li_GAMMA_TRUE[[r]] %*% li_PHI_TRUE[[r]] %*% t(li_GAMMA_TRUE[[r]])
#       # print(eigen(Psi_loc, only.values = TRUE, symmetric = TRUE)$values)
#       return(Psi_loc)
# })

li_Phi_full_diag_TRUE <- lapply(1:R_true, function(r) {
      # ATTENTION LES VAR SONT PAS DANS l ORDRE!! C EST ICI dans l ordre exo puis endo: arbitraire
      Phi_r <- rbind(
            cbind(li_PHI_TRUE[[r]], li_PHI_TRUE[[r]] %*% t(li_GAMMA_TRUE[[r]]) %*% t(solve(diag(NROW(li_BETA_TRUE[[r]])) - li_BETA_TRUE[[r]]))),
            cbind(
                  solve(diag(NROW(li_BETA_TRUE[[r]])) - li_BETA_TRUE[[r]]) %*% li_GAMMA_TRUE[[r]] %*% li_PHI_TRUE[[r]],
                  li_R22_TRUE[[r]]
            )
      )
      # print(Phi_r)
      return(Phi_r)
}) # R devient Phi_full

P_TRUE <- matrix(0, nrow = length(index_var) * R_true, ncol = length(index_var) * R_true)

# On remplit les blocs non diagonaux
for (r in 1:R_true) {
      for (i in 1:length(index_var)) {
            for (j in 1:length(index_var)) {
                  if (i <= n_exogenes) {
                        index_i <- index_of_exo[i]
                  } else {
                        index_i <- index_of_endo[i - n_exogenes]
                  }
                  if (j <= n_exogenes) {
                        index_j <- index_of_exo[j]
                  } else {
                        index_j <- index_of_endo[j - n_exogenes]
                  }
                  if (i != j) {
                        P_TRUE[R_true * (index_i - 1) + r, R_true * (index_j - 1) + r] <- li_Phi_full_diag_TRUE[[r]][i, j]
                  } else {
                        P_TRUE[R_true * (index_i - 1) + r, R_true * (index_j - 1) + r] <- 1
                        # P_TRUE[R_true * (i - 1) + r, R_true * (j - 1) + r] <- li_Phi_full_diag_TRUE[[r]][i, j]
                  }
            }
      }
}

# print(P_TRUE)

# print(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values)

least_val_propre <- min(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values) * 0.3
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}


# On complexifie les blocs diagonaux pour ne pas avoir indépendance entre variables latentes d'un meme bloc (tout en gardant la positivité de la matrice)

bloc_perturbation_TRUE <- matrix(least_val_propre, nrow = R_true, ncol = R_true, byrow = TRUE)
diag(bloc_perturbation_TRUE) <- 0
perturbation <- Matrix::bdiag(lapply(1:length(is_endogenes), function(i) bloc_perturbation_TRUE)) # 6 blocs identiques pour lambda


P_TRUE <- P_TRUE + perturbation # On ajoute une petite perturbation pour que la matrice soit positive definie

least_val_propre <- min(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values)
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre apres perturbation de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}

# print(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values)


# Il faut remplir les blocs diagonaux sans perdre la positivité de la matrice.

# Matrice de corr quelconque identique pour chaque var latente


# li_R2_1_TRUE <- lapply(1:R_true, function(r) 1 - li_PSI_TRUE[[r]][1, 1])
# li_R2_2_TRUE <- lapply(1:R_true, function(r) 1 - li_PSI_TRUE[[r]][2, 2])
li_R2_TRUE <- lapply(1:R_true, function(r) {
      vec_R2 <- sapply(1:n_endogenes, function(i) 1 - li_PSI_TRUE[[r]][i, i])
      return(vec_R2)
}) # coeff R^2 for each endogenous latent variable

if (do_composite) {
      SIGMA11_TRUE <- SIGMA22_TRUE <- SIGMA33_TRUE <- SIGMA44_TRUE <- matrix(c(
            1, .3, .4,
            .3, 1, .5,
            .4, .5, 1
      ), 3, 3)

      stop("il faut définir par listes les composites composantes")

      # suite en dimension 1
      w_exo_1 <- w_exo_2 <- rep(1, 3) / drop(sqrt(t(rep(1, 3)) %*% SIGMA11 %*% rep(1, 3)))
      w_exo_3 <- w_exo_4 <- (1:3) / drop(sqrt(t(1:3) %*% SIGMA11 %*% (1:3)))

      l1 <- l2 <- SIGMA11 %*% w_exo_1
      l3 <- l4 <- SIGMA11 %*% w_exo_3

      omega <- list(w_exo_1, w_exo_2, w_exo_3, w_exo_4)
} else {
      li_lambda_TRUE <- lapply(1:length(is_endogenes), function(i) {
            first_lambda_TRUE <- rep(.7, vec_indicators_per_bloc[i]) # vecteur de lambda pour le bloc i
            # construire une matrice de R_true colonnes contenant en premier vcetur bloc_lambda puis que des vecteurs orthogonaux a bloc lambda de norme comparable a bloc_lambda
            norme_cible_TRUE <- sqrt(sum(first_lambda_TRUE^2))
            A <- matrix(rnorm(length(first_lambda_TRUE) * R_true), nrow = length(first_lambda_TRUE), ncol = R_true)
            A[, 1] <- first_lambda_TRUE
            Q <- qr.Q(qr(A))
            Q[, 1] <- first_lambda_TRUE / norme_cible_TRUE

            # Normes des autres colonnes : proches de la premiere,
            # avec une petite perturbation aléatoire
            epsilon <- 0.05
            normes <- abs(norme_cible_TRUE * (1 + rnorm(R_true - 1, 0, epsilon)))

            # Construire la matrice finale
            bloc_lambda_TRUE <- Q
            bloc_lambda_TRUE[, -1] <- sweep(Q[, -1, drop = FALSE], 2, normes, "*")
            # Remplacer la premiere colonne par sa valeur originale
            bloc_lambda_TRUE[, 1] <- first_lambda_TRUE

            for (i in 1:ncol(bloc_lambda_TRUE)) {
                  if (bloc_lambda_TRUE[1, i] < 0) {
                        bloc_lambda_TRUE[, i] <- -bloc_lambda_TRUE[, i]
                  }
            }
            return(bloc_lambda_TRUE)
      })

      Lambda_TRUE <- Matrix::bdiag(li_lambda_TRUE)
}




Sigma_no_perturbation_TRUE <- Lambda_TRUE %*% P_TRUE %*% t(Lambda_TRUE) # covariance imlpied by the model
Theta_TRUE <- abs(0.1 + rnorm(1, mean = 0, sd = 0.05)) * diag(Sigma_no_perturbation_TRUE) # variance residuelle des indicateurs
SIGMA_TRUE <- Sigma_no_perturbation_TRUE + diag(Theta_TRUE) # covariance imlpied by the model avec variance residuelle

# SIGMA11_TRUE <- SIGMA_TRUE[1:3, 1:3]
# SIGMA22_TRUE <- SIGMA_TRUE[4:6, 4:6]
# SIGMA33_TRUE <- SIGMA_TRUE[7:9, 7:9]
# SIGMA44_TRUE <- SIGMA_TRUE[10:12, 10:12]
# SIGMA55_TRUE <- SIGMA_TRUE[13:15, 13:15]
# SIGMA66_TRUE <- SIGMA_TRUE[16:18, 16:18]


# Pour uncomment la suite il faudra augmenter les dimensions!


# true_param_with_S <- c(
#       l1, l2, l3, l4, l5, l6,
#       R[1:4, 1:4][upper.tri(R[1:4, 1:4])],
#       GAMMA[1, 1:2], GAMMA[2, 3:4],
#       BETA[1, 2], BETA[2, 1],
#       R[5, 6],
#       SIGMA11[upper.tri(SIGMA11, diag = TRUE)],
#       SIGMA22[upper.tri(SIGMA22, diag = TRUE)],
#       SIGMA33[upper.tri(SIGMA33, diag = TRUE)],
#       SIGMA44[upper.tri(SIGMA44, diag = TRUE)],
#       1 - l5^2,
#       1 - l6^2
# )
