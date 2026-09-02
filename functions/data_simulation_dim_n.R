######################################
####       data simulation        ####
#### true value of the parameters ####
######################################
library(Matrix)

R_true <- 2
do_composite <- FALSE

li_BETA_TRUE <- lapply(1:R_true, function(r) {
      matrix(c(
            0, 0.25,
            0.5, 0
      ), 2, 2, byrow = TRUE)
})


li_GAMMA_TRUE <- lapply(1:R_true, function(r) {
      matrix(c(
            -0.30, 0.5, 0, 0,
            0, 0, 0.5, 0.25
      ), 2, 4, byrow = TRUE)
})


li_R22_TRUE <- lapply(1:R_true, function(r) {
      matrix(c(1, sqrt(1 / 2), sqrt(1 / 2), 1), 2, 2) # Phi endo
})

li_PHI_TRUE <- lapply(1:R_true, function(r) {
      matrix(c(
            1, .5, .5, .5,
            .5, 1, .5, .5,
            .5, .5, 1, .5,
            .5, .5, .5, 1
      ), 4, 4, byrow = TRUE) # Phi_exo
})

li_PSI_TRUE <- lapply(1:R_true, function(r) {
      Psi_loc <- (diag(2) - li_BETA_TRUE[[r]]) %*% li_R22_TRUE[[r]] %*% t(diag(2) - li_BETA_TRUE[[r]]) - li_GAMMA_TRUE[[r]] %*% li_PHI_TRUE[[r]] %*% t(li_GAMMA_TRUE[[r]])
      return(Psi_loc)
})

li_Phi_full_diag_TRUE <- lapply(1:R_true, function(r) {
      Phi_r <- rbind(
            cbind(li_PHI_TRUE[[r]], li_PHI_TRUE[[r]] %*% t(li_GAMMA_TRUE[[r]]) %*% t(solve(diag(NROW(li_BETA_TRUE[[r]])) - li_BETA_TRUE[[r]]))),
            cbind(
                  solve(diag(NROW(li_BETA_TRUE[[r]])) - li_BETA_TRUE[[r]]) %*% li_GAMMA_TRUE[[r]] %*% li_PHI_TRUE[[r]],
                  li_R22_TRUE[[r]]
            )
      )
      return(Phi_r)
}) # R devient Phi_full

P_TRUE <- matrix(0, nrow = 6 * R_true, ncol = 6 * R_true)

# On remplit les blocs non diagonaux
for (r in 1:R_true) {
      for (i in 1:6) {
            for (j in 1:6) {
                  if (i != j) {
                        P_TRUE[R_true * (i - 1) + r, R_true * (j - 1) + r] <- li_Phi_full_diag_TRUE[[r]][i, j]
                  } else {
                        P_TRUE[R_true * (i - 1) + r, R_true * (j - 1) + r] <- 1
                        # P_TRUE[R_true * (i - 1) + r, R_true * (j - 1) + r] <- li_Phi_full_diag_TRUE[[r]][i, j]
                  }
            }
      }
}

least_val_propre <- min(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values) * 0.7
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}

# On complexifie les blocs diagonaux pour ne pas avoir indépendance entre variables latentes d'un meme bloc (tout en gardant la positivité de la matrice)

bloc_perturbation_TRUE <- matrix(least_val_propre, nrow = R_true, ncol = R_true, byrow = TRUE)
diag(bloc_perturbation_TRUE) <- 0
perturbation <- Matrix::bdiag(lapply(1:6, function(i) bloc_perturbation_TRUE)) # 6 blocs identiques pour lambda


P_TRUE <- P_TRUE + perturbation # On ajoute une petite perturbation pour que la matrice soit positive definie

least_val_propre <- min(eigen(P_TRUE, only.values = TRUE, symmetric = TRUE)$values) * 0.7
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre apres perturbation de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}


# Il faut remplir les blocs diagonaux sans perdre la positivité de la matrice.

# Matrice de corr quelconque identique pour chaque var latente


li_R2_1_TRUE <- lapply(1:R_true, function(r) 1 - li_PSI_TRUE[[r]][1, 1])
li_R2_2_TRUE <- lapply(1:R_true, function(r) 1 - li_PSI_TRUE[[r]][2, 2])
li_R2_TRUE <- lapply(1:R_true, function(r) c(li_R2_1_TRUE[[r]], li_R2_2_TRUE[[r]])) # coeff R^2 for each endogenous latent variable

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
      first_lambda_TRUE <- rep(.7, 3)
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

      Lambda_TRUE <- Matrix::bdiag(lapply(1:6, function(i) bloc_lambda_TRUE)) # 6 blocs identiques pour lambda
}



Sigma_no_perturbation_TRUE <- Lambda_TRUE %*% P_TRUE %*% t(Lambda_TRUE) # covariance imlpied by the model
Theta_TRUE <- abs(0.1 + rnorm(1, mean = 0, sd = 0.01)) * diag(Sigma_no_perturbation_TRUE) # variance residuelle des indicateurs
SIGMA_TRUE <- Sigma_no_perturbation_TRUE + diag(Theta_TRUE) # covariance imlpied by the model avec variance residuelle

SIGMA11_TRUE <- SIGMA_TRUE[1:3, 1:3]
SIGMA22_TRUE <- SIGMA_TRUE[4:6, 4:6]
SIGMA33_TRUE <- SIGMA_TRUE[7:9, 7:9]
SIGMA44_TRUE <- SIGMA_TRUE[10:12, 10:12]
SIGMA55_TRUE <- SIGMA_TRUE[13:15, 13:15]
SIGMA66_TRUE <- SIGMA_TRUE[16:18, 16:18]

li_lambda_TRUE <- lapply(1:6, function(i) bloc_lambda_TRUE)


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
