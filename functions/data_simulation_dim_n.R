######################################
####       data simulation        ####
#### true value of the parameters ####
######################################
library(Matrix)

R_true <- 2
do_composite <- FALSE

li_BETA <- lapply(1:R_true, function(r) {
      (0.95)^(r - 1) * matrix(c(
            0, 0.25,
            0.5, 0
      ), 2, 2, byrow = TRUE)
})


li_GAMMA <- lapply(1:R_true, function(r) {
      (0.95)^(r - 1) * matrix(c(
            -0.30, 0.5, 0, 0,
            0, 0, 0.5, 0.25
      ), 2, 4, byrow = TRUE)
})



li_R22 <- lapply(1:R_true, function(r) {
      matrix(c(1, sqrt(1 / 2), sqrt(1 / 2), 1), 2, 2) # Phi endo
})

li_PHI <- lapply(1:R_true, function(r) {
      matrix(c(
            1, .5, .5, .5,
            .5, 1, .5, .5,
            .5, .5, 1, .5,
            .5, .5, .5, 1
      ), 4, 4, byrow = TRUE) # Phi_exo
})

li_PSI <- lapply(1:R_true, function(r) {
      (diag(2) - li_BETA[[r]]) %*% li_R22[[r]] %*% t(diag(2) - li_BETA[[r]]) - li_GAMMA[[r]] %*% li_PHI[[r]] %*% t(li_GAMMA[[r]])
})

li_Phi_full_diag <- lapply(1:R_true, function(r) {
      Phi_r <- rbind(
            cbind(li_PHI[[r]], li_PHI[[r]] %*% t(li_GAMMA[[r]]) %*% t(solve(diag(NROW(li_BETA[[r]])) - li_BETA[[r]]))),
            cbind(
                  solve(diag(NROW(li_BETA[[r]])) - li_BETA[[r]]) %*% li_GAMMA[[r]] %*% li_PHI[[r]],
                  li_R22[[r]]
            )
      )
      return(Phi_r)
}) # R devient Phi_full

P <- matrix(0, nrow = 6 * R_true, ncol = 6 * R_true)

# On remplit les blocs non diagonaux
for (r in 1:R_true) {
      for (i in 1:6) {
            for (j in 1:6) {
                  if (i != j) {
                        P[R_true * (i - 1) + r, R_true * (j - 1) + r] <- li_Phi_full_diag[[r]][i, j]
                  } else {
                        P[R_true * (i - 1) + r, R_true * (j - 1) + r] <- 1
                        # P[R_true * (i - 1) + r, R_true * (j - 1) + r] <- li_Phi_full_diag[[r]][i, j]
                  }
            }
      }
}

least_val_propre <- min(eigen(P, only.values = TRUE, symmetric = TRUE)$values) * 0.7
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}

# On complexifie les blocs diagonaux pour ne pas avoir indépendance entre variables latentes d'un meme bloc (tout en gardant la positivité de la matrice)

bloc_perturbation <- matrix(least_val_propre, nrow = R_true, ncol = R_true, byrow = TRUE)
diag(bloc_perturbation) <- 0
perturbation <- Matrix::bdiag(lapply(1:6, function(i) bloc_perturbation))


P <- P + perturbation # On ajoute une petite perturbation pour que la matrice soit positive definie

least_val_propre <- min(eigen(P, only.values = TRUE, symmetric = TRUE)$values) * 0.7
if (least_val_propre <= 0) {
      print(paste("La plus petite valeur propre apres perturbation de la matrice de correlation est :", least_val_propre))
      stop("La matrice de correlation entre les variables latentes n'est pas positive definie. Il faut revoir la definition des parametres du modele.")
}


# Il faut remplir les blocs diagonaux sans perdre la positivité de la matrice.

# Matrice de corr quelconque identique pour chaque var latente




li_R2_1 <- lapply(1:R_true, function(r) 1 - li_PSI[[r]][1, 1])
li_R2_2 <- lapply(1:R_true, function(r) 1 - li_PSI[[r]][2, 2])
li_R2 <- lapply(1:R_true, function(r) c(li_R2_1[[r]], li_R2_2[[r]])) # coeff R^2 for each endogenous latent variable

if (do_composite) {
      SIGMA11 <- SIGMA22 <- SIGMA33 <- SIGMA44 <- matrix(c(
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
      first_lambda <- rep(.7, 3)
      # construire une matrice de R_true colonnes contenant en premier vcetur bloc_lambda puis que des vecteurs orthogonaux a bloc lambda de norme comparable a bloc_lambda
      norme_cible <- sqrt(sum(first_lambda^2))
      A <- matrix(rnorm(length(first_lambda) * R_true), nrow = length(first_lambda), ncol = R_true)
      A[, 1] <- first_lambda
      Q <- qr.Q(qr(A))
      Q[, 1] <- first_lambda / norme_cible

      # Normes des autres colonnes : proches de la premiere,
      # avec une petite perturbation aléatoire
      epsilon <- 0.05
      normes <- abs(norme_cible * (1 + rnorm(R_true - 1, 0, epsilon)))

      # Construire la matrice finale
      bloc_lambda <- Q
      bloc_lambda[, -1] <- sweep(Q[, -1, drop = FALSE], 2, normes, "*")
      # Remplacer la premiere colonne par sa valeur originale
      bloc_lambda[, 1] <- first_lambda

      Lambda <- Matrix::bdiag(lapply(1:6, function(i) bloc_lambda)) # 6 blocs identiques pour lambda
}


Sigma_no_perturbation <- Lambda %*% P %*% t(Lambda) # covariance imlpied by the model
Theta <- abs(0.1 + rnorm(1, mean = 0, sd = 0.01)) * diag(Sigma_no_perturbation) # variance residuelle des indicateurs
SIGMA <- Sigma_no_perturbation + diag(Theta) # covariance imlpied by the model avec variance residuelle

SIGMA11 <- SIGMA[1:3, 1:3]
SIGMA22 <- SIGMA[4:6, 4:6]
SIGMA33 <- SIGMA[7:9, 7:9]
SIGMA44 <- SIGMA[10:12, 10:12]
SIGMA55 <- SIGMA[13:15, 13:15]
SIGMA66 <- SIGMA[16:18, 16:18]

li_lambda <- lapply(1:6, function(i) bloc_lambda)


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
