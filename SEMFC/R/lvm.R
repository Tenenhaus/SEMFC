lvm <- function(R, C) {
  gr <- igraph::graph_from_adjacency_matrix(C)
  which_exo_endo <- ind_exo_endo(C)
  rownames(C) <- colnames(C) <- rownames(R)
  xx <- list()
  xy <- list()
  bg <- list()
  b <- list()
  g <- list()
  tC <- t(C)
  H <- which_exo_endo$ind_exo
  J <- which_exo_endo$ind_endo

  BETA <- as.matrix(tC[J, J, drop = F])
  GAMMA <- as.matrix(tC[J, H, drop = F])
  PSI <- matrix(0, NCOL(BETA), NCOL(BETA))

  for (i in 1:length(which_exo_endo$Ji)) {
    Ji <- which_exo_endo$Ji[[i]] # explicatives endo
    Hi <- which_exo_endo$Hi[[i]] # explicatives exo
    l_j <- length(Ji)
    l_h <- length(Hi)
    l <- length(c(Hi, Ji))

    config1 <- any(Ji != 0) & any(Hi != 0) # Hi ou Ji = 0 veut dire qu'aucune variable exo resp endo n'explique la i-eme endogene
    config2 <- any(Ji != 0) & any(Hi == 0)
    config3 <- any(Ji == 0) & any(Hi != 0)
    config4 <- any(Ji == 0) & any(Hi == 0)

    if (!igraph::is_dag(gr) & config1) {
      # Calcul de la cov des variables explicatives pour la i-eme endogene
      xx <- matrix(NA, l, l)
      xx[1:l_j, 1:l_j] <- R[Ji, H] %*% solve(R[H, H]) %*% R[H, Ji] # explicatives endo
      xx[1:l_j, (l_j + 1):l] <- R[Ji, Hi]
      xx[(l_j + 1):l, 1:l_j] <- R[Hi, Ji]
      xx[(l_j + 1):l, (l_j + 1):l] <- R[Hi, Hi]

      xy <- c(R[Ji, H] %*% solve(R[H, H]) %*% R[H, J[[i]]], R[Hi, J[[i]]])
      bg[[i]] <- solve(xx) %*% xy # on applique la formule 22: regression lineaire
    }

    if (!igraph::is_dag(gr) & config2) {
      xx <- R[Ji, H] %*% solve(R[H, H]) %*% R[H, Ji]
      xy <- R[Ji, H] %*% solve(R[H, H]) %*% R[H, J[[i]]]
      bg[[i]] <- solve(xx) %*% xy
    }

    if (!igraph::is_dag(gr) & config3) {
      xx <- R[Hi, Hi]
      xy <- R[Hi, J[[i]]]
      bg[[i]] <- solve(xx) %*% xy
    }

    if (!igraph::is_dag(gr) & config4) {
      print(paste("the ", i, "th structural model is not properly specified",
        sep = ""
      ))
      break
    }

    # Si le graphe est dag, on decide que Psi sera diagonale et on applique les formules d'un modele recursif

    if (igraph::is_dag(gr) & config1) {
      xx <- matrix(NA, l, l)
      xx[1:l_j, 1:l_j] <- R[Ji, Ji]
      xx[1:l_j, (l_j + 1):l] <- R[Ji, Hi]
      xx[(l_j + 1):l, 1:l_j] <- R[Hi, Ji]
      xx[(l_j + 1):l, (l_j + 1):l] <- R[Hi, Hi]
      xy <- c(R[Ji, J[i]], R[Hi, J[i]])
      bg[[i]] <- solve(xx) %*% xy
    }

    if (igraph::is_dag(gr) & config2) {
      xx <- R[Ji, Ji]
      xy <- R[Ji, J[i]]
      bg[[i]] <- solve(xx) %*% xy
    }

    if (igraph::is_dag(gr) & config3) {
      xx <- R[Hi, Hi]
      xy <- R[Hi, J[i]]
      bg[[i]] <- solve(xx) %*% xy
    }


    if (igraph::is_dag(gr) & config4) {
      print(paste("the ", i, "th structural model is not properly specified",
        sep = ""
      ))
      break
    }

    ifelse(Ji == 0,
      yes = {
        b[[i]] <- 0
      },
      no = {
        b[[i]] <- list(bg[[i]][1:length(Ji)])
      }
    )

    ifelse(Hi == 0,
      yes = {
        g[[i]] <- 0
      },
      no = {
        ifelse(Ji == 0,
          yes = {
            g[[i]] <- list(bg[[i]])
          },
          no = {
            g[[i]] <- list(bg[[i]][seq(length(bg[[i]]))[-seq(length(Ji))]])
          }
        )
      }
    )
    # remplir les lignes qui viennent d'etre calculees dans les matrices B et GAMMA
    BETA[i, ] <- replace(BETA[i, ], BETA[i, ] == 1, b[[i]][[1]])
    GAMMA[i, ] <- replace(GAMMA[i, ], GAMMA[i, ] == 1, g[[i]][[1]])

    if (igraph::is_dag(gr)) {
      # le terme diagonal est calcule par la formule du residu d'une variable gaussienne conditionnelle a ses explicatives (avec bruit independant): formule OLS
      PSI[i, i] <- 1 - t(bg[[i]]) %*% R[c(Ji, Hi), c(Ji, Hi)] %*% bg[[i]]
    }
  }

  if (!igraph::is_dag(gr)) {
    PSI <- (diag(NROW(BETA)) - BETA) %*% R[J, J] %*% t(diag(NROW(BETA)) - BETA) -
      GAMMA %*% R[H, H] %*% t(GAMMA)
  }

  R2 <- 1 - diag(PSI)

  PI <- solve(diag(NROW(BETA)) - BETA) # (I - B)^-1

  R_LVM <- matrix(0, NCOL(C), NCOL(C))

  if (!igraph::is_dag(gr)) {
    R_LVM[H, H] <- R[H, H]
    R_LVM[J, J] <- R[J, J]
    R_LVM[H, J] <- R[H, H] %*% t(GAMMA) %*% t(PI)
    R_LVM[J, H] <- PI %*% GAMMA %*% R[H, H]
  } else {
    R_LVM[H, H] <- R[H, H]
    R_LVM[H, J] <- R[H, H] %*% t(GAMMA) %*% t(PI)
    R_LVM[J, H] <- PI %*% GAMMA %*% R[H, H]
    R_LVM[J, J] <- PI %*% (GAMMA %*% R[H, H] %*% t(GAMMA) + PSI) %*% t(PI)
  }


  return(list(
    gr = gr,
    BETA = BETA, GAMMA = GAMMA,
    PSI = PSI,
    R2 = R2,
    R_LVM = R_LVM
  ))
}
