lvm <- function(R, li_C) {
  # print(R < -1e-5)
  li_gr <- lapply(li_C, function(C) igraph::graph_from_adjacency_matrix(C))
  li_which_exo_endo <- lapply(li_C, function(C) {
    out <- ind_exo_endo(C)
    return(out)
  })
  R_try <- length(li_C)

  li_P_r <- lapply(1:R_try, function(r) {
    J_nb_var <- ncol(li_C[[r]])
    P_r <- matrix(0, nrow = J_nb_var, ncol = J_nb_var)
    for (i in 1:J_nb_var) {
      for (j in 1:J_nb_var) {
        P_r[i, j] <- R[(i - 1) * R_try + r, (j - 1) * R_try + r]
      }
    }
    return(P_r)
  })


  li_BETA <- list()
  li_GAMMA <- list()
  li_PSI <- list()
  li_R2 <- list()
  li_P_induced <- list()

  for (r in 1:R_try) {
    C <- li_C[[r]]
    gr <- li_gr[[r]]
    P_r <- li_P_r[[r]]
    which_exo_endo <- li_which_exo_endo[[r]]


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
        xx[1:l_j, 1:l_j] <- P_r[Ji, H] %*% solve(P_r[H, H]) %*% P_r[H, Ji] # explicatives endo
        xx[1:l_j, (l_j + 1):l] <- P_r[Ji, Hi]
        xx[(l_j + 1):l, 1:l_j] <- P_r[Hi, Ji]
        xx[(l_j + 1):l, (l_j + 1):l] <- P_r[Hi, Hi]

        xy <- c(P_r[Ji, H] %*% solve(P_r[H, H]) %*% P_r[H, J[[i]]], P_r[Hi, J[[i]]])
        bg[[i]] <- solve(xx) %*% xy # on applique la formule 22: regression lineaire
      }

      if (!igraph::is_dag(gr) & config2) {
        xx <- P_r[Ji, H] %*% solve(P_r[H, H]) %*% P_r[H, Ji]
        xy <- P_r[Ji, H] %*% solve(P_r[H, H]) %*% P_r[H, J[[i]]]
        bg[[i]] <- solve(xx) %*% xy
      }

      if (!igraph::is_dag(gr) & config3) {
        xx <- P_r[Hi, Hi]
        xy <- P_r[Hi, J[[i]]]
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
        xx[1:l_j, 1:l_j] <- P_r[Ji, Ji]
        xx[1:l_j, (l_j + 1):l] <- P_r[Ji, Hi]
        xx[(l_j + 1):l, 1:l_j] <- P_r[Hi, Ji]
        xx[(l_j + 1):l, (l_j + 1):l] <- P_r[Hi, Hi]
        xy <- c(P_r[Ji, J[i]], P_r[Hi, J[i]])
        bg[[i]] <- solve(xx) %*% xy
      }

      if (igraph::is_dag(gr) & config2) {
        xx <- P_r[Ji, Ji]
        xy <- P_r[Ji, J[i]]
        bg[[i]] <- solve(xx) %*% xy
      }

      if (igraph::is_dag(gr) & config3) {
        xx <- P_r[Hi, Hi]
        xy <- P_r[Hi, J[i]]
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
      # print(paste("i", i))
      # print(length(b))
      BETA[i, ] <- replace(BETA[i, ], BETA[i, ] == 1, b[[i]][[1]])
      GAMMA[i, ] <- replace(GAMMA[i, ], GAMMA[i, ] == 1, g[[i]][[1]])

      if (igraph::is_dag(gr)) {
        # le terme diagonal est calcule par la formule du residu d'une variable gaussienne conditionnelle a ses explicatives (avec bruit independant): formule OLS
        PSI[i, i] <- 1 - t(bg[[i]]) %*% P_r[c(Ji, Hi), c(Ji, Hi)] %*% bg[[i]]
      }
    }

    if (!igraph::is_dag(gr)) {
      PSI <- (diag(NROW(BETA)) - BETA) %*% P_r[J, J] %*% t(diag(NROW(BETA)) - BETA) -
        GAMMA %*% P_r[H, H] %*% t(GAMMA)
    }


    R2 <- 1 - diag(PSI)

    PI <- solve(diag(NROW(BETA)) - BETA) # (I - B)^-1

    R_LVM <- matrix(0, NCOL(C), NCOL(C))


    if (!igraph::is_dag(gr)) {
      R_LVM[H, H] <- P_r[H, H]
      R_LVM[J, J] <- P_r[J, J]
      R_LVM[H, J] <- P_r[H, H] %*% t(GAMMA) %*% t(PI)
      R_LVM[J, H] <- PI %*% GAMMA %*% P_r[H, H]
    } else {
      R_LVM[H, H] <- P_r[H, H]
      R_LVM[H, J] <- P_r[H, H] %*% t(GAMMA) %*% t(PI)
      R_LVM[J, H] <- PI %*% GAMMA %*% P_r[H, H]
      R_LVM[J, J] <- PI %*% (GAMMA %*% P_r[H, H] %*% t(GAMMA) + PSI) %*% t(PI)
    }

    li_BETA[[r]] <- BETA
    li_GAMMA[[r]] <- GAMMA
    li_PSI[[r]] <- PSI
    li_R2[[r]] <- R2
    li_P_induced[[r]] <- R_LVM


    # print(abs(li_Phi_full_diag_TRUE[[r]] - R_LVM) > 1e-8)
  }

  P_induced <- matrix(0, nrow = NCOL(R), ncol = NCOL(R))

  for (i in 1:NCOL(li_P_induced[[r]])) {
    for (j in 1:NCOL(li_P_induced[[r]])) {
      if (i != j) {
        for (r in 1:R_try) {
          P_induced[(i - 1) * R_try + r, (j - 1) * R_try + r] <- li_P_induced[[r]][i, j]
        }
      } else {
        P_induced[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- R[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)]
      }
    }
  }

  # print(abs(P_induced - R) > 1e-8)

  ecart_created <- d_LS(P_induced, R)
  ecart_relat <- ecart_created / norm(R, "F")
  print(paste("Ration modif P:", round(ecart_relat, 4)))


  return(list(
    li_gr = li_gr,
    li_BETA = li_BETA,
    li_GAMMA = li_GAMMA,
    li_PSI = li_PSI,
    li_R2 = li_R2,
    li_P_induced = li_P_induced,
    P_induced = P_induced
  ))
}
