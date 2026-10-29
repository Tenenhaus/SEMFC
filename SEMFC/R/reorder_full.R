reorder_full <- function(li_S, li_Lambda_star, li_Lambda, li_vec_norm, li_C_jj_hat, R_try, J, ortho) {
    # lister les permutations de 1:R_try
    mat_permutations <- gtools::permutations(n = R_try, r = R_try, v = 1:R_try)
    # print(mat_permutations)
    n_permut <- nrow(mat_permutations)
    best_score <- Inf
    cartesian_product <- expand.grid(rep(list(1:n_permut), J - 1))
    for (i in 1:nrow(cartesian_product)) {
        score <- 0
        for (j_1 in 1:J) {
            for (j_2 in 1:J) {
                if (j_1 != j_2) {
                    if (j_1 != 1) {
                        choice_j1 <- cartesian_product[i, j_1 - 1]
                        permut_j1 <- mat_permutations[choice_j1, ]
                    } else {
                        permut_j1 <- 1:R_try
                    }
                    if (j_2 != 1) {
                        choice_j2 <- cartesian_product[i, j_2 - 1]
                        permut_j2 <- mat_permutations[choice_j2, ]
                    } else {
                        permut_j2 <- 1:R_try
                    }
                    if (ortho) {
                        P_cross <- diag(1 / li_vec_norm[[j_1]][permut_j1], nrow = R_try, ncol = R_try) %*% t(li_Lambda_star[[j_1]][, permut_j1]) %*% li_S[[j_1]][[j_2]] %*% li_Lambda_star[[j_2]][, permut_j2] %*% diag(1 / li_vec_norm[[j_2]][permut_j2], nrow = R_try, ncol = R_try)
                    } else {
                        Lambda_i <- li_Lambda[[j_1]][, permut_j1]
                        Lambda_j <- li_Lambda[[j_2]][, permut_j2]
                        min_eigen_i <- min(eigen(t(Lambda_i) %*% Lambda_i)$values)
                        min_eigen_j <- min(eigen(t(Lambda_j) %*% Lambda_j)$values)
                        if (min_eigen_i < 1e-8 | min_eigen_j < 1e-8) {
                            print(paste0("Warning: near singular matrix in reorder_full for j_1 = ", j_1, " and j_2 = ", j_2, "values: ", min_eigen_i, " and ", min_eigen_j))
                        }
                        pseudi_i <- solve(t(Lambda_i) %*% Lambda_i) %*% t(Lambda_i)
                        pseudi_j <- solve(t(Lambda_j) %*% Lambda_j) %*% t(Lambda_j)
                        diag_i <- diag(1 / li_vec_norm[[j_1]][permut_j1], nrow = R_try, ncol = R_try)
                        diag_j <- diag(1 / li_vec_norm[[j_2]][permut_j2], nrow = R_try, ncol = R_try)
                        P_cross <- diag_i %*% pseudi_i %*% li_S[[j_1]][[j_2]] %*% t(pseudi_j) %*% diag_j
                    }
                    diag(P_cross) <- 0
                    score <- score + sum(abs(P_cross))
                }
            }
            if (score < best_score) {
                best_score <- score
                best_vec <- as.numeric(cartesian_product[i, ])
            }
        }
    }
    for (j in 2:J) {
        choice_j <- best_vec[j - 1]
        permut_j <- mat_permutations[choice_j, ]
        li_Lambda_star[[j]] <- li_Lambda_star[[j]][, permut_j]
        li_Lambda[[j]] <- li_Lambda[[j]][, permut_j]
        li_vec_norm[[j]] <- li_vec_norm[[j]][permut_j]
        li_C_jj_hat[[j]] <- li_C_jj_hat[[j]][permut_j, permut_j]
    }

    return(list(
        li_Lambda_star = li_Lambda_star,
        li_Lambda = li_Lambda,
        li_vec_norm = li_vec_norm,
        li_C_jj_hat = li_C_jj_hat
    ))
}
