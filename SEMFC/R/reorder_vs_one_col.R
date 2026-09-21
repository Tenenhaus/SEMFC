reorder_vs_one_col <- function(li_S, li_Lambda_star, li_Lambda, li_vec_norm, li_C_jj_hat, R_try, J) {
    for (j in 2:J) {
        Sigma_1j <- li_S[[1]][[j]]
        svd_fit <- svd(Sigma_1j, nu = R_try, nv = R_try)

        left <- svd_fit$u
        right <- svd_fit$v
        Lambda_1_star <- li_Lambda_star[[1]]
        # retrouver les colonnes colineaires entre elles
        mat_ps <- t(t(left) %*% Lambda_1_star)
        where_one <- abs(mat_ps) > 0.5
        # print(mat_ps)
        sum_per_col <- apply(where_one, 1, sum)
        if (any(sum_per_col != 1)) {
            print(paste("Pb reordering 1 bloc", j))
            print("values svd")
            print(svd_fit$d)
            print("mat ps")
            print(mat_ps)
        }


        # vec_reorder_1 <- apply(mat_ps, 1, function(x) which.max(abs(x))[1])
        vec_reorder_1 <- as.integer(clue::solve_LSAP(abs(mat_ps), maximum = TRUE)) # deja un vecteur d'entiers, plus besoin de as.integer si vous castez apres
        good_right <- right[, vec_reorder_1]
        # print(left[, vec_reorder_1])
        # Lambda_j dans le bon ordre mais un peu mal estimed

        mat_ps <- t(t(li_Lambda_star[[j]]) %*% good_right) # POURQUOI CELA MARCHE?
        where_one <- abs(mat_ps) > 0.5
        sum_per_col <- apply(where_one, 1, sum)
        # print(mat_ps)
        if (any(sum_per_col != 1)) {
            print(paste("Pb reordering 2 bloc", j))
            print("values svd")
            print(svd_fit$d)
            print("mat ps")
            print(mat_ps)
        }

        # vec_reorder_2  <- apply(mat_ps, 1, function(x) which.max(abs(x))[1])
        vec_reorder_2 <- as.integer(clue::solve_LSAP(abs(mat_ps), maximum = TRUE))

        if (any(duplicated(vec_reorder_1)) | any(duplicated(vec_reorder_2))) {
            stop("Duplicated elements in vec_reorder_1 or vec_reorder_2")
        }
        li_Lambda_star[[j]] <- li_Lambda_star[[j]][, vec_reorder_2]
        li_Lambda[[j]] <- li_Lambda[[j]][, vec_reorder_2]
        li_vec_norm[[j]] <- li_vec_norm[[j]][vec_reorder_2]
        li_C_jj_hat[[j]] <- li_C_jj_hat[[j]][vec_reorder_2, vec_reorder_2]

        # print(round(t(li_lambda_TRUE[[j]]) %*% li_Lambda[[j]]), 4)
    }
    return(list(
        li_Lambda_star = li_Lambda_star,
        li_Lambda = li_Lambda,
        li_vec_norm = li_vec_norm,
        li_C_jj_hat = li_C_jj_hat
    ))
}
