#' Structural equation models with factors and composites with svd-SEM
#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{X_1, X_2, ..., X_J}.
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances.
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param mode a vector of lenght J indication the mode for each block formative or reflective.
#' @return \item{a}{A list of \eqn{J} elements. Each element of \eqn{a} is a matrix that contains the first singular vector for each block.}
#' @references Tenenhaus A., Tenenhaus M. and Dijkstra (2025),
#' @title Structural Equation Modeling with Factors and Composites (svdSEM)
#' @examples
#' #############
#' # Example 1 #
#' #############
#' ECSI <- read.table("mobil.txt") / 10
#' A <- list(
#'     CUSTOMER_E = ECSI[, c("CUEX1", "CUEX2", "CUEX3")],
#'     PERC_QUAL = ECSI[, c(
#'         "PERQ1", "PERQ2", "PERQ3", "PERQ4",
#'         "PERQ5", "PERQ6", "PERQ7"
#'     )],
#'     PERC_VALUE = ECSI[, c("PERV1", "PERV2")],
#'     CUSTOMER_S = ECSI[, c("CUSA1", "CUSA2", "CUSA3")],
#'     CUSTOMER_L = ECSI[, c("CUSL1", "CUSL2", "CUSL3")]
#' )
#'
#' C <- matrix(c(
#'     0, 0, 0, 0, 0,
#'     1, 0, 0, 0, 0,
#'     1, 1, 0, 0, 0,
#'     1, 1, 1, 0, 0,
#'     0, 0, 0, 1, 0
#' ), 5, 5, byrow = FALSE)
#'
#' fit <- svdSEM(A, C,
#'     scale = FALSE,v2
#'     mode = rep("reflective", length(A)),
#'     bias = FALSE
#' )
#'
#' @export svdSEM_multi
svdSEM_multi <- function(A, li_C, scale = TRUE,
                         mode = rep("reflective", length(A)),
                         bias = FALSE) {
    if (any(mode == "formative")) {
        stop("svdSEM_multi() does not support yet formative blocks. Please use svdSEM() instead.")
    }
    pjs <- sapply(A, NCOL)
    nb_row <- NROW(A[[1]])

    li_Lambda <- list()
    nb_ind <- NROW(A[[1]])
    J <- length(A)
    R_try <- length(li_C)

    cumsum_pjs <- cumsum(c(0, pjs))

    if(!triche){
        li_S <- lapply(A, function(x){
                lapply(A, function(y) cov2(x, y, bias = bias)) 
            })
        } else{

        

        ### On triche: SIGMA connu

        li_S <- lapply(1:J, function(j) {
            lapply(1:J, function(k) {
                S_jk <- as.matrix(SIGMA_TRUE[(cumsum_pjs[j] + 1):(cumsum_pjs[j + 1]), (cumsum_pjs[k] + 1):(cumsum_pjs[k + 1])])
            })
        })
    }


    S_full <- matrix(0, nrow = sum(pjs), ncol = sum(pjs))
    for (i in 1:J) {
        for (j in 1:J) {
            S_full[(cumsum_pjs[i] + 1):(cumsum_pjs[i + 1]), (cumsum_pjs[j] + 1):(cumsum_pjs[j + 1])] <- li_S[[i]][[j]]
        }
    }

    # print("init diff")
    # print(round(as.matrix(abs(SIGMA_TRUE)), 3))
    # print(round((as.matrix(abs(S_full - SIGMA_TRUE)/abs(SIGMA_TRUE + 1e-16))), 3))

    #-------------------------------------------------------
    blocks <- A
    if (scale) {
        A <- lapply(A, function(x) scale2(x, bias = bias))
    } else {
        A <- lapply(A, function(x) scale2(x, bias = bias, scale = FALSE))
    }

    

    if(ortho_solve_SVD){
        # Extract first singular vector for each block. Nouvelle methode: celle du papier
        li_Lambda_star <- lapply(
            1:J,
            function(j) {
                li_sigma_jk <- lapply(1:J, function(k) {
                    Sigma_jk <- li_S[[j]][[k]]
                    if (j == k) {
                        Sigma_jk <- matrix(0, nrow = pjs[j], ncol = pjs[k])
                    }
                    prod <- Sigma_jk %*% t(Sigma_jk)
                    return(prod)
                })
                mat_sum <- Reduce("+", li_sigma_jk)
                svd.fit <- svd(mat_sum, nu = R_try, nv = 1)
                Lambda_star <- svd.fit$u
                #print(svd.fit$d)
                return(Lambda_star)

                # svd(t(A[[x]]) %*% Reduce("cbind", A[-x]),
                #     nu = 1, nv = 1
                # )$u
            }
        )


        # check for sign inversion
        li_Lambda_star <- lapply(li_Lambda_star, function(Lambda_star) {
            for (r in 1:R_try) {
                if (Lambda_star[1, r] < 0 ) {
                    Lambda_star[, r] <- -Lambda_star[, r]
                }
            }
            return(Lambda_star)
        })
        names(li_Lambda_star) <- names(A)
    } else{
        #-------------------------------------------------------
        # 1. Reduction to R x R matrices
        #-------------------------------------------------------

        # Singular vectors associated with the R largest singular values
        # of the concatenation of all cross-block covariance matrices.

        li_U <- lapply(1:J, function(i) {

            C_i <- do.call(cbind, lapply(
                setdiff(1:J, i),
                function(j) li_S[[i]][[j]]
            ))

            svd.fit <- svd(C_i, nu = R_try, nv = 0)

            svd.fit$u[, 1:R_try, drop = FALSE]
        })


        # Reduced cross-block covariance matrices:
        # W_ij = U_i' Sigma_ij U_j

        li_W <- lapply(1:J, function(i) {
            lapply(1:J, function(j) {

                if (i == j) {
                    return(matrix(0, nrow = R_try, ncol = R_try))
                }

                t(li_U[[i]]) %*% li_S[[i]][[j]] %*% li_U[[j]]
            })
        })


        # Ratio between the smallest and largest singular values.
        # A small ratio indicates a poorly conditioned matrix.

        svd_ratio <- function(X) {

            d <- svd(X, nu = 0, nv = 0)$d

            if (length(d) == 0 || !all(is.finite(d)) || max(d) == 0) {
                return(0)
            }

            min(d) / max(d)
        }


        #-------------------------------------------------------
        # 2. Estimate the directions block by block
        #-------------------------------------------------------

        # Number of random (u, v) configurations tested per pivot
        n_uv <- 10

        # Number of candidate linear combinations tested
        n_alpha <- 100

        # Numerical thresholds
        tol_W <- 1e-6
        tol_M <- 1e-6
        tol_imag <- 1e-6


        li_Lambda_star <- lapply(1:J, function(a) {

            # Store the candidate matrices F and their quality scores.
            # Each F is normalized before aggregation.

            li_F <- list()
            quality <- numeric(0)


            #---------------------------------------------------
            # 2.1 Construct candidate matrices F for target a
            #---------------------------------------------------

            for (b in setdiff(1:J, a)) {

                js <- setdiff(1:J, c(a, b))

                # M_j = W_aj W_bj^{-1}
                # Skip this pivot if one of the W_bj is
                # numerically too close to being singular.

                li_M <- vector("list", length(js))
                pivot_ok <- TRUE

                for (k in seq_along(js)) {

                    j <- js[k]

                    W_bj <- li_W[[b]][[j]]

                    if (svd_ratio(W_bj) < tol_W) {
                        pivot_ok <- FALSE
                        break
                    }

                    inv_W_bj <- tryCatch(
                        solve(W_bj),
                        error = function(e) NULL
                    )

                    if (is.null(inv_W_bj)) {
                        pivot_ok <- FALSE
                        break
                    }

                    li_M[[k]] <- li_W[[a]][[j]] %*% inv_W_bj
                }

                if (!pivot_ok) {
                    next
                }


                # Generate random Gaussian vectors u and v.
                # Normalize them to avoid arbitrary scale effects.

                for (iter in 1:n_uv) {

                    u <- rnorm(length(js))
                    v <- rnorm(length(js))

                    u <- u / sqrt(sum(u^2))
                    v <- v / sqrt(sum(v^2))


                    # M(u) = sum_j u_j M_j
                    # M(v) = sum_j v_j M_j

                    M_u <- matrix(0, nrow = R_try, ncol = R_try)
                    M_v <- matrix(0, nrow = R_try, ncol = R_try)

                    for (k in seq_along(js)) {
                        M_u <- M_u + u[k] * li_M[[k]]
                        M_v <- M_v + v[k] * li_M[[k]]
                    }


                    # Reject poorly conditioned M(v).

                    cond_v <- svd_ratio(M_v)

                    if (!is.finite(cond_v) || cond_v < tol_M) {
                        next
                    }


                    # F = M(u) M(v)^{-1}

                    F <- tryCatch(
                        M_u %*% solve(M_v),
                        error = function(e) NULL
                    )

                    if (is.null(F) || !all(is.finite(F))) {
                        next
                    }

                    norm_F <- sqrt(sum(F^2))

                    if (!is.finite(norm_F) || norm_F < 1e-12) {
                        next
                    }


                    # Normalize F and retain the conditioning score.

                    li_F[[length(li_F) + 1L]] <- F / norm_F
                    quality <- c(quality, cond_v)
                }
            }


            # Stop if no usable candidate was found.

            K <- length(li_F)

            if (K == 0) {
                stop(
                    paste(
                        "No usable candidate matrix for block", a,
                        "- check the covariance matrices and numerical thresholds."
                    )
                )
            }


            #---------------------------------------------------
            # 2.2 Aggregate the candidate matrices
            #---------------------------------------------------

            # Better-conditioned configurations receive larger
            # weights. We then test random combinations of these
            # matrices and retain the one with the largest minimum
            # relative eigenvalue gap.

            quality <- pmax(quality, 1e-8)

            best_score <- -Inf
            best_eig <- NULL
            best_imag <- Inf


            for (iter in 1:n_alpha) {

                if (iter == 1) {

                    # First try a positive, quality-weighted average.

                    alpha <- quality

                } else {

                    # Then explore random linear combinations.

                    alpha <- rnorm(K) * sqrt(pmin(quality, 0.01))
                }

                # Normalize the coefficients.

                if (!all(is.finite(alpha)) || sum(abs(alpha)) < 1e-12) {
                    next
                }

                alpha <- alpha / sum(abs(alpha))


                # F_bar = sum_c alpha_c F_c

                F_bar <- matrix(0, nrow = R_try, ncol = R_try)

                for (k in seq_len(K)) {
                    F_bar <- F_bar + alpha[k] * li_F[[k]]
                }


                # Compute the eigendecomposition.

                eig <- tryCatch(
                    eigen(F_bar),
                    error = function(e) NULL
                )

                if (is.null(eig) || !all(is.finite(eig$values))) {
                    next
                }


                # Reject combinations with materially complex
                # eigenvalues: the desired directions are real.

                vals <- eig$values

                spectral_scale <- max(abs(vals))

                if (!is.finite(spectral_scale) || spectral_scale < 1e-12) {
                    next
                }

                imag_ratio <- max(abs(Im(vals))) / spectral_scale

                if (imag_ratio > tol_imag) {
                    next
                }


                # Minimum pairwise eigenvalue gap, normalized by
                # the spectral scale.

                D <- abs(outer(vals, vals, "-"))
                gap <- min(D[upper.tri(D)])

                score <- gap / spectral_scale

                if (is.finite(score) && score > best_score) {
                    best_score <- score
                    best_eig <- eig
                    best_imag <- imag_ratio
                    best_alpha <- alpha
                }
            }


            if (is.null(best_eig)) {
                stop(
                    paste(
                        "No suitable real eigendecomposition found for block", a,
                        "- try more configurations or inspect the eigenvalue gaps."
                    )
                )
            }


            #---------------------------------------------------
            # 2.3 Recover the directions for block a
            #---------------------------------------------------

            # The eigenvectors estimate the columns of Q_a,
            # up to permutation and scaling.

            Q_hat <- Re(best_eig$vectors)

            # Normalize the eigenvectors.

            Q_hat <- sweep(
                Q_hat, 2,
                sqrt(colSums(Q_hat^2)),
                "/"
            )


            # Return to the original block:
            # Lambda_a = U_a Q_a

            Lambda_star <- li_U[[a]] %*% Q_hat


            # Normalize each column to obtain directions.

            Lambda_star <- sweep(
                Lambda_star, 2,
                sqrt(colSums(Lambda_star^2)),
                "/"
            )


            # Fix the sign of each direction for reproducibility:
            # the coordinate with the largest absolute value is positive.

            for (r in 1:R_try) {

                #idx <- which.max(abs(Lambda_star[, r])) #mieux mais incoherent avec simu
                idx <- 1

                if (Lambda_star[idx, r] < 0) {
                    Lambda_star[, r] <- -Lambda_star[, r]
                }
            }

            Lambda_star
        })


        names(li_Lambda_star) <- names(A)

        # print("guess")
        # print(li_Lambda_star)
        li_Lambda_star_TRUE <- lapply(1:J, function(j){
             lambda_star_TRUE <- sweep(
                li_lambda_TRUE[[j]], 2,
                sqrt(colSums(li_lambda_TRUE[[j]]^2)),
                "/"
            )
            return(lambda_star_TRUE)
        })          
        # print("TRUE")
        # print(li_Lambda_star_TRUE)
    
    }

    li_vec_norm_and_Cjj <- lapply(
        1:J,
        function(j) {
            Lambda_star <- li_Lambda_star[[j]]
            pj <- pjs[j]

            # Within-block covariance matrix
            Sjj <- li_S[[j]][[j]]

            # Masking matrix Mj
            Mj <- matrix(1, nrow = pj, ncol = pj)
            diag(Mj) <- 0

            #------------------------------------------------
            # Equation (Eqx) of the paper
            #
            # Cjj_hat =
            #
            # [ (Lambda* %x% Lambda*)'
            #   diag(vec(Mj))
            #   (Lambda* %x% Lambda*) ]^(-1)
            #
            #   (Lambda* %x% Lambda*)'
            #   vec(Mj o Sjj)
            #------------------------------------------------

            Kjj <- kronecker(
                Lambda_star,
                Lambda_star
            )

            Gjj <- diag(as.vector(Mj)) %*% Kjj

            bjj <- as.vector(Mj * Sjj)

            #------------------------------------------------
            # Check whether Gjj has full rank
            #------------------------------------------------
            Gjj_cov <- t(Gjj) %*% Gjj
            rank_Gjj_cov <- qr(Gjj_cov)$rank
            full_rank <- rank_Gjj_cov == ncol(Gjj)


            if (!full_rank) {
                # Not full rank: non-unique least-squares solution

                print(paste(
                    "Block '", names(A)[j], "': Cjj estimation is not unique ",
                    "(rank = ", rank_Gjj_cov, ", full rank = ", ncol(Gjj), "). ",
                    "Regularizing the problem.",
                    sep = ""
                ))


                eigenvalues <- eigen(Gjj_cov, only.values = TRUE, symmetric = TRUE)$values
                max_eigen <- max(eigenvalues)
                eps <- 1e-8
                stop("error not unique")
            } else {
                eps <- 0
                max_eigen <- 1
            }


            A_cons <- matrix(0, nrow = ncol(Gjj), ncol = R_try)
            for (i in 1:R_try) {
                A_cons[(i - 1) * R_try + i, i] <- 1
            }
            b_0_cons <- rep(1e-10, R_try)
            # Il sera bon de checker la formule simple d'inversion si elle existe pour voir si on n'a bien aucun negatif sur la diagonale
            quadprog_fit <- quadprog::solve.QP(Gjj_cov + eps * max_eigen * diag(ncol(Gjj)), as.vector(t(Gjj) %*% bjj), A_cons, b_0_cons)
            Cjj_vec <- quadprog_fit$solution

            #------------------------------------------------
            # Cjj_hat is obtained by reshaping
            #------------------------------------------------

            Cjj_hat <- matrix(
                Cjj_vec,
                nrow = R_try,
                ncol = R_try,
                byrow = FALSE
            )

            if (any(diag(Cjj_hat) < 0)) {
                print(diag(Cjj_hat))
                stop("Negative diagonal elements in Cjj_hat. Please check the data and the model specification.")
            }
            vec_norm <- sqrt(diag(Cjj_hat))
            return(list(vec_norm = vec_norm, Cjj_hat = Cjj_hat))
        }
    )

    li_vec_norm <- lapply(li_vec_norm_and_Cjj, function(x) x$vec_norm)
    li_C_jj_hat <- lapply(li_vec_norm_and_Cjj, function(x) x$Cjj_hat)


    li_Lambda <- lapply(
        1:J,
        function(j) {
            Lambda_star <- li_Lambda_star[[j]]
            vec_norm <- li_vec_norm[[j]]
            Lambda <- sweep(Lambda_star, 2, vec_norm, FUN = "*")
            return(Lambda)
        }
    )


    # Pour le bon ordre : plus simple mais ne change rien after
    if(R_try > 1 & R_try == R_true){
        Lambda_1_TRUE <- li_lambda_TRUE[[1]]
        mat_ps <- t(t(li_Lambda_star[[1]]) %*% Lambda_1_TRUE)
        vec_reorder <- as.integer(clue::solve_LSAP(abs(mat_ps), maximum = TRUE))
        li_Lambda_star[[1]] <- li_Lambda_star[[1]][, vec_reorder]
        li_Lambda[[1]] <- li_Lambda[[1]][, vec_reorder]
        li_vec_norm[[1]] <- li_vec_norm[[1]][vec_reorder]
        li_C_jj_hat[[1]] <- li_C_jj_hat[[1]][vec_reorder, vec_reorder]
    }



    if(R_try > 1){
        if(do_full_reorder){
            reorder_full_fit <- reorder_full(li_S, li_Lambda_star, li_Lambda, li_vec_norm, li_C_jj_hat, R_try, J, ortho = ortho_solve_SVD)
            li_Lambda_star <- reorder_full_fit$li_Lambda_star
            li_Lambda <- reorder_full_fit$li_Lambda
            li_vec_norm <- reorder_full_fit$li_vec_norm
            li_C_jj_hat <- reorder_full_fit$li_C_jj_hat
        } else{
            reorder_vs_one_col_fit <- reorder_vs_one_col(li_S, li_Lambda_star, li_Lambda, li_vec_norm, li_C_jj_hat, R_try, J, ortho = ortho_solve_SVD)
            li_Lambda_star <- reorder_vs_one_col_fit$li_Lambda_star
            li_Lambda <- reorder_vs_one_col_fit$li_Lambda
            li_vec_norm <- reorder_vs_one_col_fit$li_vec_norm
            li_C_jj_hat <- reorder_vs_one_col_fit$li_C_jj_hat
        }
    }



    # print(li_Lambda)
    # # print("oui")
    # print(li_lambda_TRUE)
    # stop("ok")

    # print(lapply(1:J, function(j) round(li_Lambda[[j]] - li_lambda_TRUE[[j]], 4)))






    P_tilde <- matrix(0, nrow = R_try * J, ncol = R_try * J)
    for (i in 1:J) {
        for (j in i:J) {
            if (i == j) {
                D_j_inv <- diag(1 / li_vec_norm[[j]], nrow = R_try, ncol = R_try)
                P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- D_j_inv %*% li_C_jj_hat[[j]] %*% D_j_inv
            } else {
                D_i_inv <- diag(1 / li_vec_norm[[i]], nrow = R_try, ncol = R_try)
                D_j_inv <- diag(1 / li_vec_norm[[j]], nrow = R_try, ncol = R_try)
                S_ij <- li_S[[i]][[j]]
                if(ancien_P){
                    if(ortho_solve_SVD){
                        P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- D_i_inv %*% t(li_Lambda_star[[i]]) %*% S_ij %*% li_Lambda_star[[j]] %*% D_j_inv
                    } else{
                        p_inv_i <- solve(t(li_Lambda_star[[i]]) %*% li_Lambda_star[[i]]) %*% t(li_Lambda_star[[i]])
                        p_inv_j <- solve(t(li_Lambda_star[[j]]) %*% li_Lambda_star[[j]]) %*% t(li_Lambda_star[[j]])
                        P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- D_i_inv %*% p_inv_i %*% S_ij %*% t(p_inv_j) %*% D_j_inv
                    }
                } else{
                    P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- calc_P_ij_tilde(S_ij, li_Lambda_star, li_Lambda, li_vec_norm, i, j, R_try, ortho_solve_SVD, S_ii = li_S[[i]][[i]], S_jj = li_S[[j]][[j]])
                }
            }
        }
    }

    # P_tilde <- (P_tilde + t(P_tilde)) / 2 # Ensure symmetry

    # symmetrize
    for(i in 1:J){
        for(j in 1:J){
            if(i != j){
                if(j < i){
                    P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- t(P_tilde[((j - 1) * R_try + 1):(j * R_try), ((i - 1) * R_try + 1):(i * R_try)])
                }
            }
        }
    }



    Lambda_mat <- as.matrix(Matrix::bdiag(li_Lambda))

    ## Checker l'estimation at that stage


    Estim_Sigma <- Lambda_mat %*% P_tilde %*% t(Lambda_mat)

    diag(Estim_Sigma) <-  diag(S_full)




    # print("error sur P")
    # print(round(abs(as.matrix(P_tilde - P_TRUE)), 2))
    # print(round(as.matrix(P_TRUE), 2))
    # print(round(as.matrix(P_tilde), 2))

    
    # erreur_abs <- norm(Estim_Sigma - as.matrix(SIGMA_TRUE), "F")
    # error_SVD <- erreur_abs 
    # print(paste("Absolute error measurement model:", round(error_SVD, 4)))
    # print(li_Lambda)
    # print(li_lambda_TRUE)


    ## fin check estim

    var_MVs <- lapply(1:J, function(j) diag(li_S[[j]][[j]]))
    lv <- lvm(P_tilde, li_C)
    P_induced <- lv$P_induced
    # print(P_induced)
    SIGMA_IMPLIED_no_diag <- Lambda_mat %*% P_induced %*% t(Lambda_mat)
    SIGMA_IMPLIED <- SIGMA_IMPLIED_no_diag
    diag(SIGMA_IMPLIED) <- Reduce("cbind", var_MVs)
    Theta <- diag(SIGMA_IMPLIED) - diag(SIGMA_IMPLIED_no_diag)

    li_theta <- lapply(1:J, function(j) {
        theta_j <- Theta[(cumsum_pjs[j] + 1):(cumsum_pjs[j + 1])]
        #print(theta_j)
        return(theta_j)
    })


    # standardized loadings
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    if(R_try >1){
    li_std_Lambda <- lapply(1:J, function(j) {
        std_Lambda_j <- sweep(li_Lambda[[j]], 1, sqrt(var_MVs[[j]]), FUN = "/")
        return(std_Lambda_j)
    })
    } else{
        li_std_Lambda <- lapply(1:J, function(j) {
            std_Lambda_j <- li_Lambda[[j]] / sqrt(var_MVs[[j]])
            return(std_Lambda_j)
        })
    }


    # Measure of goodness-of(fit)

    T_LS <- d_LS(S_full, SIGMA_IMPLIED) # ecart quadratique entre la matrice de covariance empirique et la matrice de covariance impliquee par le modele
    # print(Lambda_mat)

    # print(P_tilde)
    # print(P_induced)


    if(!no_svd_refining){
        Sigma_implied <- SIGMA_IMPLIED
    } else{
        Sigma_implied <- Estim_Sigma
    }


    out <- list(
        li_Lambda_star = li_Lambda_star,
        li_Lambda = li_Lambda,
        li_std_Lambda = li_std_Lambda,
        Theta = Theta,
        li_gr = lv$li_gr,
        li_beta = lv$li_BETA,
        li_gamma = lv$li_GAMMA,
        li_R2 = lv$li_R2,
        li_PSI = lv$li_PSI,
        li_vec_norm = li_vec_norm,
        P_tilde = P_tilde,
        P_IMPLIED = P_induced,
        SIGMA_IMPLIED = Sigma_implied, 
        T_LS = T_LS,
        li_C = li_C,
        mode = mode,
        bias = bias,
        scale = scale,
        blocks = blocks,
        li_theta = li_theta,
        li_phi_endo = lv$li_phi_endo,
        li_phi_exo = lv$li_phi_exo,
        li_S = li_S
    )

    return(out)
}
