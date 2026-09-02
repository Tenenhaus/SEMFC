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


    # li_S <- lapply(A, function(x){
    #         lapply(A, function(y) cov2(x, y, bias = bias)) 
    #      })

    # On triche: SIGMA connu
    li_S <- lapply(1:J, function(j) {
        lapply(1:J, function(k) {
            S_jk <- as.matrix(SIGMA_TRUE[(cumsum_pjs[j] + 1):(cumsum_pjs[j + 1]), (cumsum_pjs[k] + 1):(cumsum_pjs[k + 1])])
        })
    })

    S_full <- matrix(0, nrow = sum(pjs), ncol = sum(pjs))
    for (i in 1:J) {
        for (j in 1:J) {
            S_full[(cumsum_pjs[i] + 1):(cumsum_pjs[i + 1]), (cumsum_pjs[j] + 1):(cumsum_pjs[j + 1])] <- li_S[[i]][[j]]
        }
    }

    #-------------------------------------------------------
    blocks <- A
    if (scale) {
        A <- lapply(A, function(x) scale2(x, bias = bias))
    } else {
        A <- lapply(A, function(x) scale2(x, bias = bias, scale = FALSE))
    }




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
            Lambda_star <- svd(mat_sum, nu = R_try, nv = 1)$u
            return(Lambda_star)

            # svd(t(A[[x]]) %*% Reduce("cbind", A[-x]),
            #     nu = 1, nv = 1
            # )$u
        }
    )

    # check for sign inversion
    li_Lambda_star <- lapply(li_Lambda_star, function(Lambda_star) {
        for (r in 1:R_try) {
            if (Lambda_star[1, r] < 0) {
                Lambda_star[, r] <- -Lambda_star[, r]
            }
        }
        return(Lambda_star)
    })
    names(li_Lambda_star) <- names(A)


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

    li_Lambda <- lapply(
        1:J,
        function(j) {
            Lambda_star <- li_Lambda_star[[j]]
            vec_norm <- li_vec_norm_and_Cjj[[j]]$vec_norm
            Lambda <- sweep(Lambda_star, 2, vec_norm, FUN = "*")
            return(Lambda)
        }
    )

    P_tilde <- matrix(0, nrow = R_try * J, ncol = R_try * J)
    for (i in 1:J) {
        for (j in 1:J) {
            if (i == j) {
                D_j_inv <- diag(1 / li_vec_norm_and_Cjj[[j]]$vec_norm, nrow = R_try, ncol = R_try)
                P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- D_j_inv %*% li_vec_norm_and_Cjj[[j]]$Cjj_hat %*% D_j_inv
            } else {
                D_i_inv <- diag(1 / li_vec_norm_and_Cjj[[i]]$vec_norm, nrow = R_try, ncol = R_try)
                D_j_inv <- diag(1 / li_vec_norm_and_Cjj[[j]]$vec_norm, nrow = R_try, ncol = R_try)
                S_ij <- li_S[[i]][[j]]
                P_tilde[((i - 1) * R_try + 1):(i * R_try), ((j - 1) * R_try + 1):(j * R_try)] <- D_i_inv %*% t(li_Lambda_star[[i]]) %*% S_ij %*% li_Lambda_star[[j]] %*% D_j_inv
                # print(t(li_Lambda_star[[i]]) %*% S_ij %*% li_Lambda_star[[j]])
            }
        }
    }

    P_tilde <- (P_tilde + t(P_tilde)) / 2 # Ensure symmetry
    Lambda_mat <- as.matrix(Matrix::bdiag(li_Lambda))

    ## Checker l'estimation at that stage


    Estim_Sigma <- Lambda_mat %*% P_tilde %*% t(Lambda_mat)

    diag(Estim_Sigma) <-  diag(S_full)



    erreur_abs <- d_LS(Estim_Sigma, SIGMA_TRUE)
    ratio_error_SVD <- erreur_abs / sum(diag(as.matrix(SIGMA_TRUE^2)))
    print(paste("Ratio of absolute error measurement model:", round(ratio_error_SVD, 4)))

    ## fin check estim

    var_MVs <- lapply(1:J, function(j) diag(li_S[[j]][[j]]))
    lv <- lvm(P_tilde, li_C)
    P_induced <- lv$P_induced
    SIGMA_IMPLIED_no_diag <- Lambda_mat %*% P_induced %*% t(Lambda_mat)
    SIGMA_IMPLIED <- SIGMA_IMPLIED_no_diag
    diag(SIGMA_IMPLIED) <- Reduce("cbind", var_MVs)
    Theta <- diag(SIGMA_IMPLIED) - diag(SIGMA_IMPLIED_no_diag)


    # standardized loadings
    # Expected value : -1 <= cor(y_jh, eta_j) <=1
    li_std_Lambda <- lapply(1:J, function(j) {
        std_Lambda_j <- sweep(li_Lambda[[j]], 1, sqrt(var_MVs[[j]]), FUN = "/")
        return(std_Lambda_j)
    })

    # Measure of goodness-of(fit)

    T_LS <- d_LS(S_full, SIGMA_IMPLIED) # ecart quadratique entre la matrice de covariance empirique et la matrice de covariance impliquee par le modele

    # print(P_tilde)
    # print(P_induced)


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
        SIGMA_IMPLIED = SIGMA_IMPLIED,
        T_LS = T_LS,
        li_C = li_C,
        mode = mode,
        bias = bias,
        scale = scale,
        blocks = blocks
    )

    return(out)
}
