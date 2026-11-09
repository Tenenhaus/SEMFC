# =============================================================================
# Co-diagonalisation non orthogonale de matrices NON symetriques

# Probleme : on a K matrices R x R  F_1, ..., F_K  telles que
#     F_c = A  Diag(d_c)  A^{-1}        (A inversible, non orthogonale)
# avec du bruit. On cherche A (colonnes = vecteurs propres communs).

# Methode : iterations de type "FFDIAG" adaptees au cas non symetrique.
#   On maintient W ~ A^{-1} et C_c = W F_c W^{-1}.
#   A chaque iteration on cherche E (diag nulle) tel que
#       (I+E) C_c (I+E)^{-1}  ~  C_c + E C_c - C_c E
#   ait des termes hors diagonale nuls. Au premier ordre, le terme (i,j) de
#   E C - C E vaut E_ij (d_j - d_i), et les equations (i,j) et (j,i) sont
#   decouplees (cas non symetrique), d'ou une formule explicite :
#       E_ij = - sum_c w_c C_c[i,j] (d_c[j]-d_c[i]) / sum_c w_c (d_c[j]-d_c[i])^2
#   Cout par iteration : O(K R^3). Convergence en quelques dizaines d'iterations.
# =============================================================================

ancien <- FALSE
# seed_used <- 2

# --- critere : somme ponderee des normes^2 des parties hors diagonale ---------
if (ancien) {
    .off_crit <- function(Cs, w) {
        sum(vapply(seq_along(Cs), function(c) {
            C <- Cs[[c]]
            diag(C) <- 0
            w[c] * sum(C^2)
        }, numeric(1)))
    }

    # --- initialisation : vecteurs propres (reels) d'une combinaison de F_c -------
    .init_A <- function(Fs, w) {
        K <- length(Fs)
        R <- nrow(Fs[[1]])
        alpha <- w * rnorm(K)
        F <- Reduce(`+`, Map(function(M, a) a * M, Fs, alpha))
        e <- eigen(F)
        V <- e$vectors
        A <- matrix(0, R, R)
        r <- 1
        while (r <= R) {
            if (abs(Im(e$values[r])) < 1e-10 * max(1, abs(e$values[r])) || r == R) {
                A[, r] <- Re(V[, r])
                r <- r + 1
            } else { # paire conjuguee : Re et Im donnent une base reelle
                A[, r] <- Re(V[, r])
                A[, r + 1] <- Im(V[, r])
                r <- r + 2
            }
        }
        A
    }

    # --- fonction principale -------------------------------------------------------
    # Fs      : liste de K matrices R x R
    # weights : poids des matrices (defaut 1)
    # A0      : initialisation (colonnes = directions). Si NULL, vecteurs propres d'une
    #           combinaison aleatoire.
    # Retourne : A (colonnes de norme 1), D (K x R, valeurs propres jointes),
    #            crit (critere final), iter, crit_path
    joint_diag_nonorth <- function(Fs, weights = NULL, A0 = NULL,
                                   maxit = 200, tol = 1e-12, ridge = 1e-10,
                                   verbose = FALSE) {
        K <- length(Fs)
        R <- nrow(Fs[[1]])
        w <- if (is.null(weights)) rep(1, K) else weights
        A <- if (is.null(A0)) .init_A(Fs, w) else A0
        A <- sweep(A, 2, sqrt(colSums(A^2)), "/")
        W <- solve(A)
        conj <- function(W) lapply(Fs, function(F) W %*% F %*% solve(W))
        Cs <- conj(W)
        crit <- .off_crit(Cs, w)
        path <- crit

        for (it in seq_len(maxit)) {
            # diagonales des C_c : matrice K x R
            Dm <- t(vapply(Cs, diag, numeric(R)))
            E <- matrix(0, R, R)
            for (i in 1:R) {
                for (j in 1:R) {
                    if (i != j) {
                        dd <- Dm[, j] - Dm[, i]
                        num <- sum(w * vapply(Cs, function(C) C[i, j], numeric(1)) * dd)
                        den <- sum(w * dd^2)
                        E[i, j] <- -num / (den + ridge * (1 + den))
                    }
                }
            }
            # recherche de pas (backtracking) pour garantir la decroissance du critere
            step <- 1
            ok <- FALSE
            while (step > 1e-4) {
                Wn <- (diag(R) + step * E) %*% W
                An <- tryCatch(solve(Wn), error = function(e) NULL)
                if (!is.null(An)) {
                    An <- sweep(An, 2, sqrt(colSums(An^2)), "/") # colonnes de norme 1
                    Wn <- solve(An)
                    Cn <- conj(Wn)
                    cn <- .off_crit(Cn, w)
                    if (is.finite(cn) && cn <= crit) {
                        ok <- TRUE
                        break
                    }
                }
                step <- step / 2
            }
            if (!ok) break
            rel <- (crit - cn) / max(crit, .Machine$double.eps)
            W <- Wn
            Cs <- Cn
            crit <- cn
            path <- c(path, crit)
            if (verbose) cat(sprintf("it %3d  crit %.3e  step %.3f\n", it, crit, step))
            if (rel < tol || crit < 1e-28) break
        }
        A <- solve(W)
        A <- sweep(A, 2, sqrt(colSums(A^2)), "/")
        list(
            A = A, D = t(vapply(conj(solve(A)), diag, numeric(R))),
            crit = crit, iter = length(path) - 1, crit_path = path
        )
    }
} else {
    joint_diag_nonorth <- function(Fs, weights = NULL, A0 = NULL,
                                   maxit = 100, tol = 1e-12,
                                   ridge = 1e-10, alpha = 0.25,
                                   nstart = 5, real_tol = 1e-10,
                                   verbose = FALSE) {
        # ----------------------------------------------------------
        # 0. Input checks
        # ----------------------------------------------------------

        K <- length(Fs)

        if (K < 2) {
            stop("At least two matrices are required.")
        }

        if (!all(vapply(Fs, is.matrix, logical(1)))) {
            stop("All elements of Fs must be matrices.")
        }

        if (any(vapply(Fs, is.complex, logical(1)))) {
            stop("Complex input matrices are not supported.")
        }

        R <- nrow(Fs[[1]])

        if (R != ncol(Fs[[1]]) || R < 1) {
            stop("Matrices must be square.")
        }

        for (k in seq_len(K)) {
            if (nrow(Fs[[k]]) != R ||
                ncol(Fs[[k]]) != R) {
                stop("All matrices must have the same dimensions.")
            }

            if (any(!is.finite(Fs[[k]]))) {
                stop("Matrices must contain finite values.")
            }

            storage.mode(Fs[[k]]) <- "double"
        }

        if (maxit < 1 || nstart < 1) {
            stop("maxit and nstart must be positive.")
        }

        if (alpha <= 0 || alpha >= 1) {
            stop("alpha must belong to (0, 1).")
        }

        if (ridge < 0 || tol <= 0 || real_tol <= 0) {
            stop("ridge must be non-negative; tolerances must be positive.")
        }

        # weights is intentionally not used.

        # ----------------------------------------------------------
        # 1. Utility functions
        # ----------------------------------------------------------

        normalize_columns <- function(A) {
            norms <- sqrt(colSums(A^2))

            if (any(!is.finite(norms)) ||
                any(norms < .Machine$double.eps)) {
                stop("A contains invalid or zero columns.")
            }

            sweep(A, 2, norms, "/")
        }

        off_matrix <- function(M) {
            diag(M) <- 0
            M
        }

        criterion <- function(Gs) {
            sum(vapply(Gs, function(G) {
                sum(off_matrix(G)^2)
            }, numeric(1)))
        }

        safe_inverse <- function(B) {
            if (any(!is.finite(B))) {
                return(NULL)
            }

            rc <- tryCatch(
                rcond(B),
                error = function(e) NA_real_
            )

            if (!is.finite(rc) || rc < .Machine$double.eps) {
                return(NULL)
            }

            tryCatch(
                solve(B),
                error = function(e) NULL
            )
        }

        transform_matrices <- function(B) {
            Binv <- safe_inverse(B)

            if (is.null(Binv)) {
                return(NULL)
            }

            Gs <- lapply(Fs, function(F) {
                B %*% F %*% Binv
            })

            if (any(!vapply(Gs, function(G) {
                all(is.finite(G))
            }, logical(1)))) {
                return(NULL)
            }

            Gs
        }

        # ----------------------------------------------------------
        # 2. Random real initialization
        # ----------------------------------------------------------

        random_initialization <- function() {
            for (attempt in seq_len(100)) {
                coeff <- rnorm(K)

                M <- matrix(0.0, R, R)

                for (k in seq_len(K)) {
                    M <- M + coeff[k] * Fs[[k]]
                }

                eig <- tryCatch(
                    eigen(M),
                    error = function(e) NULL
                )

                if (is.null(eig)) {
                    next
                }

                # Check eigenvalues.

                eigval_scale <- max(1, max(abs(Re(eig$values))))

                eigval_imag <- max(abs(Im(eig$values)))

                if (eigval_imag > real_tol * eigval_scale) {
                    next
                }

                # Check eigenvectors.

                eigvec_scale <- max(1, max(abs(Re(eig$vectors))))

                eigvec_imag <- max(abs(Im(eig$vectors)))

                if (eigvec_imag > real_tol * eigvec_scale) {
                    next
                }

                # Keep only numerically real eigenvectors.

                A <- Re(eig$vectors)

                if (any(!is.finite(A))) {
                    next
                }

                A <- tryCatch(
                    normalize_columns(A),
                    error = function(e) NULL
                )

                if (is.null(A)) {
                    next
                }

                B <- safe_inverse(A)

                if (!is.null(B)) {
                    return(A)
                }
            }

            stop(
                paste0(
                    "Unable to find a well-conditioned real eigenbasis. ",
                    "Try supplying a real A0 or increasing the number ",
                    "of random initializations."
                )
            )
        }

        # ----------------------------------------------------------
        # 3. One JDTE-EP optimization run
        # ----------------------------------------------------------

        run_one <- function(A_init) {
            if (is.complex(A_init)) {
                stop("A0 must be real.")
            }

            A_init <- normalize_columns(A_init)

            B <- safe_inverse(A_init)

            if (is.null(B)) {
                stop("Initial matrix is singular or ill-conditioned.")
            }

            Gs <- transform_matrices(B)

            if (is.null(Gs)) {
                stop("Initial matrix transformation failed.")
            }

            crit <- criterion(Gs)

            if (!is.finite(crit)) {
                stop("Initial criterion is not finite.")
            }

            crit_path <- crit
            iter <- 0

            # ------------------------------------------------------
            # Main loop
            # ------------------------------------------------------

            for (it in seq_len(maxit)) {
                # Exterior penalty factor.

                mu <- alpha^(it - 1)

                Lambda <- lapply(Gs, diag)
                Off <- lapply(Gs, off_matrix)

                # First-order correction matrix.

                Z <- matrix(0.0, R, R)

                for (i in seq_len(R)) {
                    for (j in seq_len(R)) {
                        if (i == j) {
                            next
                        }

                        delta <- vapply(seq_len(K), function(k) {
                            Lambda[[k]][i] - Lambda[[k]][j]
                        }, numeric(1))

                        off_ij <- vapply(seq_len(K), function(k) {
                            Off[[k]][i, j]
                        }, numeric(1))

                        numerator <- sum(delta * off_ij)

                        denominator <- sum(delta^2) + mu + ridge

                        Z[i, j] <- numerator / denominator
                    }
                }

                # --------------------------------------------------
                # Line search
                # --------------------------------------------------

                step <- 1
                accepted <- FALSE

                for (ls in seq_len(15)) {
                    X <- diag(1.0, R) + step * Z

                    B_new <- X %*% B

                    Gs_new <- transform_matrices(B_new)

                    if (!is.null(Gs_new)) {
                        crit_new <- criterion(Gs_new)

                        if (is.finite(crit_new) && crit_new < crit) {
                            accepted <- TRUE
                            break
                        }
                    }

                    step <- step / 2
                }

                if (!accepted) {
                    if (verbose) {
                        cat(
                            "No improving step at iteration",
                            it, "\n"
                        )
                    }

                    break
                }

                # Update.

                relative_change <- abs(crit - crit_new) /
                    max(crit, .Machine$double.eps)

                B <- B_new
                Gs <- Gs_new
                crit <- crit_new

                iter <- it
                crit_path <- c(crit_path, crit)

                if (verbose && (it == 1 || it %% 10 == 0)) {
                    cat(
                        "Iteration:", it,
                        "| Criterion:", signif(crit, 8),
                        "| Mu:", signif(mu, 5),
                        "| Step:", signif(step, 4),
                        "\n"
                    )
                }

                if (relative_change < tol) {
                    break
                }
            }

            # ------------------------------------------------------
            # Recover A and joint eigenvalues
            # ------------------------------------------------------

            Binv <- safe_inverse(B)

            if (is.null(Binv)) {
                stop("Final matrix B is singular.")
            }

            A <- normalize_columns(Binv)

            # Strict real-valued output.

            if (is.complex(A)) {
                stop("Internal error: complex matrix encountered.")
            }

            D <- matrix(0.0, K, R)

            for (k in seq_len(K)) {
                D[k, ] <- diag(Gs[[k]])
            }

            list(
                A = A,
                D = D,
                criterion = crit,
                iterations = iter,
                criterion_path = crit_path
            )
        }

        # ----------------------------------------------------------
        # 4. Multiple initializations
        # ----------------------------------------------------------

        best <- NULL

        total_starts <- nstart + if (!is.null(A0)) 1 else 0

        for (s in seq_len(total_starts)) {
            if (!is.null(A0) && s == 1) {
                if (!is.matrix(A0) ||
                    nrow(A0) != R ||
                    ncol(A0) != R) {
                    stop("A0 must be an R x R matrix.")
                }

                if (is.complex(A0)) {
                    stop("A0 must be real.")
                }

                if (any(!is.finite(A0))) {
                    stop("A0 must contain finite values.")
                }

                A_init <- A0
            } else {
                A_init <- random_initialization()
            }

            if (verbose) {
                cat("\nStarting run", s, "of", total_starts, "\n")
            }

            current <- tryCatch(
                run_one(A_init),
                error = function(e) {
                    message(
                        "Initialization ", s, " failed: ",
                        conditionMessage(e)
                    )
                    NULL
                }
            )

            if (is.null(current)) {
                next
            }

            if (is.null(best) ||
                current$criterion < best$criterion) {
                best <- current
            }
        }

        # ----------------------------------------------------------
        # 5. Final result
        # ----------------------------------------------------------

        if (is.null(best)) {
            stop("All optimization runs failed.")
        }

        if (verbose) {
            cat("\nOptimization finished.\n")
            cat("Final criterion:", best$criterion, "\n")
            cat("Iterations:", best$iterations, "\n")
        }

        best
    }
}
# # =============================================================================
# # TOY EXAMPLE  (Rscript joint_diag_nonorth.R   ou   source(...); toy())
# #  1) bruit homogene : une seule matrice (eigen) vs shear-rotation sur K matrices
# #  2) bruit heterogene : sans poids vs poids w_n = 1 / sigma_n^2
# #  3) sans bruit : depart par defaut, puis depart A0 = diag(R) (celui du papier)
# # =============================================================================
toy <- function(R = 4, K = 20, noise = 0.02, nrep = 200, seed = 123, no_weight = TRUE) {
    set.seed(seed)
    cn <- function(M) sweep(M, 2, sqrt(colSums(M^2)), "/")
    # erreur moyenne 1-|cos| sur les colonnes, a permutation et signe pres
    col_err <- function(Ahat, A) {
        S <- abs(crossprod(cn(A), cn(Ahat)))
        err <- numeric(R)
        for (k in 1:R) { # appariement glouton
            ij <- which(S == max(S), arr.ind = TRUE)[1, ]
            err[k] <- 1 - S[ij[1], ij[2]]
            S[ij[1], ] <- -1
            S[, ij[2]] <- -1
        }
        mean(err)
    }
    rand_A <- function() matrix(rnorm(R * R), R, R) + 0.8 * outer(rnorm(R), rnorm(R)) + diag(R)
    make_F <- function(A, sd_c) {
        Ai <- solve(A)
        lapply(seq_along(sd_c), function(c) {
            F <- A %*% diag(rnorm(R)) %*% Ai
            F <- F / sqrt(sum(F^2))
            F + sd_c[c] * matrix(rnorm(R * R), R, R) / R
        })
    }

    # ---- 1) bruit homogene
    res1 <- matrix(NA, nrep, 2, dimnames = list(NULL, c("une_matrice", "shear_rotation")))
    li_A <- lapply(1:nrep, function(rep) rand_A())
    li_Fs <- lapply(li_A, function(A) make_F(A, rep(noise, K)))
    for (rep in 1:nrep) {
        A <- li_A[[rep]]
        Fs <- li_Fs[[rep]]
        res1[rep, ] <- c(
            col_err(Re(eigen(Fs[[1]])$vectors), A),
            col_err(joint_diag_nonorth(Fs)$A, A)
        )
        # print(joint_diag_nonorth(Fs)$A)
    }
    cat(sprintf("1) Bruit homogene : R=%d, K=%d, sigma=%.3f, %d repetitions\n", R, K, noise, nrep))
    cat("   Erreur 1-|cos| :\n")
    print(rbind(mediane = apply(res1, 2, median), moyenne = colMeans(res1)))
    if (no_weight) {
        return("done")
    }

    # ---- 2) bruit heterogene : moitie des matrices avec sigma x 20
    sd_c <- rep(c(noise, 20 * noise), each = K / 2)
    res2 <- matrix(NA, nrep, 2, dimnames = list(NULL, c("sans_poids", "poids_1/sigma2")))
    for (rep in 1:nrep) {
        A <- rand_A()
        Fs <- make_F(A, sd_c)
        res2[rep, ] <- c(
            col_err(joint_diag_nonorth(Fs)$A, A),
            col_err(joint_diag_nonorth(Fs, weights = 1 / sd_c^2)$A, A)
        )
    }
    cat("\n2) Bruit heterogene (moitie des matrices avec sigma x20) :\n")
    print(rbind(mediane = apply(res2, 2, median), moyenne = colMeans(res2)))

    # ---- 3) sans bruit
    A <- rand_A()
    Fs <- make_F(A, rep(0, K))
    f1 <- joint_diag_nonorth(Fs)
    f2 <- joint_diag_nonorth(Fs, A0 = diag(R))
    cat("\n3) Sans bruit (attendu ~ 1e-15)\n")
    cat("   depart par defaut  : erreur =", col_err(f1$A, A), "| sweeps =", f1$iter, "\n")
    cat("   depart A0 = diag(R): erreur =", col_err(f2$A, A), "| sweeps =", f2$iter, "\n")
    cat("   critere par sweep (depart A0 = diag(R)) :\n")
    print(signif(f2$crit_path, 3))
    invisible(list(res1 = res1, res2 = res2))
}

# toy(seed = seed_used)
ancien <- NULL
