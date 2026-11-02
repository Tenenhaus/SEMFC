# =============================================================================
# Co-diagonalisation non orthogonale de matrices NON symetriques
#
# Probleme : on a K matrices R x R  F_1, ..., F_K  telles que
#     F_c = A  Diag(d_c)  A^{-1}        (A inversible, non orthogonale)
# avec du bruit. On cherche A (colonnes = vecteurs propres communs).
#
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

# --- critere : somme ponderee des normes^2 des parties hors diagonale ---------
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
    set.seed(1) # reproductible ; la combinaison aleatoire evite les valeurs propres confondues
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

# =============================================================================
# TOY EXAMPLE  (executer le fichier avec Rscript, ou sourcer et lancer toy())
# =============================================================================
toy <- function(R = 4, K = 20, noise = 0.02, nrep = 200, seed = 123) {
    set.seed(seed)
    # distance entre colonnes estimees et vraies, a permutation et signe pres
    col_err <- function(Ahat, A) {
        Ahat <- sweep(Ahat, 2, sqrt(colSums(Ahat^2)), "/")
        A <- sweep(A, 2, sqrt(colSums(A^2)), "/")
        S <- abs(crossprod(A, Ahat)) # |cos| (vrai x estime)
        err <- numeric(R)
        for (k in 1:R) { # appariement glouton
            ij <- which(S == max(S), arr.ind = TRUE)[1, ]
            err[k] <- 1 - S[ij[1], ij[2]]
            S[ij[1], ] <- -1
            S[, ij[2]] <- -1
        }
        mean(err)
    }
    res <- matrix(NA, nrep, 2, dimnames = list(NULL, c("une_combinaison", "joint_diag")))
    for (rep in 1:nrep) {
        # A volontairement non orthogonale, assez mal conditionnee
        A <- matrix(rnorm(R * R), R, R) + 0.8 * outer(rnorm(R), rnorm(R)) + diag(R)
        Ainv <- solve(A)
        Fs <- lapply(1:K, function(c) {
            d <- rnorm(R)
            F <- A %*% diag(d) %*% Ainv
            F <- F / sqrt(sum(F^2)) # normalisation de Frobenius
            F + noise * matrix(rnorm(R * R), R, R) / R # bruit
        })
        # (1) methode actuelle : UNE seule paire (u, v) => F = M(u) M(v)^{-1}
        #     ici : F_1 seule, diagonalisee par eigen (equivalent a K = 1)
        e <- eigen(Fs[[1]])
        A1 <- Re(e$vectors)
        # (2) co-diagonalisation des K matrices
        fit <- joint_diag_nonorth(Fs)
        res[rep, ] <- c(col_err(A1, A), col_err(fit$A, A))
    }
    cat(sprintf("R = %d, K = %d, bruit = %.3f, %d repetitions\n", R, K, noise, nrep))
    cat("Erreur moyenne (1 - |cos|) sur les colonnes, mediane / moyenne :\n")
    print(rbind(mediane = apply(res, 2, median), moyenne = colMeans(res)))
    cat("\nSans bruit (doit etre ~1e-15) :\n")
    A <- matrix(rnorm(R * R), R, R) + diag(R)
    Ainv <- solve(A)
    Fs <- lapply(1:K, function(c) A %*% diag(rnorm(R)) %*% Ainv)
    fit <- joint_diag_nonorth(Fs)
    cat(
        "  erreur colonnes :", col_err(fit$A, A), " | critere final :", fit$crit,
        " | iterations :", fit$iter, "\n"
    )
    invisible(res)
}

# toy()
