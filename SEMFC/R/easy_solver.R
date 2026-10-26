easy_solver <- function(
  pars,
  fun,
  eqfun = NULL,
  eqB = NULL,
  ineqLB = NULL,
  ...,
  maxit = 1000,
  tol = 1e-8,
  kkt_tol = 1e-5,
  stepsize = 1e-3,
  stepsize_min = 1e-9,
  stepsize_max = 1,
  max_step = 1,
  backtrack = 0.5,
  f_min = -1e-10,
  eq_tol = 1e-6,
  grad_eps = 1e-6,
  active_tol = 1e-6,
  max_eq_correction = 20,
  trace = 1,
  use_bfgs = TRUE, # TRUE : direction quasi-Newton (BFGS inverse)
  bfgs_curv_tol = 1e-8, # seuil de la condition de courbure s'y > tol*|s||y|
  bfgs_scale = TRUE,
  n_cores = 10, # > 1 : evaluations des differences finies en parallele (Linux/macOS)
  cl = NULL
) {
    stop("fail: implementation in R instead of C")
    x <- as.numeric(pars)
    n <- length(x)

    if (is.null(eqfun)) {
        eqfun <- function(x, ...) numeric(0)
    }

    if (is.null(eqB)) {
        eqB <- numeric(0)
    }

    if (is.null(ineqLB)) {
        ineqLB <- rep(-Inf, n)
    }

    if (length(ineqLB) != n) {
        stop("ineqLB must have the same length as pars")
    }

    if (length(eqB) > 0) {
        if (length(eqfun(x, ...)) != length(eqB)) {
            stop("eqB and eqfun have incompatible lengths")
        }
    }


    # ------------------------------------------------------------
    # Objective
    # ------------------------------------------------------------

    f <- function(x) {
        fun(x, ...)
    }


    # ------------------------------------------------------------
    # Equality constraints
    # ------------------------------------------------------------

    h <- function(x) {
        eqfun(x, ...)
    }


    # ------------------------------------------------------------
    # Numerical gradient of the objective
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Differences finies centrees, paralleles, avec cache
    # ------------------------------------------------------------

    par_lapply <- function(X, FUN, ...) {
        res <- if (!is.null(cl)) {
            parallel::parLapply(cl, X, FUN, ...)
        } else if (n_cores > 1 && .Platform$OS.type != "windows") {
            parallel::mclapply(X, FUN, ..., mc.cores = n_cores)
        } else {
            lapply(X, FUN, ...)
        }

        if (any(vapply(res, inherits, logical(1), "try-error"))) {
            stop("Error during parallel finite differences")
        }

        res
    }

    fd_task <- function(j, x, do_f, do_h) {
        dx <- grad_eps * max(1, abs(x[j]))

        xp <- x
        xm <- x
        xp[j] <- xp[j] + dx
        xm[j] <- xm[j] - dx

        out <- list(dx = dx)

        if (do_f) {
            out$fp <- f(xp)
            out$fm <- f(xm)
        }

        if (do_h) {
            out$hp <- h(xp)
            out$hm <- h(xm)
        }

        out
    }

    compute_derivs <- function(x, do_f, do_h) {
        res <- par_lapply(seq_len(n), fd_task, x = x, do_f = do_f, do_h = do_h)

        dx <- vapply(res, function(r) r$dx, numeric(1))
        g <- NULL
        J <- NULL

        if (do_f) {
            fp <- vapply(res, function(r) r$fp, numeric(1))
            fm <- vapply(res, function(r) r$fm, numeric(1))

            if (any(!is.finite(c(fp, fm)))) {
                stop("Non finite objective during gradient calculation")
            }

            g <- (fp - fm) / (2 * dx)
        }

        if (do_h) {
            J <- do.call(
                cbind,
                lapply(res, function(r) (r$hp - r$hm) / (2 * r$dx))
            )
            J <- matrix(J, ncol = n)
        }

        list(g = g, J = J)
    }

    cache <- new.env()
    cache$x <- NULL
    cache$g <- NULL
    cache$J <- NULL

    reset_cache <- function(x) {
        cache$x <- x
        cache$g <- NULL
        cache$J <- NULL
    }

    same_x <- function(x) !is.null(cache$x) && identical(cache$x, x)

    # calcule gradient ET jacobienne en un seul passage parallel
    numerical_gradient <- function(x) {
        if (!same_x(x)) reset_cache(x)

        if (is.null(cache$g)) {
            d <- compute_derivs(x, TRUE, length(eqB) > 0)
            cache$g <- d$g
            if (!is.null(d$J)) cache$J <- d$J
        }

        cache$g
    }

    # jacobienne mise en cache
    numerical_jacobian <- function(x) {
        if (length(eqB) == 0) {
            return(matrix(0, 0, n))
        }

        if (!same_x(x)) reset_cache(x)

        if (is.null(cache$J)) {
            cache$J <- compute_derivs(x, FALSE, TRUE)$J
        }

        cache$J
    }

    # jacobienne exacte SANS cache (utilisee par equality_correction)
    jacobian_only <- function(x) {
        compute_derivs(x, FALSE, TRUE)$J
    }

    # ------------------------------------------------------------
    # Euclidean projection on lower bounds
    # ------------------------------------------------------------

    project_lb <- function(x) {
        pmax(x, ineqLB)
    }


    # ------------------------------------------------------------
    # Solve least squares problem
    #
    # Returns x solving approximately A x = b
    #
    # SVD is used only for robustness.
    # No Hessian is involved.
    # ------------------------------------------------------------

    least_squares <- function(A, b) {
        if (length(b) == 0) {
            return(numeric(0))
        }

        if (nrow(A) == 0 || ncol(A) == 0) {
            return(numeric(ncol(A)))
        }

        s <- svd(A)

        if (length(s$d) == 0) {
            return(rep(0, ncol(A)))
        }

        dmax <- max(s$d)

        if (dmax == 0) {
            return(rep(0, ncol(A)))
        }

        keep <- s$d > dmax * 1e-10

        if (!any(keep)) {
            return(rep(0, ncol(A)))
        }

        V <- s$v[, keep, drop = FALSE]
        U <- s$u[, keep, drop = FALSE]
        Dinv <- diag(1 / s$d[keep], nrow = sum(keep))

        as.numeric(V %*% Dinv %*% t(U) %*% b)
    }


    # ------------------------------------------------------------
    # Projection of a vector onto the null space of A
    #
    # Computes
    #
    #     d = v - A' (A A')^-1 A v
    #
    # which satisfies approximately
    #
    #     A d = 0
    #
    # ------------------------------------------------------------

    null_projection <- function(v, A) {
        if (nrow(A) == 0) {
            return(v)
        }

        lambda <- least_squares(A %*% t(A), A %*% v)

        v - as.numeric(t(A) %*% lambda)
    }

    null_projection_H <- function(v, A, H) {
        w <- as.numeric(H %*% v)

        if (nrow(A) == 0) {
            return(w)
        }

        AH <- A %*% H
        lambda <- least_squares(AH %*% t(A), as.numeric(AH %*% v))

        w - as.numeric(H %*% t(A) %*% lambda)
    }


    # ------------------------------------------------------------
    # Compute a feasible tangent direction
    #
    # Equality constraints:
    #
    #     J d = 0
    #
    # Active lower bounds:
    #
    #     d_i >= 0
    #
    # We use a simple active set strategy.
    #
    # If a direction wants to leave an active lower bound,
    # that variable is temporarily fixed with d_i = 0.
    # ------------------------------------------------------------

    feasible_direction <- function(x, g, J = NULL) {
        if (is.null(J)) J <- numerical_jacobian(x)

        active <- which(
            is.finite(ineqLB) &
                (x <= ineqLB + active_tol)
        )

        fixed <- integer(0) # <-- etait: fixed <- active

        for (iter in seq_len(n + 1)) {
            A <- J

            if (length(fixed) > 0) {
                E <- matrix(0, length(fixed), n)

                for (k in seq_along(fixed)) {
                    E[k, fixed[k]] <- 1
                }

                A <- rbind(A, E)
            }

            d <- if (use_bfgs) {
                null_projection_H(-g, A, H)
            } else {
                null_projection(-g, A)
            }

            if (length(active) == 0) {
                break
            }

            bad <- active[d[active] < -1e-14] # <-- etait: < -active_tol

            new_bad <- setdiff(bad, fixed)

            if (length(new_bad) == 0) {
                break
            }

            fixed <- c(fixed, new_bad)
        }

        nd <- sqrt(sum(d^2))

        if (!is.finite(nd) || nd == 0) {
            return(list(
                direction = rep(0, n),
                active = active,
                fixed = fixed,
                jacobian = J
            ))
        }

        if (nd > max_step) {
            d <- d * max_step / nd
        }

        list(
            direction = d,
            active = active,
            fixed = fixed,
            jacobian = J
        )
    }


    # ------------------------------------------------------------
    # Compute KKT residual
    #
    # We solve for lambda using only free variables:
    #
    #     g_free + J_free' lambda = 0
    #
    # Then
    #
    #     r = g + J' lambda
    #
    # For free variables:
    #
    #     r_i = 0
    #
    # For active lower bounds:
    #
    #     r_i >= 0
    #
    # ------------------------------------------------------------

    kkt_measure <- function(x, g = NULL) {
        if (is.null(g)) {
            g <- numerical_gradient(x)
        }

        hx <- h(x)

        if (length(hx) == 0) {
            J <- matrix(0, 0, n)
            lambda <- numeric(0)
            r <- g
        } else {
            J <- numerical_jacobian(x)

            active <- which(
                is.finite(ineqLB) &
                    (x <= ineqLB + active_tol)
            )

            free <- setdiff(seq_len(n), active)

            if (length(free) == 0) {
                lambda <- numeric(length(hx))
                r <- g
            } else {
                Jfree <- J[, free, drop = FALSE]

                A <- t(Jfree)

                lambda <- least_squares(
                    A,
                    -g[free]
                )

                r <- g + as.numeric(t(J) %*% lambda)
            }
        }

        active <- which(
            is.finite(ineqLB) &
                (x <= ineqLB + active_tol)
        )

        free <- setdiff(seq_len(n), active)

        free_residual <- if (length(free) > 0) {
            max(abs(r[free]))
        } else {
            0
        }

        active_violation <- if (length(active) > 0) {
            max(pmax(-r[active], 0))
        } else {
            0
        }

        eq_residual <- if (length(hx) > 0) {
            max(abs(hx - eqB))
        } else {
            0
        }

        # print(paste("free_residual:", free_residual, "active_violation:", active_violation, "eq_residual:", eq_residual))


        kkt <- max(
            free_residual,
            active_violation,
            eq_residual
        )

        list(
            kkt = kkt,
            residual = r,
            lambda = lambda,
            active = active,
            free = free,
            eq_residual = eq_residual,
            jacobian = J
        )
    }


    # ------------------------------------------------------------
    # Equality correction
    #
    # Given x, try to reduce
    #
    #     h(x) - eqB
    #
    # using the minimum norm linearized correction
    #
    #     J dx = -(h(x) - eqB)
    #
    # ------------------------------------------------------------

    equality_correction <- function(x, J0 = NULL) {
        if (length(eqB) == 0) {
            return(list(x = x, success = TRUE))
        }

        xcur <- x
        hcur <- h(xcur)
        r <- hcur - eqB
        J <- J0
        J_fresh <- FALSE

        for (iter in seq_len(max_eq_correction)) {
            r_norm <- max(abs(r))

            if (r_norm <= eq_tol) {
                return(list(x = xcur, success = TRUE))
            }

            if (is.null(J)) {
                J <- jacobian_only(xcur)
                J_fresh <- TRUE
            }

            dx <- least_squares(J, -r)
            ndx <- sqrt(sum(dx^2))

            if (!is.finite(ndx)) {
                return(list(x = xcur, success = FALSE))
            }

            if (ndx > max_step) {
                dx <- dx * max_step / ndx
            }

            accepted <- FALSE
            alpha <- 1

            while (alpha >= 1e-6) {
                xt <- project_lb(xcur + alpha * dx)
                ht <- h(xt)

                if (max(abs(ht - eqB)) < r_norm) {
                    accepted <- TRUE
                    break
                }

                alpha <- alpha * backtrack
            }

            if (accepted) {
                # Mise a jour de Broyden : J <- J + (y - J s) s' / s's
                s <- xt - xcur
                y <- ht - hcur
                ss <- sum(s^2)

                if (ss > 0) {
                    J <- J + outer(y - as.numeric(J %*% s), s) / ss
                }

                xcur <- xt
                hcur <- ht
                r <- ht - eqB
                J_fresh <- FALSE
            } else if (J_fresh) {
                return(list(x = xcur, success = FALSE))
            } else {
                J <- NULL # chec avec J approche : on recalcule la J exacte
            }
        }

        list(x = xcur, success = max(abs(r)) <= eq_tol)
    }

    # ------------------------------------------------------------
    # Initial point
    # ------------------------------------------------------------

    x <- project_lb(x)

    f_current <- f(x)
    # print(paste("Initial objective:", f_current))

    if (!is.finite(f_current)) {
        stop("Initial objective is not finite")
    }

    if (f_current < f_min) {
        stop("Initial objective is below f_min")
    }


    # ------------------------------------------------------------
    # Try to make the initial point feasible
    # ------------------------------------------------------------

    if (length(eqB) > 0) {
        corr <- equality_correction(x)

        if (corr$success) {
            fcorr <- f(corr$x)

            if (
                is.finite(fcorr) &&
                    fcorr >= f_min
            ) {
                x <- corr$x
                f_current <- fcorr
            }
        }
    }


    # ------------------------------------------------------------
    # History
    # ------------------------------------------------------------

    history <- data.frame(
        iter = integer(0),
        f = numeric(0),
        kkt = numeric(0),
        eq = numeric(0),
        grad_norm = numeric(0),
        step = numeric(0),
        accepted = logical(0)
    )


    converged <- FALSE
    reason <- "maxit reached"


    # ------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------

    H <- diag(n)
    H_is_identity <- TRUE
    x_old <- NULL
    g_old <- NULL
    J_old <- NULL

    for (iter in seq_len(maxit)) {
        g <- numerical_gradient(x)

        km <- kkt_measure(x, g)

        kkt <- km$kkt
        eq_norm <- km$eq_residual
        grad_norm <- sqrt(sum(g^2))


        if (
            trace > 0 &&
                (iter == 1 || iter %% trace == 0)
        ) {
            cat(
                sprintf(
                    paste(
                        "iter %4d |",
                        "f = % .6e |",
                        "kkt = %.3e |",
                        "eq = %.3e |",
                        "grad = %.3e |",
                        "step = %.3e\n"
                    ),
                    iter,
                    f_current,
                    kkt,
                    eq_norm,
                    grad_norm,
                    stepsize
                )
            )
        }


        # --------------------------------------------------------
        # KKT stopping criterion
        # --------------------------------------------------------

        if (kkt <= kkt_tol) {
            converged <- TRUE
            reason <- "KKT tolerance reached"

            break
        }

        # --------------------------------------------------------
        # Mise a jour BFGS (inverse), sur le gradient du Lagrangien
        # --------------------------------------------------------
        if (use_bfgs && !is.null(x_old)) {
            s <- x - x_old
            y <- g - g_old

            if (length(km$lambda) > 0) {
                y <- y + as.numeric(t(km$jacobian - J_old) %*% km$lambda)
            }

            sy <- sum(s * y)

            # Condition de courbure : si elle echoue (cas non convexe),
            # on saute la mise a jour pour garder H definie positive
            if (is.finite(sy) &&
                sy > bfgs_curv_tol * sqrt(sum(s^2)) * sqrt(sum(y^2))) {
                if (bfgs_scale && H_is_identity) {
                    H <- diag(sy / sum(y^2), n)
                }

                # rho <- 1 / sy
                # I_n <- diag(n)

                # H <- (I_n - rho * s %*% t(y)) %*% H %*% (I_n - rho * y %*% t(s)) +
                #     rho * s %*% t(s)

                Hy <- as.numeric(H %*% y)
                H <- H +
                    ((sy + sum(y * Hy)) / sy^2) * outer(s, s) -
                    (outer(Hy, s) + outer(s, Hy)) / sy

                H_is_identity <- FALSE
            }
        }

        fd <- feasible_direction(x, g, km$jacobian) # <-- J reused


        # --------------------------------------------------------
        # Compute feasible descent direction
        # --------------------------------------------------------

        fd <- feasible_direction(x, g)

        d <- fd$direction

        nd <- sqrt(sum(d^2))


        # --------------------------------------------------------
        # No admissible descent direction
        #
        # This is also a possible KKT point.
        # --------------------------------------------------------

        if (!is.finite(nd) || nd <= tol) {
            km <- kkt_measure(x, g)

            if (km$kkt <= kkt_tol) {
                converged <- TRUE
                reason <- "No feasible descent direction"
            } else {
                reason <- "No usable feasible descent direction"
            }

            break
        }


        # --------------------------------------------------------
        # Backtracking
        # --------------------------------------------------------

        alpha <- if (use_bfgs) stepsize_max else min(stepsize, stepsize_max)

        accepted <- FALSE

        while (alpha >= stepsize_min) {
            x_trial <- x + alpha * d

            x_trial <- project_lb(x_trial)

            f_trial <- f(x_trial)


            # ----------------------------------------------------
            # Objective barrier
            # ----------------------------------------------------

            if (
                !is.finite(f_trial) ||
                    f_trial < f_min
            ) {
                alpha <- alpha * backtrack
                next
            }


            # ----------------------------------------------------
            # Monotonicity
            # ----------------------------------------------------

            if (f_trial >= f_current) {
                alpha <- alpha * backtrack
                next
            }


            # ----------------------------------------------------
            # Equality correction
            # ----------------------------------------------------

            if (length(eqB) > 0) {
                corr <- equality_correction(x_trial)

                if (!corr$success) {
                    alpha <- alpha * backtrack
                    next
                }

                x_trial2 <- corr$x

                f_trial2 <- f(x_trial2)

                if (
                    !is.finite(f_trial2) ||
                        f_trial2 < f_min
                ) {
                    alpha <- alpha * backtrack
                    next
                }

                if (f_trial2 > f_current) {
                    alpha <- alpha * backtrack
                    next
                }

                x_trial <- x_trial2
                f_trial <- f_trial2
            }


            # ----------------------------------------------------
            # Accept
            # ----------------------------------------------------

            accepted <- TRUE
            break
        }


        # --------------------------------------------------------
        # No accepted step
        # --------------------------------------------------------


        if (!accepted) {
            history <- rbind(
                history,
                data.frame(
                    iter = iter, f = f_current, kkt = kkt, eq = eq_norm,
                    grad_norm = grad_norm, step = 0, accepted = FALSE
                )
            )

            if (use_bfgs) {
                if (H_is_identity) {
                    reason <- "Minimum stepsize reached"
                    break
                }
                H <- diag(n)
                H_is_identity <- TRUE
                x_old <- NULL
                next
            }

            stepsize <- stepsize * backtrack

            if (stepsize < stepsize_min) {
                reason <- "Minimum stepsize reached"
                break
            }

            next
        }


        # --------------------------------------------------------
        # Accept the step
        # --------------------------------------------------------

        if (use_bfgs) {
            x_old <- x
            g_old <- g
            J_old <- km$jacobian
        }

        x <- x_trial
        f_current <- f_trial


        history <- rbind(
            history,
            data.frame(
                iter = iter,
                f = f_current,
                kkt = kkt,
                eq = eq_norm,
                grad_norm = grad_norm,
                step = alpha,
                accepted = TRUE
            )
        )


        # --------------------------------------------------------
        # Very conservative stepsize increase
        # --------------------------------------------------------

        if (alpha >= 0.9 * stepsize) {
            stepsize <- min(
                stepsize * 1.1,
                stepsize_max
            )
        }

        if (alpha < 0.9 * stepsize) {
            stepsize <- max(
                stepsize * 0.7,
                stepsize_min
            )
        }

        # print(paste("chosen alpha:", alpha))
        # print(paste("new stepsize:", stepsize))
    }


    # ------------------------------------------------------------
    # Final diagnostics
    # ------------------------------------------------------------

    g_final <- numerical_gradient(x)
    km_final <- kkt_measure(x, g_final)

    list(
        pars = x,
        value = f_current,
        convergence = if (converged) 0 else 1,
        message = reason,
        gradient = g_final,
        kkt = km_final$kkt,
        kkt_residual = km_final$residual,
        lambda = km_final$lambda,
        eq = h(x),
        eq_residual = km_final$eq_residual,
        active = km_final$active,
        iterations = iter,
        stepsize = stepsize,
        history = history
    )
}
