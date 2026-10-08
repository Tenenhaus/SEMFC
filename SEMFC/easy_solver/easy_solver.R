dyn.load("easy_solver.so") # chemin complet si besoin

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
  use_bfgs = TRUE,
  bfgs_curv_tol = 1e-8,
  bfgs_scale = TRUE
) {
    x <- as.numeric(pars)
    n <- length(x)

    if (is.null(eqfun)) eqfun <- function(x, ...) numeric(0)
    if (is.null(eqB)) eqB <- numeric(0)
    if (is.null(ineqLB)) ineqLB <- rep(-Inf, n)
    eqB <- as.numeric(eqB)
    ineqLB <- as.numeric(ineqLB)
    m <- length(eqB)

    if (length(ineqLB) != n) stop("ineqLB must have the same length as pars")
    if (m > 0 && length(eqfun(x, ...)) != m) {
        stop("eqB and eqfun have incompatible lengths")
    }

    f <- function(x) fun(x, ...)
    h <- function(x) eqfun(x, ...)

    # ordre identique a l'enum P_* du C
    ctrl <- as.numeric(c(
        maxit, tol, kkt_tol, stepsize, stepsize_min, stepsize_max,
        max_step, backtrack, f_min, eq_tol, active_tol, max_eq_correction,
        trace, use_bfgs, bfgs_curv_tol, bfgs_scale, grad_eps
    ))

    res <- .Call("easy_solver_c", x, ineqLB, eqB, ctrl, f, h, environment())

    hm <- res$history
    res$history <- data.frame(
        iter = as.integer(hm[, 1]), f = hm[, 2], kkt = hm[, 3], eq = hm[, 4],
        grad_norm = hm[, 5], step = hm[, 6], accepted = hm[, 7] == 1
    )
    res
}
