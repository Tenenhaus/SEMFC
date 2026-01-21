library(PMA)





new_SPC <- function(x, sumabsv = 4, niter = 20, K = 1, orth = FALSE, trace = TRUE, v = NULL, center = TRUE, cnames = NULL, vpos = FALSE, vneg = FALSE, compute.pve = TRUE) {
  if (vpos && vneg) stop("Cannot constrain elements to be positive AND negative.")
  out <- PMA:::PMDL1L1(x, sumabsu = sumabsv, sumabsv = sumabsv, niter = niter, K = K, orth = orth, trace = trace, v = v, center = center, cnames = cnames, upos = FALSE, uneg = FALSE, vpos = vpos, vneg = vneg)
  if (compute.pve) {
    v <- matrix(out$v, ncol = K)
    ve <- NULL
    xfill <- x
    if (center) xfill <- x - PMA:::mean_na(x)
    xfill[is.na(x)] <- PMA:::mean_na(xfill)
    for (k in 1 : K) {
      vk <- matrix(v[, 1 : k], ncol = k)
      xk <- xfill %*% vk %*% solve(t(vk) %*% vk) %*% t(vk)
      svdxk <- svd(xk)
      ve <- c(ve, sum(svdxk$d ^ 2))
    }
    pve <- ve / sum(svd(xfill)$d ^ 2)
    out$prop.var.explained <- pve
  }
  out$vpos <- vpos
  out$vneg <- vneg
  class(out) <- "SPC"
  return(out)
}


new_SPC.cv <- function(x, sumabsvs = seq(1.2, 5, len = 10), nfolds = 5, niter = 5, v = NULL, trace = TRUE, orth = FALSE, center = TRUE, vpos = FALSE, vneg = FALSE) {
  if (vpos && vneg) stop("Cannot constrain elements of v to be both positive and negative.")
  if (nfolds < 2) stop("Must run at least 2 cross-validation folds.")
  percentRemove <- min(0.25, 1 / nfolds)
  call <- match.call()
  if (max(sumabsvs) > sqrt(ncol(x)) || min(sumabsvs) < 1) stop("sumabs must be between 1 and sqrt(ncol(x))")
  xfill <- x
  missing <- is.na(x)
  v <- PMA:::CheckPMDV(v, x, K = 1)
  errs <- matrix(NA, nrow = nfolds, ncol = length(sumabsvs))
  nonzerovs <- matrix(NA, nrow = nfolds, ncol = length(sumabsvs))
  rands <- matrix(runif(nrow(x) * ncol(x)), ncol = ncol(x))
  for (i in 1 : nfolds) {
    if (trace) cat(" Fold ", i, " out of ", nfolds, "\n")
    rm <- ((i - 1) * percentRemove < rands) & (rands < i * percentRemove)
    xrm <- x
    xrm[rm] <- NA
    for (j in 1 : length(sumabsvs)) {
      out <- new_SPC(xrm, sumabsv = sumabsvs[j], orth = orth, niter = niter, v = v, trace = FALSE, center = center, K = 1, vpos = vpos, vneg = vneg)
      xhat <- as.numeric(out$d) * out$u %*% t(out$v)
      errs[i, j] <- sum(((xhat - x)[rm & !missing]) ^ 2)
      nonzerovs[i, j] <- sum(out$v != 0)
    }
  }
  if (trace) cat(fill = TRUE)
  err.means <- apply(errs, 2, mean)
  err.sds <- apply(errs, 2, sd) / sqrt(nfolds)
  nonzerovs.mean <- apply(nonzerovs, 2, mean)
  bestsumabsv <- sumabsvs[which.min(err.means)]
  bestsumabsv1se <- sumabsvs[min(which(err.means < min(err.means) + err.sds[which.min(err.means)]))]
  object <- (list(cv = err.means, cv.error = err.sds, bestsumabsv = bestsumabsv, nonzerovs = nonzerovs.mean, v.init = v, call = call, sumabsvs = sumabsvs, nfolds = nfolds, bestsumabsv1se = bestsumabsv1se, vpos = vpos, vneg = vneg))
  class(object) <- "SPC.cv"
  return(object)
}






