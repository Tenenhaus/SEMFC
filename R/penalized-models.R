# Modeles syntaxiques lavaan / regsem / lslx pour des blocs de tailles libres.
# lambda : liste des vrais loadings ; seuls lengths(lambda) sont utilises ici.

regssem_model <- function(lambda, lslx_identification = "marker") {
  p <- lengths(lambda)
  J <- length(p)
  obs <- lapply(seq_len(J), function(j) paste0("X", j, seq_len(p[j])))
  names(obs) <- paste0("eta", seq_len(J))

  loading_lines <- character(J)
  regsem_lines <- character(J)
  lslx_lines <- character()
  pen_labels <- character()

  for (j in seq_len(J)) {
    eta <- paste0("eta", j)
    x <- obs[[j]]
    labels <- if (p[j] > 1) paste0("b_", j, "_", seq.int(2, p[j])) else character(0)

    # SEM classique : premier loading fixe a 1 par lavaan.
    loading_lines[j] <- paste0(eta, " =~ ", paste(x, collapse = " + "))

    # regsem : le premier loading fixe a 1, les autres ont un label unique.
    regsem_lines[j] <- paste0(
      eta, " =~ 1*", x[1],
      if (length(labels)) paste0(" + ", paste(paste0(labels, "*", x[-1]), collapse = " + "))
    )
    pen_labels <- c(pen_labels, labels)

    # Avec "marker" : meme identification regsem / lslx.
    # Avec "variance" : specification originale de ton article pour lslx.
    first <- if (lslx_identification == "marker") paste0("1*", x[1]) else x[1]
    lslx_lines <- c(lslx_lines, paste0(eta, " =~ ", first))
    if (length(labels)) {
      lslx_lines <- c(lslx_lines,
                      paste0("pen() * ", eta, " =~ ", paste(x[-1], collapse = " + ")))
    }
  }

  # Meme modele structurel que dans ton ancien benchmark (6 variables latentes).
  structure <- c(
    "eta5 ~ eta1 + eta2 + eta6",
    "eta6 ~ eta3 + eta4 + eta5",
    "eta5 ~~ eta6"
  )
  variances_lslx <- if (lslx_identification == "variance") {
    paste0("eta", seq_len(J), " ~~ 1*eta", seq_len(J))
  } else character(0)

  list(
    sem.model = paste(c(loading_lines, structure), collapse = "\n"),
    sem.model.lslx = paste(c(lslx_lines, structure, variances_lslx), collapse = "\n"),
    sem.model.lsem = paste(c(loading_lines, structure), collapse = "\n"),
    sem.model.regsem = paste(c(regsem_lines, structure), collapse = "\n"),
    pen_labels = pen_labels,
    indicators = obs
  )
}
