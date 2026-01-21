







print_link <- function(df, type, op) {

  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  # 2. Affichage de l'en-tête
  cat(type, ":\n")
  cat(sprintf("%-14s %8s %8s %8s %8s\n",
              "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))

  latents <- unique(df$lhs)

  for (lat in latents) {
    # Affichage du nom de la variable latente (ex: CUSTOMER_E =~)
    cat(paste0("  ", lat, " ", op, "\n"))

    # On prend les lignes correspondant à cette variable
    group <- df[df$lhs == lat, ]

    for (i in 1:nrow(group)) {
      row <- group[i, ]

      # Nom de l'indicateur (rhs)
      rhs <- row$rhs

      # Estimation formatée à 3 décimales
      est <- sprintf("%.3f", row$est)

      # Gestion des paramètres fixés : si SE est 0 ou NA, on laisse vide
      if (is.na(row$se) || row$se == 0) {
        se <- ""
        z_val <- ""
        p_val <- ""
      } else {
        se <- sprintf("%.3f", row$se)
        z_val <- sprintf("%.3f", row$z)
        p_val <- sprintf("%.3f", row$pvalue)
      }
      cat(sprintf("    %-10s %8s %8s %8s %8s\n",
                  rhs, est, se, z_val, p_val))
    }
  }
}
print_variance <- function(df, type) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  cat("\n", type, ":\n")
  cat(sprintf("%-14s %8s %8s %8s %8s\n",
              "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))

  for (i in 1:nrow(df)) {
    row <- df[i, ]
    var_name <- row$lhs
    display_name <- paste0("   .", var_name) # Ajout du point
    # --- Formatage des chiffres ---
    est <- sprintf("%.3f", row$est)

    # Gestion des cas où SE est vide/zéro (paramètres fixés)
    if (is.na(row$se) || row$se == 0) {
      se <- ""
      z_val <- ""
      p_val <- ""
    } else {
      se <- sprintf("%.3f", row$se)
      z_val <- sprintf("%.3f", row$z)
      p_val <- sprintf("%.3f", row$pvalue)
    }

    cat(sprintf("%-14s %8s %8s %8s %8s\n",
                display_name, est, se, z_val, p_val))
  }
}