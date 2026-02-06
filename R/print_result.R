


get_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.1)   return(".")
  return("")
}


print_link <- function(df, type, op, standardized = FALSE) {

  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  # 2. Affichage de l'en-tête
  cat(type, ":\n")
  if (standardized) {
    cat(sprintf("%-14s %8s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)", "Std.all"))
  } else {
    cat(sprintf("%-14s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))
  }

  latents <- unique(df$lhs)

  for (lat in latents) {
    # Affichage du nom de la variable latente (ex: CUSTOMER_E =~)
    cat(paste0("  ", lat, " ", op, "\n"))

    # On prend les lignes correspondant à cette variable
    group <- df[df$lhs == lat, ]

    for (i in seq_len(nrow(group))) {
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
        std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
        stars <- ""
      } else {
        se <- sprintf("%.3f", row$se)
        z_val <- sprintf("%.3f", row$z)
        p_val <- sprintf("%.3f", row$pvalue)
        std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
        stars <- get_stars(row$pvalue)
      }
      if (standardized) {
        cat(sprintf("    %-10s %8s %8s %8s %8s %8s %-3s\n",
                    row$rhs, est, se, z_val, p_val, std_all, stars))
      } else {
        cat(sprintf("    %-10s %8s %8s %8s %8s %-3s\n",
                    row$rhs, est, se, z_val, p_val, stars))
      }
    }
  }
}




print_variance <- function(df, type, standardized = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  cat("\n", type, ":\n")
  if (standardized) {
    cat(sprintf("%-14s %8s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)", "Std.all"))
  } else {
    cat(sprintf("%-14s %8s %8s %8s %8s\n",
                "", "Estimate", "Std.Err", "z-value", "P(>|z|)"))
  }

  for (i in seq_len(nrow(df))) {
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
      std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
      stars <- ""
    } else {
      se <- sprintf("%.3f", row$se)
      z_val <- sprintf("%.3f", row$z)
      p_val <- sprintf("%.3f", row$pvalue)
      std_all <- if (standardized) sprintf("%.3f", row$std.all) else ""
      stars <- get_stars(row$pvalue)
    }

    if (standardized) {
      cat(sprintf("%-14s %8s %8s %8s %8s %8s %-3s\n",
                  display_name, est, se, z_val, p_val, std_all, stars))
    } else {
      cat(sprintf("%-14s %8s %8s %8s %8s %-3s\n",
                  display_name, est, se, z_val, p_val, stars))
    }
  }
}