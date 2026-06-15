# ====================================================================
# LE TEST
# ====================================================================

# 1. Configuration du système
P <- 4
block_sizes <- c(2, 2)
j <- 1
p_j <- 2

# 2. On génère une matrice S symétrique (4x4) et un vecteur a (2x1) aléatoires
set.seed(42)
S <- matrix(rnorm(16), P, P)
S[lower.tri(S)] <- t(S)[lower.tri(S)] # On force la symétrie
a <- rnorm(p_j)

# --- PARTIE A : CALCUL THÉORIQUE MANUEL ---

# On extrait les variables de la matrice pour nos formules
a1 <- a[1] ; a2 <- a[2]
s31 <- S[3,1] ; s41 <- S[4,1]
s32 <- S[3,2] ; s42 <- S[4,2]

# 1. Construction manuelle de J_z (2 lignes, 10 colonnes uniques)
J_z_theorique <- matrix(0, nrow = 2, ncol = 10)

# Remplissage de la Ligne 1 (pour z1)
J_z_theorique[1, 3] <- 2*a1*s31 + a2*s32   # Colonne s31
J_z_theorique[1, 4] <- 2*a1*s41 + a2*s42   # Colonne s41
J_z_theorique[1, 6] <- a2*s31              # Colonne s32
J_z_theorique[1, 7] <- a2*s41              # Colonne s42

# Remplissage de la Ligne 2 (pour z2)
J_z_theorique[2, 3] <- a1*s32              # Colonne s31
J_z_theorique[2, 4] <- a1*s42              # Colonne s41
J_z_theorique[2, 6] <- a1*s31 + 2*a2*s32   # Colonne s32
J_z_theorique[2, 7] <- a1*s41 + 2*a2*s42   # Colonne s42

# 2. Calcul manuel du vecteur z
z1 <- (s31^2 + s41^2)*a1 + (s31*s32 + s41*s42)*a2
z2 <- (s31*s32 + s41*s42)*a1 + (s32^2 + s42^2)*a2
z_theorique <- c(z1, z2)

# 3. Calcul manuel de la matrice de normalisation N_z
norme_z <- sqrt(sum(z_theorique^2))
N_z_theorique <- (diag(2) - (z_theorique %*% t(z_theorique)) / (norme_z^2)) / norme_z

# 4. Jacobienne Finale Théorique
J_finale_theorique <- N_z_theorique %*% J_z_theorique


# --- PARTIE B : CALCUL AVEC NOTRE CODE OPTIMISÉ ---

# On génère la matrice de duplication une seule fois
D_P_cache <- generer_Dp(P)

# On lance la grande machine vectorisée
J_finale_code <- jacobienne_finale_vech(block_sizes, j, S, a, D_P_cache)


# --- PARTIE C : LA COMPARAISON ---

cat("--- RÉSULTATS DU TEST ---\n")
# Arrondi pour l'affichage visuel
cat("\nJacobienne Théorique (les 7 premières colonnes) :\n")
print(round(J_finale_theorique[, 1:7], 4))

cat("\nJacobienne du Code (les 7 premières colonnes) :\n")
print(round(as.matrix(J_finale_code)[, 1:7], 4))

# Le test absolu : les matrices sont-elles numériquement identiques ?
# (La tolérance par défaut gère les micro-erreurs d'arrondi des processeurs)
sont_identiques <- all.equal(as.matrix(J_finale_theorique), as.matrix(J_finale_code))
cat("\nLe code et la théorie donnent-ils le même résultat ? ->", sont_identiques, "\n")