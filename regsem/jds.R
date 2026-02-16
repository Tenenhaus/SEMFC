library(RGCCA)
library(PMA)



source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')


Ns <- c(300, 500, 1000)

table_AUC <- matrix(NA, nrow = length(Ns), ncol = 2)
colnames(table_AUC) <- c('pmd', 'sgcca')
rownames(table_AUC) <- paste0('N=', Ns)



table_nul_FP <- matrix(NA, nrow = length(Ns), ncol = 4)


table_fort_FN <- matrix(NA, nrow = length(Ns), ncol = 4)

table_faible_FN <- matrix(NA, nrow = length(Ns), ncol = 4)


faible_range <- 41:60
fort_range <- 1:20
nul_range <- 21:40

# faible_range <- 131:160
# fort_range <- 1:30
# nul_range <- 31:130


colnames(table_nul_FP) <- colnames(table_fort_FN) <- colnames(table_faible_FN) <-
      c('pmd_optimal', 'sgcca_optimal', 'pmd_cv', 'sgcca_permutation')

rownames(table_nul_FP) <- rownames(table_fort_FN) <- rownames(table_faible_FN)  <- paste0('N=', Ns)



table_sparcity<- matrix(NA, nrow = length(Ns), ncol = 4)
colnames(table_sparcity) <- c('pmd_optimal', 'sgcca_optimal', 'pmd_cv', 'sgcca_permutation')
rownames(table_sparcity) <- paste0('N=', Ns)


for (N in Ns){
  print(N)
  source('data/data_generated_reflective.R')
  X <- X_2
  Y <- Y_2

  optsvd <- find_max_sparsity_svd(300)
  optrgcca <- find_max_sparsity_rgcca(300)

  rocsggcca <- roc(300, 'sgcca')
  rocsvd <- roc(300, 'pmd')


  table_AUC[which(Ns == N), 'pmd'] <-  rocsvd$auc
  table_AUC[which(Ns == N), 'sgcca'] <-  rocsggcca$auc

  sgcca_opt <- rgcca(Y_2, sparsity = c(optrgcca[[1]],1 ,1,1,1,1))
  ssvd_opt <- sparse_svd(Y, c(1,0,0,0,0,0), c(optsvd[[1]],0,0,0,0,0))


  sparse.svd.cv = sparse_svd.cv(Y_2, c(1,0,0,0,0,0), 90,30,40)



  perm_out = rgcca_permutation(Y_2, scheme = "factorial", par_type = "sparsity",
                             par_value = cbind(seq((1/sqrt(ncol(Y[[1]]))+0.01), 1, length = 70), 1, 1, 1, 1, 1), n_perms = 30)
  # perm_out = rgcca_permutation(Y_2, scheme = "factorial", par_type = "sparsity",
  #                              par_value = cbind(seq(0.1, 1, length = 70), 1, 1, 1, 1, 1), n_perms = 30)

  rgcca_final = rgcca(perm_out)

  fnr_faible_svd_opt = classification(ssvd_opt[[1]][faible_range], l1[faible_range])$FNR
  fnr_fort_svd_opt = classification(ssvd_opt[[1]][fort_range], l1[fort_range])$FNR
  fpr_nul_svd_opt = classification(ssvd_opt[[1]][nul_range], l1[nul_range])$FPR

  fnr_faible_sgcca_opt = classification(sgcca_opt$a[[1]][faible_range], l1[faible_range])$FNR
  fnr_fort_sgcca_opt = classification(sgcca_opt$a[[1]][fort_range], l1[fort_range])$FNR
  fpr_nul_sgcca_opt = classification(sgcca_opt$a[[1]][nul_range], l1[nul_range])$FPR

  fnr_faible_svd_cv = classification(sparse.svd.cv$a[[1]][faible_range], l1[faible_range])$FNR
  fnr_fort_svd_cv = classification(sparse.svd.cv$a[[1]][fort_range], l1[fort_range])$FNR
  fpr_nul_svd_cv = classification(sparse.svd.cv$a[[1]][nul_range], l1[nul_range])$FPR

  fnr_faible_sgcca_perm = classification(rgcca_final$a[[1]][faible_range], l1[faible_range])$FNR
  fnr_fort_sgcca_perm = classification(rgcca_final$a[[1]][fort_range], l1[fort_range])$FNR
  fpr_nul_sgcca_perm = classification(rgcca_final$a[[1]][nul_range], l1[nul_range])$FPR


  table_fort_FN[which(Ns == N), ] <- c(fnr_fort_svd_opt, fnr_fort_sgcca_opt,
                                         fnr_fort_svd_cv, fnr_fort_sgcca_perm)

  table_faible_FN[which(Ns == N), ] <- c(fnr_faible_svd_opt, fnr_faible_sgcca_opt,
                                           fnr_faible_svd_cv, fnr_faible_sgcca_perm)

  table_nul_FP[which(Ns == N), ] <- c(fpr_nul_svd_opt, fpr_nul_sgcca_opt, fpr_nul_svd_cv, fpr_nul_sgcca_perm)

  table_sparcity[which(Ns == N), ] <- c(optsvd[[1]], optrgcca[[1]]*sqrt(ncol(Y[[1]])),
                                        sparse.svd.cv$param, perm_out$best_params[[1]]*sqrt(ncol(Y[[1]])))


  print(table_AUC)
  print(table_nul_FP)
  print(table_fort_FN)
  print(table_faible_FN)
  print(table_sparcity)

}
