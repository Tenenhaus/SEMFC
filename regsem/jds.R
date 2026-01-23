library(RGCCA)
library(PMA)



source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')



for (N in c(300, 500, 1000)){
  source('data/data_generated_reflective.R')

  optsvd <- find_max_sparsity_svd(300)
  optrgcca <- find_max_sparsity_rgcca(300)

  rocsggcca <- roc(300, 'sgcca')

  rocsvd <- roc(300, 'pmd')








}
