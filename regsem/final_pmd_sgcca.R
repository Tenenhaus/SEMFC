library(RGCCA)
library(PMA)


source('data/data_generated_reflective.R')
source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')

X <- X_2
Y <- Y_2



optsvd <- find_max_sparsity_svd(300)
optrgcca <- find_max_sparsity_rgcca(300)



time_rocrgcca <- system.time(rocsggcca <- roc(300, 'sgcca'))
mean_sgcca <- colMeans(rocsggcca$tab)
print(time_rocrgcca)

time_rocsvd <- system.time(rocsvd <- roc(300, 'pmd'))
mean_svd <- colMeans(rocsvd$tab)
print(time_rocsvd)


sgcca_opt <- rgcca(Y_2, sparsity = c(optrgcca[[1]],1 ,1,1,1,1))
class_rgcca_opt <- classification(sgcca_opt$a[[1]])

ssvd_opt <- sparse_svd(Y, c(1,0,0,0,0,0), c(optsvd[[1]],0,0,0,0,0))
res_svd <- classification(ssvd_opt[[1]])


