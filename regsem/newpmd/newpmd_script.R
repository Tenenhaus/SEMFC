library(RGCCA)
library(PMA)


source('data/data_generated_reflective.R')
source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')
source('regsem/newpmd/pmd_coordonate_descente.R')


X <- X_2
Y <- Y_2

model <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen=F,
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=NULL)
model$fit_svd()

rocsggcca <- roc2(300, 'sgcca')
rocnewsvd <- roc2(300, 'newpmd')


plot_compare_metric(rocnewsvd$tab, rocsggcca$tab, 'f1')
plot_compare_metric(rocnewsvd$tab, rocsggcca$tab, 'balanced_accuracy')




lambdas_newpmd <- newpmd(Y, C, mode, c(7.56,0.632 * sqrt(100),sqrt(3),sqrt(3),sqrt(3),sqrt(3)),
                         niter=50, tol=1e-07)



modelrgcca <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='rgcca',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=0.5635770)

modelrgcca$fit_svd()

sgcca_opt <- rgcca(Y_2, sparsity = c(0.5635770,0.635 ,1,1,1,1))
plot(sgcca_opt$a[[2]])