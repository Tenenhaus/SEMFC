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

cors <- sapply(1:300, function(i) abs(cor(rocsvd$lambdas[i, ], rocsggcca$lambdas[i, ])))



sgcca_opt <- rgcca(Y_2, sparsity = c(optrgcca[[1]],1 ,1,1,1,1))
class_rgcca_opt <- classification(sgcca_opt$a[[1]])

ssvd_opt <- sparse_svd(Y, c(1,0,0,0,0,0), c(optsvd[[1]],0,0,0,0,0))
res_svd <- classification(ssvd_opt[[1]])


# val_diff <- rocsggcca$tab$val[ which(abs(rocsggcca$tab$f1 - rocsvd$tab$f1) != 0)]
#
# class_diff_svd <- rocsvd$tab[ which(abs(rocsggcca$tab$f1 - rocsvd$tab$f1) != 0), ]
# class_diff_sgcca <- rocsggcca$tab[ which(abs(rocsggcca$tab$f1 - rocsvd$tab$f1) != 0), ]
#
# rownames(class_diff_svd) <- paste0('svd', seq_len(nrow(class_diff_svd)))
# rownames(class_diff_sgcca) <- paste0('sgcca', seq_len(nrow(class_diff_sgcca)))

plot_compare_metric(rocsvd$tab, rocsggcca$tab, 'f1')




sgcca_diff <- rgcca(Y_2, sparsity = c(ind,1 ,1,1,1,1), tol = 1e-07, verbose = TRUE)
plot(sgcca_diff$a[[1]])

svd_diff <- sparse_svd(Y, c(1,0,0,0,0,0), c(ind*sqrt(180),0,0,0,0,0), trace = TRUE, v = mon_vecteur)
plot(svd_diff[[1]])

plot(abs(sgcca_diff$a[[1]] - svd_diff[[1]]))


modelrgcca <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='rgcca',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=0.5635770)

modelrgcca$fit_svd()

modelpmd <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='pmd',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=7.0628614)
modelpmd$fit_svd()


tls_sgcca <- comp_tls(300, 'sgcca')

tls_pmd <- comp_tls(300, 'pmd')

plot_compare_metric(tls_pmd, tls_sgcca, 'theorical_ls')
plot_compare_metric(tls_pmd, tls_sgcca, 'empirical_ls')





source('newpmd/pmd_coordonate_descente.R')


lambdas_newpmd <- newpmd(Y, C, mode, c(7.56,0.632 * sqrt(100),sqrt(3),sqrt(3),sqrt(3),sqrt(3)),
                         niter=50, tol=1e-07)

time_rocnewsvd <- system.time(rocnewsvd <- roc(300, 'newpmd'))

lambdas_test <- newpmd(Y, C, mode, c( 0.1488204,sqrt(3),sqrt(3),sqrt(3),sqrt(3),sqrt(3)),
                         niter=50, tol=1e-07)



plot_compare_metric(rocnewsvd$tab, rocsggcca$tab, 'f1')
plot_compare_metric(rocnewsvd$tab, rocsggcca$tab, 'balanced_accuracy')



