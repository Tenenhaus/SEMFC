library(RGCCA)
library(PMA)



source('R/SEMFC/sem_f_c.R')
source('regsem/utils_function_sparse.R')


N = 1000



source('data/data_generated_reflective.R')
X <- X_2
Y <- Y_2


L <- lapply(Y,scale)
L <- lapply(L, function (x) x/sqrt(NCOL(x)))


data <- t(t(L[[1]])%*%Reduce("cbind", L[-1]))



perm_out = rgcca_permutation(list(data, data), method = "spls")
plot(perm_out)
summary(perm_out)




out <- SPC(data, sumabsv=perm_out$best_params[[1]] *sqrt(ncol(Y[[1]])), K=1, center = FALSE, , niter = 1000)$v


plot(out)