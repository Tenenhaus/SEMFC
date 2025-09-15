library(RGCCA)
library(PMA)


source('data/data_generated_reflective.R')
source('R/SEMFC/sem_f_c.R')

pen_vars <- paste0("X1", 2:180, collapse = " + ")
pen_line <- paste0("pen() * eta1 =~ ", pen_vars)



sem.model <-paste0('
# latent variable definitions
eta1 =~ X11+',pen_vars,'
eta2 =~ X21+X22+X23
eta3 =~ X31+X32+X33
eta4 =~ X41+X42+X43
eta5 =~ X51+X52+X53
eta6 =~ X61+X62+X63

# Regressions
eta5 ~ eta1 + eta2 + eta6
eta6 ~ eta3 + eta4 + eta5

# residual covariances
eta5 ~~ eta6')

X = X_2
Y = Y_2

lav <- sem(sem.model, X)
estimate = parameterEstimates(lav, standardized = TRUE)

#formate lavaan stand lambda
lambda_lavaan  = estimate[estimate$op == "=~",c(1:3,11) ]



modelpmd <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='pmd',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=7)
modelpmd$fit_svd()

pensvd = mapply(function(x,y) sqrt(1-(x/y)), unlist(modelpmd$svd_parameters$residual_variance),diag(modelpmd$cov_S))


modelpmdcv <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='pmd.cv',
                      len_seq = 50,
                      nfold=20,
                      niter=30,
                      sparse_val=NULL)
modelpmdcv$fit_svd()

pensvdcv = mapply(function(x,y) sqrt(1-(x/y)), unlist(modelpmdcv$svd_parameters$residual_variance),diag(modelpmdcv$cov_S))



modelrgcca <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='rgcca',
                      len_seq = NULL,
                      nfold=NULL,
                      niter=NULL,
                      sparse_val=0.57)

modelrgcca$fit_svd()

penrgcca = mapply(function(x,y) sqrt(1-(x/y)), unlist(modelrgcca$svd_parameters$residual_variance),diag(modelrgcca$cov_S))



modelrgccacv <- SemFC$new(data=Y, relation_matrix = C, mode=mode, scale=F, bias=F, pen='rgcca.cv',
                      len_seq = 50,
                      nfold=NULL,
                      niter=10,
                      sparse_val=NULL)
modelrgccacv$fit_svd()
penrgccacv = mapply(function(x,y) sqrt(1-(x/y)), unlist(modelrgccacv$svd_parameters$residual_variance),diag(modelrgccacv$cov_S))




lambda_th = c(l1,l2,l3,l4,l5,l6)

table_lambda_empirical = cbind(lambda_th, lambda_lavaan, pensvd, pensvdcv,penrgcca, penrgccacv)



plot(table_lambda_empirical$pensvdcv)





perm_out = rgcca_permutation(Y_2, scheme = "factorial", par_type = "sparsity",
                             par_value = cbind(seq(0.1, 1, length = 50), 1, 1, 1, 1, 1), n_perms = 10)
rgcca_final = rgcca(perm_out)
plot(rgcca_final$a[[1]])




sparse_svd.cv <- function(L, pen, len_seq, nfold, niter){

  res <-  list()

  for (x in 1:length(L)){
    if (pen[[x]] == 1){
      data <- t(L[[x]])%*%Reduce("cbind", L[-x])%*%t(t(L[[x]])%*%Reduce("cbind", L[-x]))
      cv.out <- SPC.cv(data, sumabsvs=seq(1, sqrt(ncol(data)), len=len_seq), nfold = nfold, niter =niter)
      print(cv.out$bestsumabsv)
      out <- SPC(data, sumabsv=cv.out$bestsumabsv, K=1, v=cv.out$v.init)$v

      res[[x]] <- out

    }

  }

  return (res)

}


sparse_svd <- function(L, pen, values){

  res <-  list()

  for (x in 1:length(L)){
    if (pen[[x]] == 1){
      data <- t(L[[x]])%*%Reduce("cbind", L[-x])%*%t(t(L[[x]])%*%Reduce("cbind", L[-x]))
      out <- SPC(data, sumabsv=values[[x]], K=1)$v

      res[[x]] <- out

    }
  }

  return (res)

}




sparse.svd = sparse_svd(Y_2, c(1,1,0,0,0,0), c(8.3,sqrt(3),0,0,0,0))
plot(sparse.svd[[1]])

sparse.svd.cv = sparse_svd.cv(Y_2, c(1,0,0,0,0,0), 50,20,30)
plot(sparse.svd.cv[[1]])




yy = rgcca(Y_2, sparsity = c(0.6,1 ,1,1,1,1))
plot(yy$a[[1]])







