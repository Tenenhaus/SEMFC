library(RGCCA)
library(PMA)


source('data/data_generated_reflective.R')


pen_vars <- paste0("X1", 2:20, collapse = " + ")
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



perm_out = rgcca_permutation(Y_2, scheme = "factorial", par_type = "sparsity", par_value = cbind(seq(0.1, 1, length = 50), 1, 1, 1, 1, 1), n_perms = 10)
xx = rgcca(perm_out)
plot(xx$a[[1]])




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




res = sparse_svd(Y_2, c(1,1,0,0,0,0), c(7,sqrt(3),0,0,0,0))


res.cv = sparse_svd.cv(Y_2, c(1,0,0,0,0,0), 50,20,30)