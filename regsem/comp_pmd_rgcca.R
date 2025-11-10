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
                      sparse_val=7.0628614)
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
                      sparse_val=0.5635770)

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



plot(table_lambda_empirical$penrgccacv)









sparse_svd.cv <- function(L, pen, len_seq, nfold, niter){

  res <-  list()

  for (x in 1:length(L)){
    if (pen[[x]] == 1){
      # data <- t(L[[x]])%*%Reduce("cbind", L[-x])%*%t(t(L[[x]])%*%Reduce("cbind", L[-x]))
      data <- t(t(L[[x]])%*%Reduce("cbind", L[-x]))
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
  L = lapply(L,scale)
  L = lapply(L, function (x) x/sqrt(NCOL(x)))

  for (x in 1:length(L)){
    if (pen[[x]] == 1){
      data <- t(t(L[[x]])%*%Reduce("cbind", L[-x]))

      out <- SPC(data, sumabsv=values[[x]], K=1, center = FALSE)$v

      res[[x]] <- out

    }
  }

  return (res)

}




sparse.svd = sparse_svd(Y_2, c(1,1,0,0,0,0), c(optsvd[[1]],sqrt(3),0,0,0,0))
plot(-sparse.svd[[1]])

sparse.svd.cv = sparse_svd.cv(Y_2, c(1,0,0,0,0,0), 50,20,30)
plot(sparse.svd.cv[[1]])


perm_out = rgcca_permutation(Y_2, scheme = "factorial", par_type = "sparsity",
                             par_value = cbind(seq(0.1, 1, length = 70), 1, 1, 1, 1, 1), n_perms = 30)
rgcca_final = rgcca(perm_out)
plot(rgcca_final$a[[1]])


yy = rgcca(Y_2, sparsity = c(optrgcca[[1]],1 ,1,1,1,1))
plot(yy$a[[1]])

clas_opt_svd = unlist(classification(sparse.svd[[1]]))
clas_cv_svd = unlist(classification(sparse.svd.cv[[1]]))
class_opt_rgcca = unlist(classification(yy$a[[1]]))
class_perm_rgcca = unlist(classification(rgcca_final$a[[1]]))


performance = rbind(clas_opt_svd, clas_cv_svd, class_opt_rgcca,class_perm_rgcca)

classification <- function(pred){
  true_val <- as.integer(l1 !=0)
  pred_val <- as.integer(pred !=0)
  true_val <- factor(true_val, levels = c(0, 1))
  pred_val <- factor(pred_val, levels = c(0, 1))

  matrice_confusion <- table(Observation = true_val, Prediction = pred_val)

  TP <- matrice_confusion["1","1"]
  FP <- matrice_confusion["0","1"]
  FN <- matrice_confusion["1","0"]
  TN <- matrice_confusion["0","0"]

  total <- TP + FP + TN + FN
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  balanced_accuracy <- (recall + specificity) / 2
  accuracy <- (TP + TN) / total
  FPR <- FP / (FP + TN)
  FNR <- 1 - recall
  f1 <- 2 * (precision * recall) / (precision + recall)
  youden <- recall + specificity - 1

  return (list(
    precision = precision,
    recall = recall,
    f1 = f1,
    FPR = FPR,
    FNR = FNR,
    specificity = specificity,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    youden = youden
  ))


}






n_rep = 10

tab_svd <- matrix(NA, nrow = n_rep, ncol = 3)
colnames(tab_svd) <- c("precision", "recall", "f1")

tab_rgcca <- matrix(NA, nrow = n_rep, ncol = 4)
colnames(tab_rgcca) <- c("precision", "recall", "f1", 'val')




for (i in 1:n_rep){
  print('svd iteration :')
  print(i)

  sparse.svd.cv = sparse_svd.cv(Y, c(1,0,0,0,0,0), 50,20,30)
  res_svd = classification(sparse.svd.cv[[1]])
  tab_svd[i,] <- unlist(res_svd)

  print('rgcca iteration :')
  print(i)



  perm_out = rgcca_permutation(Y, scheme = "factorial", par_type = "sparsity",
                               par_value = cbind(seq(0.1, 1, length = 70), 1, 1, 1, 1, 1), n_perms = 30)
  rgcca_final = rgcca(perm_out)

  res_rgcca <- classification(rgcca_final$a[[1]])
  tab_rgcca[i,] <- c(unlist(res_rgcca), perm_out$best_params[[1]])




}



find_max_sparsity_svd <-function(len){
  tab<- matrix(nrow=0, ncol = 2)
  colnames(tab) <- c("val","f1")


  for (val in seq(1, sqrt(ncol(Y[[1]])), length = len)){
    ssvd = sparse_svd(Y, c(1,0,0,0,0,0), c(val,0,0,0,0,0))
    res_svd = classification(ssvd[[1]])
    f1 = res_svd[[3]]
    tab = rbind(tab,c(val, f1))


  }

  f1max = max(tab[,2])
  imax= which.max(tab[,2])
  valmax= tab[imax, 1]

  return(c(valmax,f1max))


}


find_max_sparsity_rgcca <-function(len){
  tab<- matrix(nrow=0, ncol = 2)
  colnames(tab) <- c("val","f1")


  for (val in seq(1/sqrt(ncol(Y[[1]])), 1, length = len)){
    sgcca = rgcca(Y_2, sparsity = c(val,1 ,1,1,1,1))
    res_sgcca = classification(sgcca$a[[1]])
    f1 = res_sgcca[[3]]
    tab = rbind(tab,c(val, f1))


  }

  f1max = max(tab[,2])
  imax= which.max(tab[,2])
  valmax= tab[imax, 1]

  return(c(valmax,f1max))


}


optsvd = find_max_sparsity_svd(300)
optrgcca = find_max_sparsity_rgcca(300)




roc = function(len, method){
  tab<- matrix(nrow=0, ncol = 4)
  colnames(tab) <- c("val","recall", "fpr", "f1")


  for (val in seq(1/sqrt(ncol(Y[[1]])), 1, length = len)){


      if (method=='sgcca'){
        fit = rgcca(Y_2, sparsity = c(val,1 ,1,1,1,1))
        serie = fit$a[[1]]
      }
      else if (method=='pmd'){
        val_pmd = val * sqrt(ncol(Y[[1]]))
        fit = sparse_svd(Y, c(1,0,0,0,0,0), c(val_pmd,0,0,0,0,0))
        serie = fit[[1]]

      }

      res = classification(serie)
      recall = res[[2]]
      fpr = res[[4]]
      f1 = res[[3]]
      tab = rbind(tab,c(val,recall, fpr, f1))


    }

  tab <- as.data.frame(tab)

  # Tri par FPR croissant (nécessaire pour l’intégration)
  tab <- tab[order(tab$fpr), ]

  # Calcul de l’AUC via la règle du trapèze
  auc <- sum(diff(tab$fpr) * (head(tab$recall, -1) + tail(tab$recall, -1)) / 2)


  # Tracé de la ROC
  plot(tab$fpr, tab$recall, type = "l", col = "blue", lwd = 2,
         xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = paste("ROC Curve", method))
  abline(0, 1, col = "red", lty = 2)  # Ligne diagonale

  return(list(tab = tab, auc = auc))

}

rocsggcca = roc(300, 'sgcca')
rocsvd = roc(300, 'pmd')

time_rocrgcca <- system.time(rocsggcca <- roc(300, 'sgcca'))
print(time_rocrgcca)

rocsggcca <- roc(300, 'sgcca')

time_rocsvd <- system.time(rocsvd <- roc(300, 'pmd'))
print(time_rocsvd)


sgcca_opt = rgcca(Y_2, sparsity = c(optrgcca[[1]],1 ,1,1,1,1))
class_rgcca_opt = classification(sgcca_opt$a[[1]])

ssvd_opt = sparse_svd(Y, c(1,0,0,0,0,0), c(optsvd[[1]],0,0,0,0,0))
res_svd = classification(ssvd_opt[[1]])


val_diff = rocsggcca$tab$val[ which(abs(rocsggcca$tab$f1 - rocsvd$tab$f1) != 0)]
ind = val_diff[1]

sgcca_diff = rgcca(Y_2, sparsity = c(ind,1 ,1,1,1,1))
plot(sgcca_diff$a[[1]])

svd_diff = sparse_svd(Y, c(1,0,0,0,0,0), c(ind*sqrt(180),0,0,0,0,0))
plot(svd_diff[[1]])

plot(abs(sgcca_diff$a[[1]] - svd_diff[[1]]))


for (val in rocsggcca$tab$val){
  sgcca_diff = rgcca(Y_2, sparsity = c(val,1 ,1,1,1,1))
  plot(sgcca_diff$a[[1]])
}