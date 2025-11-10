


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



find_max_sparsity_svd <-function(len){
  tab<- matrix(nrow=0, ncol = 2)
  colnames(tab) <- c("val","f1")


  for (val in seq(1, sqrt(ncol(Y[[1]])), length = len)){
    ssvd <- sparse_svd(Y, c(1,0,0,0,0,0), c(val,0,0,0,0,0))
    res_svd <- classification(ssvd[[1]])
    f1 <- res_svd[[3]]
    tab <- rbind(tab,c(val, f1))


  }

  f1max <- max(tab[,2])
  imax <- which.max(tab[,2])
  valmax <- tab[imax, 1]

  return(c(valmax,f1max))


}


find_max_sparsity_rgcca <-function(len){
  tab<- matrix(nrow=0, ncol = 2)
  colnames(tab) <- c("val","f1")


  for (val in seq(1/sqrt(ncol(Y[[1]])), 1, length = len)){
    sgcca <- rgcca(Y_2, sparsity = c(val,1 ,1,1,1,1))
    res_sgcca <- classification(sgcca$a[[1]])
    f1 <- res_sgcca[[3]]
    tab <- rbind(tab,c(val, f1))


  }

  f1max <- max(tab[,2])
  imax <- which.max(tab[,2])
  valmax <- tab[imax, 1]

  return(c(valmax,f1max))


}



roc <- function(len, method){
  first_cols <- c("val", "recall", "fpr", "f1")
  other_metrics <- c("precision", "FNR", "specificity", "accuracy", "balanced_accuracy", "youden")
  all_cols <- c(first_cols, other_metrics)

  tab <- matrix(nrow=0, ncol = length(all_cols))
  colnames(tab) <- all_cols




  for (val in seq(1/sqrt(ncol(Y[[1]])), 1, length = len)){


      if (method=='sgcca'){
        fit <- rgcca(Y_2, sparsity = c(val,1 ,1,1,1,1))
        serie <- fit$a[[1]]
      }
      else if (method=='pmd'){
        val_pmd <- val * sqrt(ncol(Y[[1]]))
        fit <- sparse_svd(Y, c(1,0,0,0,0,0), c(val_pmd,0,0,0,0,0))
        serie <- fit[[1]]

      }

      res <- classification(serie)
      vals <- c(val,
                res$recall, res$FPR, res$f1,
                res$precision, res$FNR, res$specificity, res$accuracy, res$balanced_accuracy, res$youden)
      tab <- rbind(tab, vals)

    }

  tab <- as.data.frame(tab)

  # Tri par FPR croissant (nécessaire pour l’intégration)
  tab_roc <- tab[order(tab$fpr), ]

  # Calcul de l’AUC via la règle du trapèze
  auc <- sum(diff(tab_roc$fpr) * (head(tab_roc$recall, -1) + tail(tab_roc$recall, -1)) / 2)


  # Tracé de la ROC
  plot(tab_roc$fpr, tab_roc$recall, type = "l", col = "blue", lwd = 2,
         xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = paste("ROC Curve", method))
  abline(0, 1, col = "red", lty = 2)  # Ligne diagonale

  return(list(tab = tab, auc = auc))

}
