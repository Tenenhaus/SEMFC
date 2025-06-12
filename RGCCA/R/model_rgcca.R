
#' Model RGCCA
#'
#' @param fit A fitted rgccac object
#' @param labels  boolean , if labels = TRUE, add the correlation
#' between latents and manifests variables 
#' and the betas/gammas between the latents nodes
#' 
#' 
#'
#' @return
#' a string which will be used as a basis for the creation of the graph
#' @export
#'
#' @examples
#' model_rgcca(fit,labels = TRUE)


model_rgcca <- function(fit,labels = labels){
  mat_gamma_beta = cbind(fit$gamma, fit$beta)

  
  idf = sapply(names(A), function(x) {
    which(colnames(mat_gamma_beta) == x)})
  mat_gamma_beta = mat_gamma_beta[ , idf]
  mat2 = t(mat_gamma_beta)
  Cb = fit$C
  Cb[Cb == 1] = mat2[mat2 != 0]
  corr_result = round(Cb[which(Cb != 0, arr.ind = TRUE)], 2)
  
  
  fit_cor_y_eta = lapply(as.vector(as.vector(Reduce("c", fit$cor_y_eta))), 
                         round, 2)
  corr_result_round = lapply(corr_result, round, 2)
  inner = matrix(names(A)[which(fit$C != 0, arr.ind = TRUE)], ncol = 2)
  LV = rep(names(A), sapply(fit$cor_y_eta, length))
  MV = Reduce("c", sapply(fit$cor_y_eta, rownames))
  #tau_rep = rep(fit$call$tau, sapply(fit$cor_y_eta, length))
  tau_rep=rep(1,sum(sapply(fit$cor_y_eta, length)))
  if (labels == T){
    INNER = paste(inner[,1]," -> ",inner[,2], "[label = '", corr_result_round,
                "']",collapse = " ")
    LM = paste(sapply(1:length(tau_rep),
                      function(x)
                        ifelse(tau_rep[x] == 1,
                               paste(LV[x], " -> ", MV[x],"[label = '",
                                     fit_cor_y_eta[x] , "']", collapse = "   "),
                               paste(MV[x], " -> ", LV[x],"[label = '",
                                     fit_cor_y_eta[x] , "']", collapse = "   ")
                        )
    ), collapse = " ")
  }
  else {
    INNER = paste(inner[,1]," -> ", inner[ , 2], collapse = " ")
    LM = paste(sapply(1:length(tau_rep),
                      function(x)
                        ifelse(tau_rep[x] == 1,
                               paste(LV[x], " -> ", MV[x], collapse="   "),
                               paste(MV[x], " -> ", LV[x], collapse="   ")
                        )
    ), collapse = " ")
  }
  paste(c(INNER,LM), collapse = " ")
  
}







