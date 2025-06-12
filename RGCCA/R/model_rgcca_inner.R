#' Title
#'
#' @param fit A fitted rgccac object
#' @param labels  boolean , if labels = TRUE, add 
#' the betas/gammas between the latents nodes
#'
#' @return
#' @export
#'
#' @examples
#' model_rgcca_inner(fit,labels = TRUE)
  model_rgcca_inner <- function(fit,labels = labels){
  mat_gamma_beta = cbind(fit$gamma, fit$beta)
  idf = sapply(names(A), function(x) {
    which(colnames(mat_gamma_beta) == x)})
 
   mat_gamma_beta = mat_gamma_beta[ , idf]
  mat2 = t(mat_gamma_beta)
  Cb = fit$C
  Cb[Cb == 1] = mat2[mat2 != 0]
  corr_result = round(Cb[which(Cb != 0, arr.ind = TRUE)], 2)
  corr_result_round = lapply(corr_result, round, 2)
  inner = matrix(names(A)[which(fit$C != 0, arr.ind = TRUE)], ncol = 2)
  
  if (labels == T){
    INNER = paste(inner[,1]," -> ",inner[,2], "[label = '", corr_result_round,
                  "']",collapse = " ")
    
  }
  else {
    INNER = paste(inner[,1]," -> ", inner[ , 2], collapse = " ")
    
  }
  paste(c(INNER), collapse = " ")
  
}