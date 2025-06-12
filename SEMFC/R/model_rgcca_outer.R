#' Outer model RGCCA
#' 
#' 
#'
#' @param fit A fitted rgccac object
#' @param labels  boolean , if labels = TRUE, add the correlation
#' between latents and manifests variables 
#' @param num string or integer. Display the outer graph of the num block or
#' all the blocks if num  = "all"
#'
#' @return
#' @export
#'
#' @examples
#' model_rgcca_outer(fit,labels = TRUE, num = 1)
model_rgcca_outer <- function(fit, labels = labels, num = num){
  
  
  
  fit_cor_y_eta = lapply(as.vector(as.vector(Reduce("c", fit$cor_y_eta))), 
                         round, 2)
  
  inner = matrix(names(A)[which(fit$C != 0, arr.ind = TRUE)], ncol = 2)
  LV = rep(names(A), sapply(fit$cor_y_eta, length))
  MV = Reduce("c", sapply(fit$cor_y_eta, rownames))
  tau_rep = rep(1,sum(sapply(fit$cor_y_eta, length)))
  if (num == "all"){
    
    if (labels == T){
      
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
      
      LM = paste(sapply(1:length(tau_rep),
                        function(x)
                          ifelse(tau_rep[x] == 1,
                                 paste(LV[x], " -> ", MV[x], collapse="   "),
                                 paste(MV[x], " -> ", LV[x], collapse="   ")
                          )
      ), collapse = " ")
    }
  }else
  {
    LV2 = LV[which(LV == names(fit$d[num]))]
    MV2 = MV[which(LV == names(fit$d[num]))]
    fit_cor_y_eta2 = fit_cor_y_eta[which(LV == names(fit$d[num]))]
    tau_rep = rep(1,length(LV2))
    
    
    if (labels == T){
      
      LM = paste(sapply(1:length(tau_rep),
                        function(x)
                          ifelse(tau_rep[x] == 1,
                                 paste(LV2[x], " -> ", MV2[x],"[label = '",
                                       fit_cor_y_eta2[x] , "']", collapse = "   "),
                                 paste(MV2[x], " -> ", LV2[x],"[label = '",
                                       fit_cor_y_eta2[x] , "']", collapse = "   ")
                          )
      ), collapse = " ")
    }
    else {
      
      LM = paste(sapply(1:length(tau_rep),
                        function(x)
                          ifelse(tau_rep[x] == 1,
                                 paste(LV2[x], " -> ", MV2[x], collapse="   "),
                                 paste(MV2[x], " -> ", LV2[x], collapse="   ")
                          )
      ), collapse = " ")
    }
  }
  paste(c(LM), collapse = " ")
  
}