#' Dynamic Network Visualization for fitted RGCCAC object 
#' 
#' Gives a dynamic graph where you can change the colors and the arrangement 
#' of the nodes 
#'
#' @param fit a fitted rgccac object 
#' @param color.background string's vector of the color of the node's 
#' background of the latents and manifests variables.
#' Default to c("white", "white")
#' @param color.border string's vector of the color of the 
#' node's border of the latents and manifests variables.
#' Default to c("black", "black")
#' @param shape string's vector of the shape of the nodes
#' latents variables and manifests variables.
#' Default to c("oval", "box")
#' @param color.edge string's vector of the color of the 
#' edges between LV and MV and between LV.
#' Default to c("black", "black")
#' @param display_numbers boolean if display_numbers = TRUE, 
#' display the cor_y_eta between latents and manifests nodes 
#' and the betas/gammas between the latents nodes
#' @param configure boolean If true, 
#' open a control panel that allows you to modify some parameters. 
#' 
#'
#' @return dynamic graph
#' @export
#'
#' @examples
#' VisNetworkRGCCA(fit.rgccac,color.background = c("red", "blue"), 
#' color.border = c("red","blue"), shape = c("oval", "box"), 
#' color.edge = c("black","black"), display_numbers = T, configure = T)



VisNetworkRGCCA = function(fit, color.background = c("white", "white"),
                         color.border = c("black", "black"), 
                         shape = c("oval", "box"), 
                         color.edge = c("black", "black"),
                         display_numbers = T, configure = T){
  
  #tau_rep  = rep(fit$call$tau, sapply(fit$a, length))
  tau_rep = rep(1,sum(sapply(fit$cor_y_eta, length)))
  name_block = names(fit$d) 
  # Give the "from" and the "to" of the latent variables
  inner = matrix(name_block[which(fit$C !=0, arr.ind = TRUE)], ncol = 2)
  mat_gamma_beta = cbind(fit$gamma, fit$beta)
  idf = sapply(rownames(fit$Phat), function(x) {
    which(colnames(mat_gamma_beta) == x)})
  mat_gamma_beta <- mat_gamma_beta[ , idf]
  mat2 = t(mat_gamma_beta)
  Cb = fit$C
  Cb[Cb == 1] = mat2[mat2 != 0]
  Cb
  corr_results = round(Cb[which(Cb != 0, arr.ind = TRUE)], 2)
  
  
  #Correlation between connected LVs
  #corr_matrix= cor(Reduce ("cbind", fit$Y))*C
  
  #corr_results=corr_matrix[which(C!=0, arr.ind = TRUE)]
  
  fit_cor_y_eta = as.vector(as.vector(Reduce("c", fit$cor_y_eta)))
  fusion = round(c(corr_results,fit_cor_y_eta), 2)
  # We round to 2 digits after the comma for more readability and we give 
  # the name "label" so that it appears on the edges
  
  
  ManifVar = as.vector(unlist(sapply(fit$cor_y_eta,rownames)))
  nomdf = c(name_block, ManifVar)
  
  #On rajoute une colonne pour la couleur et une autre pour la forme
  #Add a column for the color and another one for the shape
  df = data.frame(id = nomdf,
                label = nomdf,
                shape = rep(shape, c(length(name_block), length(ManifVar))),
                color.background = rep(color.background, 
                                       c(length(name_block), length(ManifVar))),
                color.border = rep(color.border, 
                                   c(length(name_block), length(ManifVar)))
  )
  
  #names of the "manifests" variables
  MV = Reduce("c", sapply(fit$cor_y_eta, rownames))
  LV = rep(rownames(fit$Phat), sapply(fit$cor_y_eta, length))
  MVcopy = MV
  LVcopy = LV
  for (x in which(tau_rep == 1)){
    LV[x] = MVcopy[x]
    MV[x] = LVcopy[x]
  }
  # "from" nodes
  from = c(inner[, 1], MV)
  # "to" nodes
  to = c(inner[, 2], LV)
  color.edge=rep(color.edge,c(length(inner[ , 1]), length(MV)))
  # data.frame which make the link between the nodes with optional labels
  #Smooth = T prevents overlapping of the labels on the edges
  if(display_numbers)
  link <- data.frame(from = from  , to = to, 
                       color = color.edge, label = fusion, smooth = TRUE) else 
  link <- data.frame(from = from, color = color.edge, to = to, smooth = FALSE)
  
  
  visNetwork(df, link, width = "100%")%>%
    visEdges(arrows="to",smooth = T)%>%
    visInteraction(multiselect = T)%>%
    visConfigure(enabled = configure)%>%
    visPhysics(enabled = T)
}