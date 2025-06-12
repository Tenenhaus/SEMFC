#'  Dynamic Inner Network Visualization for fitted RGCCAC object 
#'
#'Gives a dynamic inner graph where you can change the colors and the arrangement 
#' of the nodes 
#'
#'  @param fit a fitted rgccac object 
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
#' @param display_numbers boolean. if display_numbers = TRUE, 
#' display the cor_y_eta between latents and manifests nodes 
#' and the betas/gammas between the latents nodes
#' @param configure boolean. If true, 
#' open a control panel that allows you to modify some parameters. 
#'
#' @return
#' inner dynamic graph
#' @export
#'
#' @examples
#' innerplot(fit.rgccac, color.background = "white",
#'color.border = "black", 
#'shape = "oval", 
#'color.edge = "black",
#'display_numbers = TRUE, configure = TRUE)
innerplot = function(fit, color.background = "white",
                     color.border = "black", 
                     shape = "oval", 
                     color.edge = "black",
                     display_numbers = TRUE, configure = FALSE){
  
  name_block = names(fit$d) 
  inner = matrix(name_block[which(C !=0, arr.ind = TRUE)], ncol = 2)
  mat_gamma_beta = cbind(fit$gamma, fit$beta)
  idf=sapply(rownames(fit$Phat), function(x) {
    which(colnames(mat_gamma_beta) == x)})
  mat_gamma_beta <- mat_gamma_beta[ , idf]
  mat2 = t(mat_gamma_beta)
  Cb = fit$C
  Cb[Cb == 1] = mat2[mat2 != 0]
  Cb
  corr_results=round(Cb[which(Cb != 0, arr.ind = TRUE)], 2)
  
  
  
  
  
  nomdf = name_block
  
  df = data.frame(id = nomdf,
                  label = nomdf,
                  shape = rep(shape, c(length(name_block))),
                  color.background = rep(color.background, 
                                         c(length(name_block))),
                  color.border = rep(color.border, 
                                     c(length(name_block)))
  )
  
  
  fusion = round(corr_results, 2)
  # "from" nodes
  from = inner[, 1]
  # "to" nodes
  to = inner[, 2]
  color.edge=rep(color.edge, c(length(inner[ , 1])))
  # data.frame which make the link between the nodes with optional labels
  #Smooth = T prevents overlapping of the labels on the edges
  if(display_numbers)
    link <- data.frame(from = from  , to = to, 
                       color = color.edge, label = fusion, smooth = TRUE) else 
                         link <- data.frame(from = from, color = color.edge, 
                                            to = to, smooth = FALSE)
  
  valeur = ifelse(display_numbers == TRUE, FALSE, TRUE)
  
  visNetwork(df, link, width = "100%")%>%
    visEdges(arrows = "to",smooth = valeur)%>%
    visInteraction(multiselect = TRUE)%>%
    visConfigure(enabled = configure)%>%
    visPhysics(enabled = valeur)
}