#' Dynamic Outer Network Visualization for fitted RGCCAC object 
#'
#'Gives a dynamic outer graph where you can change the colors and the arrangement 
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
#' @param num string or integer. Display the outer graph of the num block or
#' all the blocks if num  = "all"
#'
#' @return
#' Outer dynamic graph 
#' @export
#'
#' @examples
#' outerplot (fit.rgccac, color.background = c("white", "white"),
#'color.border = c("black", "black"), 
#'shape = c("oval", "box"), 
#'color.edge = c("black", "black"),
#'display_numbers = TRUE, configure = TRUE, num = "all)
outerplot = function(fit, color.background = c("white", "white"),
                     color.border = c("black", "black"), 
                     shape = c("oval", "box"), 
                     color.edge = c("black", "black"),
                     display_numbers = TRUE, configure = FALSE, num = "all"){
  
  tau_rep = rep(1,sum(sapply(fit$cor_y_eta, length)))
  
  name_block = names(fit$d) 
  inner = matrix(name_block[which(fit$C !=0, arr.ind = TRUE)], ncol = 2)
  
  
  fit_cor_y_eta = as.vector(as.vector(Reduce("c", fit$cor_y_eta)))
  
  fusion = round(fit_cor_y_eta, 2)
  
  
  ManifVar = as.vector(unlist(sapply(fit$cor_y_eta,rownames)))
  nomdf = c(name_block, ManifVar)
  
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
  
  from = MV
  
  # "to" nodes
  
  to = LV
 
  color.edge =  rep(color.edge,1, length(MV))
                   
  
  # data.frame which make the link between the nodes with optional labels
  #Smooth = T prevents overlapping of the labels on the edges
  if(display_numbers)
    link <- data.frame(from = from  , to = to, 
                       color = color.edge, label = fusion, 
                       smooth = FALSE) else 
    link <- data.frame(from = from, color = color.edge, 
                       to = to, smooth = FALSE)
  
  q=list()
  m=list()
  n=list()
  
  for(i in 1:length(fit$d)){
    q[[i]]=link[which(link[,1]==names(fit$d[i])),]
    m[[i]]=df[which(df[,2]%in%q[[i]][,2]),]
    n[[i]]=rbind(df[i,],m[[i]])
  }
  if(num == "all"){
    visNetwork(df, link, width = "100%")%>%
      visEdges(arrows="to", smooth = FALSE)%>%
      visInteraction(multiselect = TRUE)%>%
      visConfigure(enabled = configure)%>%
      visPhysics(enabled = FALSE)
  }else{
    visNetwork(n[[num]], q[[num]], width = "100%")%>%
      visEdges(arrows="to",smooth = FALSE)%>%
      visInteraction(multiselect = TRUE)%>%
      visConfigure(enabled = configure)%>%
      visPhysics(enabled = FALSE)
  }
}
