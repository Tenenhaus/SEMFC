#' Outer Static Visualization Network RGCCAC
#'
#'path_diagram_outer allows to display an outer graphic of an RGCCAC object 
#'
#' @param model fitted rgccac object
#' @param isplay_numbers boolean, default to TRUE, if display_numbers = TRUE, 
#' display the cor_eta_y between latents and manifests nodes 
#' and the betas/gammas between the latents nodes 
#' @param graph_options string,default to "dot". The layout among "dot and "fdp". 
#' The rankdir sets direction of graph layout. Default to "LR" 
#' Concentrate allows to merges multiedges into a single edge.
#' This feature is not  available outside of "dot" layout
#'  and display_number = TRUE. Default to TRUE
#' @param dge_options string. The color of the edges and shape of the arrow
#' Default to list(color = "black", arrowhead = "vee")
#' @param node_options_LV string. The shape and the color of the latents nodes
#' Default to list(shape = "oval", color = "black")
#' @param num string or integer. Display the outer graph of the num's block or
#' all the blocks if num  = "all"
#'
#' @return
#' Outer static graph for rgcca object 
#' @export
#'
#' @examples
#' path_diagram_outer (fit.rgccac, display_numbers = TRUE,
#'graph_options = list(overlap = F, layout = "dot", 
#'                     rankdir = "UD", concentrate = T),
#'edge_options = list(color = "black", arrowhead = "vee"),
#'node_options_LV = list(shape = "oval",  color = "black"),
#'node_options_MV = list(shape = "box",  color = "black"), num = "all")
path_diagram_outer = function(model, display_numbers = TRUE,
                        graph_options = list(overlap = F, layout = "dot", 
                                             rankdir = "UD", concentrate = T),
                        edge_options = list(color = "black", arrowhead = "vee"),
                        node_options_LV = list(shape = "oval",  color = "black"),
                        node_options_MV = list(shape = "box",  color = "black")
                        ,num = num)
{
  
  plotCall <- BuildCall_outer(model = model, display_numbers = display_numbers,
                              graph_options = graph_options,
                              edge_options = edge_options,
                              node_options_LV = node_options_LV,
                              node_options_MV = node_options_MV,num = num
  )
  
  DiagrammeR:::grViz(plotCall)
}