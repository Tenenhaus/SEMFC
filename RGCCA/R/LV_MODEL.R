#' Title
#' Function that returns betas and gammas from a C matrix
#' @param C a matrix filled with 0 and 1 with the same number of row and column 
#' without any 1 on the diag
#' @param R a Rho matrix
#'
#' @return list of beta matrix and gamma matrix
#' @export
#'
#' @examples
#' C= matrix(c(0, 0, 0, 0, 0,
#'             1, 0, 0, 0, 0,
#'             1, 1, 0, 0, 0,
#'             1, 1, 1, 0, 0,
#'             0, 0, 0, 1, 0),5, 5)
#'             
#'R = matrix(c(1, 0.8893068, 0.6090810, 0.8331237, 0.6754145,
#'             0.8893068, 1, 0.6785475, 0.9545197, 0.6750227, 
#'             0.6090810, 0.6785475, 1, 0.7529390, 0.6878589,
#'             0.8331237, 0.9545197, 0.7529390, 1, 0.8617895,
#'             0.6754145, 0.6750227, 0.6878589, 0.8617895, 1),
#'             5, 5)
#'             
#' LV_MODEL(R,C)
#' 
LV_MODEL=function(R, C){
  g = igraph::graph_from_adjacency_matrix(C)
  which_exo_endo = ind_exo_endo(C)
  rownames(C) = colnames(C) = rownames(R)
  xx = list()
  xy = list()
  beta_gamma = list()
  beta = list()
  gamma = list()
  tC = t(C)
  H = which_exo_endo$ind_exo
  J = which_exo_endo$ind_endo
  
  mat_beta = as.matrix(tC[J, J, drop = F])
  mat_gamma = as.matrix(tC[J, H, drop = F])
  
  for (i in 1:length(which_exo_endo$Ji)){
    Ji = which_exo_endo$Ji[[i]]
    Hi = which_exo_endo$Hi[[i]]
    
    #For non-recursive structural model
    if(!igraph::is_dag(g) & any(Ji != 0)){
      #TSLS
      xx[[i]] = ifelse(length(Ji) > 1 & length(Hi) == 1,
                  ##YES##
                  list(
                  rbind(cbind(R[Ji, H]%*%solve(R[H, H])%*%R[H, Ji], R[Ji, Hi]),
                        c(R[Hi, Ji], R[Hi, Hi]))
                  ),
        
                  ##NO##
                  ifelse(length(Hi)>1 & length(Ji) == 1, 
                  #YES
                  list(
                  rbind(c(R[Ji, H]%*%solve(R[H, H])%*%R[H, Ji], R[Ji, Hi]),
                        cbind(R[Hi, Ji], R[Hi, Hi]))),
                  #NO
                  list(
                  rbind(cbind(R[Ji, H]%*%solve(R[H, H])%*%R[H, Ji], R[Ji, Hi]),
                        cbind(R[Hi, Ji], R[Hi, Hi]))))
                  )
      
      xy[[i]] = c(R[Ji, H]%*%solve(R[H, H])%*%R[H, J[[i]]], R[Hi, J[[i]]])
      
      beta_gamma[[i]] = solve(xx[[i]][[1]])%*%xy[[i]]
      
    }else{
      #OLS
      xx[[i]] = ifelse(length(Ji) > 1 & length(Hi) == 1,
                       ##YES##
                       list(
                       rbind(cbind(R[Ji, Ji], R[Ji, Hi]),
                             c(R[Hi, Ji], R[Hi, Hi]))),
          
                       ##NO##
                       ifelse(length(Hi)>1 & length(Ji) == 1,
                       #YES#
                       list(
                       rbind(c(R[Ji, Ji], R[Ji, Hi]),
                             cbind(R[Hi, Ji], R[Hi, Hi]))),
                       #NO#
                       list(
                       rbind(cbind(R[Ji, Ji], R[Ji, Hi]),
                             cbind(R[Hi, Ji], R[Hi, Hi]))))
                       )
      
      xy[[i]] = c(R[Ji, J[i]], R[Hi, J[i]])
      beta_gamma[[i]] = solve(xx[[i]][[1]])%*%xy[[i]]
      
    }

    # beta[[i]] = ifelse(Ji == 0, 
    #                    0,
    #                    beta_gamma[[i]][1:length(Ji)]
    #                    )
    # 
    # gamma[[i]] = ifelse(Hi == 0, 
    #                     0,
    #                     ifelse(Ji == 0,
    #                     beta_gamma[[i]], 
    #                     list(beta_gamma[[i]][seq(length(beta_gamma[[i]]))[-seq(length(Ji))]]))
    #                     )
    
    ifelse(Ji == 0, 
           yes = {beta[[i]] = 0},
           no = {beta[[i]] = list(beta_gamma[[i]][1:length(Ji)])} 
    )
    
    ifelse(Hi == 0, 
      yes = {gamma[[i]] = 0},
      no = {ifelse(Ji == 0,
            yes = {gamma[[i]] = list(beta_gamma[[i]])},
            no = {gamma[[i]] = list(beta_gamma[[i]][seq(length(beta_gamma[[i]]))[-seq(length(Ji))]])})
           })
    
     #beta_gamma[[i]][!beta_gamma[[i]]%in%beta[[i]]]
    
    mat_beta[i, ] = replace(mat_beta[i, ], mat_beta[i, ] == 1, beta[[i]][[1]])
    mat_gamma[i, ] = replace(mat_gamma[i, ], mat_gamma[i , ] == 1, gamma[[i]][[1]])
  }
  
  return(list( mat_beta = mat_beta, mat_gamma = mat_gamma))
}
