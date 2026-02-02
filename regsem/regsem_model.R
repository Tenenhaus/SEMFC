


regssem_model <- function(len_block){

  pen_vars <- paste0("X1", 2:len_block, collapse = " + ")
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


  sem.model.lslx <-  paste0('
  # latent variable definitions
  eta1 =~ X11
  eta2 =~ X21+X22+X23
  eta3 =~ X31+X32+X33
  eta4 =~ X41+X42+X43
  eta5 =~ X51+X52+X53
  eta6 =~ X61+ X62+X63

  # penalisation
  ',pen_line,'

  # Regressions
  eta5 ~ eta1 + eta2 + eta6
  eta6 ~ eta3 + eta4 + eta5

  #variance

  eta1 ~~ 1*eta1
  eta2 ~~ 1*eta2
  eta3 ~~ 1*eta3
  eta4 ~~ 1*eta4
  eta5 ~~ 1*eta5
  eta6 ~~ 1*eta6




  # residual covariances
  eta5 ~~ eta6


  ')


  sem.model.lsem <-paste0('
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
  eta5 ~~ eta6

  ')


  return(list(sem.model=sem.model,
              sem.model.lslx=sem.model.lslx,
              sem.model.lsem=sem.model.lsem))
  }