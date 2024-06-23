#' Estimate C matrix with bassPCA as a function of t
#'
#' Closed form estimator of the C(t) matrix using a BASS model
#'
#' @param modPCA a fitted model of class bassBasis from bassPCA function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param func.use a vector of points of the functional variable to use
#' @return A list returning the (posterior samples?) of the C matrix for each point specified in func.use
#' @details This function works by converting the linear combination of bass models to a single bass model. See C_bass for more details
#' @export
C_bassPCA <- function(modPCA, prior=NULL, mcmc.use=NULL, func.use=NULL){
  if(class(modPCA) != "bassBasis") stop("modPCA must be an object with class bassBasis")

  # Get weights matrix
  phi <- modPCA$dat$basis
  nfunc <- nrow(phi)
  if(is.null(func.use)) func.use <- 1:nfunc

  # Get model list
  mod_list <- modPCA$mod.list

  Ct <- list()
  for(t in seq_along(func.use)){
    tt <- func.use[t]
    fit_lc_curr <- lcbass2bass(mod_list, weights=phi[tt,])
    Ct[[t]] <- C_bass(fit_lc_curr, prior, mcmc.use)
  }

  return(Ct)
}

#' Estimate C matrix with bassPCA as a function of t
#'
#' Closed form estimator of the C(t) matrix using a BASS model.
#' An alternative approach see details.
#'
#' @param modPCA a fitted model of class bassBasis from bassPCA function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param func.use a vector of points of the functional variable to use
#' @return A list returning the (posterior samples?) of the C matrix for each point specified in func.use
#' @details This function works by decomposing the C of a linear combination into the pairwise Cfg matrices of the components. See C_bass for more details
#' @export
C_bassPCA_v2 <- function(modPCA, prior=NULL, mcmc.use=NULL, func.use=NULL){
  if(class(modPCA) != "bassBasis") stop("modPCA must be an object with class bassBasis")

  # Get weights matrix
  phi <- modPCA$dat$basis
  nfunc <- nrow(phi)
  if(is.null(func.use)) func.use <- 1:nfunc

  # Get model list
  mod_list <- modPCA$mod.list

  # Get Cij for all model pairs
  nbassmodels <- length(modPCA)
  Cij <- list()
  cnt <- 1
  for(i in 1:nbassmodels){
    for(j in 1:nbassmodels){
      if(i == j){
        Cij[[cnt]] <- C_bass(mod_list[[i]], prior, mcmc.use)
      }else{
        Cij[[cnt]] <- Cfg_bass(mod_list[[i]], mod_list[[j]], prior, mcmc.use)
      }
      cnt <- cnt + 1
    }
  }
  p <- nrow(Cij[[1]])

  # Assemble C list for each t
  Ct <- list()
  cnt <- 1
  for(t in seq_along(func.use)){
    tt <- func.use[t]
    Ctmp <- matrix(0, nrow=p, ncol=p)
    for(i in 1:nbassmodels){
      for(j in 1:nbassmodels){
        Ctmp <- Ctmp + phi[tt,i]*phi[tt,j]*Cij[[cnt]]
        cnt <- cnt + 1
      }
    }
    Ct[[t]] <- Ctmp
  }
  return(Ct)
}

#' Estimate Cfg matrix with bassPCA as a function of t
#'
#' Closed form estimator of the Cfg(t) matrix using a BASS model
#'
#' @param modPCA1 a fitted model of class bassBasis from bassPCA function
#' @param modPCA2 a fitted model of class bassBasis from bassPCA function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param func.use a vector of points of the functional variable to use
#' @return A list returning the (posterior samples?) of the Cfg matrix for each point specified in func.use
#' @details This function works by converting the linear combination of bass models to a single bass model. See Cfg_bass for more details
#' @export
Cfg_bassPCA <- function(modPCA1, modPCA2, prior=NULL, mcmc.use=NULL, func.use=NULL){
  if(class(modPCA) != "bassBasis") stop("modPCA must be an object with class bassBasis")

  # Get weights matrix
  phi1 <- modPCA1$dat$basis
  phi2 <- modPCA1$dat$basis
  nfunc <- nrow(phi1)
  if(is.null(func.use)) func.use <- 1:nfunc

  # Get model list
  mod_list1 <- modPCA1$mod.list
  mod_list2 <- modPCA2$mod.list
  if(is.null(mcmc.use)) mcmc.use <- 1:length(mod_list1[[1]]$s2)

  Cfgt <- list()
  for(t in seq_along(func.use)){
    tt <- func.use[t]
    fit_lc_curr1 <- lcbass2bass(mod_list1, weights=phi1[tt,])
    fit_lc_curr2 <- lcbass2bass(mod_list2, weights=phi2[tt,])
    Cfgt[[t]] <- Cfg_bass(fit_lc_curr1, fit_lc_curr2, prior, mcmc.use)
  }
  return(Cfgt)
}


#' Estimate Cfg matrix with bassPCA as a function of t
#'
#' Closed form estimator of the Cfg(t) matrix using a BASS model.
#' An alternative approach, see details.
#'
#' @param modPCA1 a fitted model of class bassBasis from bassPCA function
#' @param modPCA2 a fitted model of class bassBasis from bassPCA function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param func.use a vector of points of the functional variable to use
#' @return A list returning the (posterior samples?) of the Cfg matrix for each point specified in func.use
#' @details This function works by decomposing the Cfg of a linear combination into the pairwise Cfigj matrices of the components. See Cfg_bass for more details
#' @export
Cfg_bassPCA_v2 <- function(modPCA1, modPCA2, prior=NULL, mcmc.use=NULL, func.use=NULL){
  if(class(modPCA1) != "bassBasis") stop("modPCA1 must be an object with class bassBasis")
  if(class(modPCA2) != "bassBasis") stop("modPCA2 must be an object with class bassBasis")

  # Get weights matrix
  phi1 <- modPCA1$dat$basis
  phi2 <- modPCA2$dat$basis
  if(any(dim(phi1) != dim(phi2))){
    stop("modPCA$dat$basis should hve the same dimension for both models.")
  }
  nfunc <- nrow(phi1)
  if(is.null(func.use)) func.use <- 1:nfunc

  # Get model list
  mod_list1 <- modPCA1$mod.list
  mod_list2 <- modPCA2$mod.list

  # Get Cij for all model pairs
  nbassmodels1 <- length(modPCA1)
  nbassmodels2 <- length(modPCA2)
  Cij <- list()
  cnt <- 1
  for(i in 1:nbassmodels1){
    for(j in 1:nbassmodels2){
      Cij[[cnt]] <- Cfg_bass(mod_list1[[i]], mod_list2[[j]], prior, mcmc.use)
      cnt <- cnt + 1
    }
  }
  p <- nrow(Cij[[1]])

  # Assemble C list for each t
  Cfgt <- list()
  cnt <- 1
  for(t in seq_along(func.use)){
    tt <- func.use[t]
    Ctmp <- matrix(0, nrow=p, ncol=p)
    for(i in 1:nbassmodels1){
      for(j in 1:nbassmodels2){
        Ctmp <- Ctmp + phi1[tt,i]*phi2[tt,j]*Cij[[cnt]]
        cnt <- cnt + 1
      }
    }
    Cfgt[[t]] <- Ctmp
  }
  return(Cfgt)
}
