
#' Estimate K matrix with bassPCA
#'
#' Closed form estimator of the K matrix using a BASS model
#'
#' @param modPCA a fitted model of class bassBasis from bassPCA function
#' @param type 1 or 2. Use C_bassPCA or C_bassPCA_v2?
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param prior_func a vector of weights to use when summing over functional variable. Should be same length as func.use.
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param func.use a vector of points of the functional variable to use
#' @return A list returning the (posterior samples?) of the C matrix for each point specified in func.use
#' @details This function works by converting the linear combination of bass models to a single bass model. See C_bass for more details
#' @export
K_bassPCA <- function(modPCA, type=1, prior=NULL, prior_func=NULL, mcmc.use=NULL, func.use=NULL){
  if(class(modPCA) != "bassBasis") stop("modPCA must be an object with class bassBasis")

  # Get weights matrix
  phi <- modPCA$dat$basis
  nfunc <- nrow(phi)
  if(is.null(func.use)) func.use <- 1:nfunc
  if(is.null(prior_func)) prior_func <- rep(1/nfunc, nfunc)

  # Call C_bassPCA
  C_bassPCA_type <- ifelse(type == 1, get("C_bassPCA"), get("C_bassPCA_v2"))
  Clist <- C_bassPCA_type(modPCA, prior, mcmc.use, func.use)

  # Get weighted sum
  Ktmp <- Clist[[1]]
  if(is.matrix(Ktmp)){
    K <- Ktmp*prior_func[1]
    for(t in seq_along(func.use)[-1]){
      K <- K + Clist[[t]]*prior_func[t]
    }
  }else{
    nmcmc <- length(Ktmp)
    K <- list()
    Ktmp <- Ktmp*prior_func[1]
    for(mm in 1:nmcmc){
      for(t in seq_along(func.use)[-1]){
        Ctmp <- Clist[[t]]
        K[[mm]] <- K[[mm]] + Ctmp[[mm]]*prior_func[t]
      }
    }
  }
  return(K)
}

Kfg_bass <- function(modPCA1, modPCA2, type=1, prior=NULL, prior_func=NULL, symmetric=FALSE, mcmc.use=NULL, func.use=NULL){
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

  # Call C_bassPCA
  Cfg_bassPCA_type <- ifelse(type == 1, get("Cfg_bassPCA"), get("Cfg_bassPCA_v2"))
  Clist <- Cfg_bassPCA_type(modPCA1, modPCA2, prior, mcmc.use, func.use)

  # Get weighted sum
  Ktmp <- Clist[[1]]
  if(is.matrix(Ktmp)){
    K <- Ktmp*prior_func[1]
    for(t in seq_along(func.use)[-1]){
      K <- K + Clist[[t]]*prior_func[t]
    }
  }else{
    nmcmc <- length(Ktmp)
    K <- list()
    Ktmp <- Ktmp*prior_func[1]
    for(mm in 1:nmcmc){
      for(t in seq_along(func.use)[-1]){
        Ctmp <- Clist[[t]]
        K[[mm]] <- K[[mm]] + Ctmp[[mm]]*prior_func[t]
      }
    }
  }

  # Symmetrize if desired
  if(symmetric){
    if(is.matrix(K)){
      K <- (K + t(K))/2
    }else{
      for(mm in 1:nmcmc){
        K[[mm]] <- (K[[mm]] + t(K[[mm]]))/2
      }
    }
  }
  return(K)
}
