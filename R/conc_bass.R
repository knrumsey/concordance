#' Concordance analysis using bass models
#'
#' Closed form estimator of the Cfg matrix using a BASS model
#'
#' @param mod1 a fitted BASS model for first function
#' @param mod2 a fitted BASS model for second function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @param type Only used if class(mod1) == class(mod2) == "bassBasis". See \code{C_bassPCA}, \code{C_bassPCA_v2}, \code{Cfg_bassPCA}, and \code{Cfg_bassPCA_v2} for details
#' @param prior_func Only used if class(mod1) == class(mod2) == "bassBasis". See \code{C_bassPCA}, \code{C_bassPCA_v2}, \code{Cfg_bassPCA}, and \code{Cfg_bassPCA_v2} for details
#' @param func.use Only used if class(mod1) == class(mod2) == "bassBasis". See \code{C_bassPCA}, \code{C_bassPCA_v2}, \code{Cfg_bassPCA}, and \code{Cfg_bassPCA_v2} for details
#' @return A list with matrices Cf, Cg, Cfg, Vfg and the concordance.
#' @details
#' When models are class bass, each field of the returned object is a list for each mcmc iteration (or a vector for concordance). If mcmc.use = NULL or length(mcmc.use) = 1, then each field is just a matrix (or a scalar for concordance).
#'
#' When models are class bassBasis (from bassPCA function), each field will be a list for each time point in func.use (func.use = NULL uses all time points in the training data by default). Each component of the list has the same structure as described above for the class == "bass" case.
#'
#' @export
conc_bass <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, type=1, prior_func=NULL, func.use=NULL){
  if(class(mod1) != class(mod2)) stop("mod1 and mod2 should have the same class")
  if(class(mod1) == "bassBasis"){
    conc_bassPCA(mod1, mod2, prior, mcmc.use, type, prior_func, func.use)
  }

  out <- list()
  out$Cf <- C_bass(mod1, prior, mcmc.use)
  out$Cg <- C_bass(mod2, prior, mcmc.use)
  out$Cfg <- Cfg_bass(mod1, mod2, prior, mcmc.use)

  if(is.matrix(out$Cf)){
    out$Vfg <- (out$Cfg + t(out$Cfg))/2
    out$conc <- tr(out$Cfg)/sqrt(tr(out$Cf)*tr(out$Cg))
  }else{
    nmcmc <- length(out$Cf)
    Vfg_list <- list()
    conc <- rep(NA, nmcmc)
    for(mm in 1:nmcmc){
      Vfg_list[[mm]] <- (out$Cfg[[mm]] + t(out$Cfg[[mm]]))/2
      conc[mm] <- tr(out$Cfg[[mm]])/sqrt(tr(out$Cf[[mm]])*tr(out$Cg[[mm]]))
    }
    out$Vfg <- Vfg_list
    out$conc <- conc
  }

  return(out)
}

conc_bassPCA <- function(modPCA1, modPCA2, prior = NULL, mcmc.use=NULL, type=1, prior_func=NULL, func.use=NULL){
  out <- list()
  if(type == 1){
    out$Cf <- C_bassPCA(modPCA1, prior, mcmc.use, func.use)
    out$Cg <- C_bassPCA(modPCA2, prior, mcmc.use, func.use)
    out$Cfg <- Cfg_bassPCA(modPCA1, modPCA2, prior, mcmc.use, func.use)
  }else{
    out$Cf <- C_bassPCA_v2(modPCA1, prior, mcmc.use, func.use)
    out$Cg <- C_bassPCA_v2(modPCA2, prior, mcmc.use, func.use)
    out$Cfg <- Cfg_bassPCA_v2(modPCA1, modPCA2, prior, mcmc.use, func.use)
  }

  nfunc <- length(func.use)
  if(is.matrix(out$Cfg[[1]])){
    Vfg <- list()
    conc <- rep(NA, nfunc)
    for(t in 1:nfunc){
      Vfg[[t]] <- (out$Cfg[[t]] + t(out$Cfg[[t]]))/2
      conc[t]  <- tr(Vfg[[t]])/sqrt(tr(out$Cf[[t]])*tr(out$Cg[[t]]))
    }
    conc_total2 <- sum(conc*prior_func)
  }else{
    nmcmc <- length(out$Cfg[[1]])
    Vfg <- list()
    conc <- matrix(NA, nrow=nfunc, ncol=nmcmc)
    for(t in 1:nfunc){
      Cf_t <- out$Cf[[t]]
      Cf_g <- out$Cg[[t]]
      Cfg_t <- out$Cfg[[t]]
      conc_t <- rep(NA, nmcmc)
      Vfg_t <- list()
      for(mm in 1:nmcmc){
        Vfg_t[[mm]] <- (Cfg_t[[mm]] + t(Cfg_t[[mm]]))/2
        conc_t[mm] <- tr(Cfg_t[[mm]])/sqrt(tr(Cf_t[[mm]])*tr(Cg_t[[mm]]))
      }
      Vfg[[t]] <- Vfg_t
      conc[t,] <- conc_t
    }
    conc_total2 <- apply(conc, 2, function(cc) sum(cc*func_prior))
  }
  out$Vfg <- Vfg
  out$conc <- conc
  out$conc_total2 <- conc_total2

  # Get K matrices
  out$Kf  <- K_bassPCA(modPCA1, type, prior, prior_func, mcmc.use, func.use)
  out$Kg  <- K_bassPCA(modPCA2, type, prior, prior_func, mcmc.use, func.use)
  out$Kfg <- K_bassPCA(modPCA1, modPCA2, type, prior, prior_func, mcmc.use, func.use)

  # Get conc_total (using K's)
  if(is.matrix(out$Kfg)){
    out$conc_total <- tr(out$Kfg)/sqrt(tr(out$Kf)*tr(out$Kg))
  }else{
    conc_total_vec <- rep(NA, nmcmc)
    for(mm in 1:nmcmc){
      conc_total_vec[mm] <- tr(out$Kfg[[mm]])/sqrt(tr(out$Kf[[mm]])*tr(out$Kg[[mm]]))
    }
    out$conc_total <- conc_total_vec
  }
  return(out)
}
