#' Convert a linear combination of BASS models to a single BASS model
#'
#' A linear combination of BASS models is also a BASS model. This function takes a list of BASS models (all with the same data matrix xx.des)
#' and returns the resulting linear combination as a new BASS model.
#' One useful application of this function is to convert bassPCA to bass for a fixed time point.
#' Does not currently work for bass models with functional or categorical inputs.
#'
#' @param mod_list A list of bass models.
#' @param weights An optional vector of weights.
#' @param yy The data vector. Optional, but useful for some bass object methods.
#' @param mcmc.use set of indices telling which mcmc draws to use.
#' @export
#' @examples
#' a <- 1
#' @export
lcbass2bass <- function(mod_list, weights=rep(1, length(mod_list)), yy=NULL, mcmc.use=NULL){
  if(!is.list(mod_list)){
    stop("mod_list must be a list")
  }
  if(!all(unlist(lapply(lapply(mod_list, class), function(cl) "bass" %in% cl)))){
    stop("Every component of mod_list should be a bass model")
  }
  if(is.null(mcmc.use)){
    mcmc.use <- seq_along(mod_list[[1]]$s2)
  }
  # Quick pass over models (nothing expensive)
  nbassmodels <- length(mod_list)
  nmcmc <- length(mcmc.use)

  # These fields are the same for each model (assumption)
  mod_curr <- mod_list[[1]]
  nobs <- length(mod_curr$y)

  p <- mod_curr$p
  pdes <- mod_curr$pdes
  xx.des <- mod_curr$xx.des
  range.des <- mod_curr$range.des
  des <- mod_curr$des
  type <- mod_curr$type

  # These fields may be different for each bass model
  # or they might change over the course of the code
  nbasis   <- matrix(mod_curr$nbasis[mcmc.use], ncol=1)
  n.models <- length(unique(mod_curr$model.lookup[mcmc.use]))
  z <- matrix(mod_curr$y, ncol=1)
  yhat <- weights[1]*mod_curr$yhat[mcmc.use,]
  yhat.mean <- weights[1]*mod_curr$yhat.mean
  lam <- matrix(mod_curr$lam[mcmc.use], ncol=1)
  s2 <- weights[1]^2*mod_curr$s2[mcmc.use]
  beta.prec <- matrix(mod_curr$beta.prec[mcmc.use], ncol=1)
  degree <- mod_curr$degree
  maxInt.des <- mod_curr$maxInt.des
  lookup <- matrix(mod_curr$model.lookup[mcmc.use], nrow=length(mcmc.use), ncol=1)

  if(nbassmodels > 1){
    for(i in 2:nbassmodels){
      mod_curr <- mod_list[[i]]
      nbasis <- cbind(nbasis, mod_curr$nbasis[mcmc.use])
      n.models <- c(n.models, length(unique(mod_curr$model.lookup[mcmc.use])))
      z <- cbind(z, mod_curr$y)
      yhat <- weights[i]*mod_curr$yhat[mcmc.use,]
      yhat.mean <- weights[i]*mod_curr$yhat.mean
      lam <- cbind(lam, mod_curr$lam)
      s2 <- weights[i]^2*mod_curr$s2
      beta.prec <- cbind(beta.prec, mod_curr$beta.prec[mcmc.use])
      degree <- c(degree, mod_curr$degree)
      maxInt.des <- c(maxInt.des, mod_curr$maxInt.des)
      lookup <- cbind(lookup, mod_curr$model.lookup[mcmc.use])
    }
  }

  # Extract information about the final model
  maxBasis <- max(apply(nbasis, 1, sum))
  maxInt <- max(maxInt.des)
  maxModels <- sum(n.models)

  n.int <- matrix(NA, nrow=maxModels, ncol=maxBasis)
  knots <- signs <- vars <- array(NA, dim=c(maxModels, maxBasis, maxInt))
  beta  <- matrix(NA, nrow=nmcmc, ncol=maxBasis + 1) # Plus 1 for intercept term

  nbasisFull  <- rep(NA, nmcmc)
  lookupFull  <- rep(NA, nmcmc)
  nmodelsFull <- 0

  # Start loop over mcmc iterations
  for(mm in seq_along(mcmc.use)){
    if(mm == 1){
      new_model_flag <- TRUE
    }else{
      if(any(lookup[mm,] != lookup[mm-1,])){
        new_model_flag <- TRUE
      }else{
        new_model_flag <- FALSE
      }
    }

    if(new_model_flag){
      # Combine models sequentially
      for(ii in 1:nbassmodels){
        mod_curr <- mod_list[[ii]]
        lookup_curr <- lookup[mm,ii] # first index of the arrays
        nb_mm_ii <- nbasis[mm,ii]
        #if(nb_mm_ii == 0) browser()
        tmp <- mod_curr$beta[mm,1:(1+nb_mm_ii)] * weights[ii]
        beta0_curr <- tmp[1]
        beta_curr <- tmp[-1]

        n.int_curr <- mod_curr$n.int.des[lookup_curr,seq_along(beta_curr),drop=FALSE]
        knots_curr <- mod_curr$knotInd.des[lookup_curr,seq_along(beta_curr),]
        signs_curr <- mod_curr$signs.des[lookup_curr,seq_along(beta_curr),]
        vars_curr <- mod_curr$vars.des[lookup_curr,seq_along(beta_curr),]

        if(ii == 1){
          n.int_cand <- n.int_curr
          knots_cand <- knots_curr
          signs_cand <- signs_curr
          vars_cand  <- vars_curr
          beta0_cand <- beta0_curr
          beta_cand  <- beta_curr
        }else{
          n.int_cand <- c(n.int_cand, n.int_curr)
          knots_cand <- rbind(knots_cand, knots_curr)
          signs_cand <- rbind(signs_cand, signs_curr)
          vars_cand  <- rbind(vars_cand, vars_curr)
          beta0_cand <- beta0_cand + beta0_curr
          beta_cand  <- c(beta_cand, beta_curr)
        }
      }

      ## Simplify model, in the case of any duplicates.
      model_cand <- cbind(knots_cand, signs_cand, vars_cand)
      duplicates <- find_duplicate_groups(model_cand)
      for(jj in seq_along(duplicates)){
        dup_jj <- duplicates[[jj]]

        # Add the coefficients together
        beta_cand[dup_jj[1]] <- sum(beta_cand[dup_jj])

        # Remove the excess rows
        to_remove <- dup_jj[-1]
        n.int_cand <- n.int_cand[-to_remove]
        knots_cand <- knots_cand[-to_remove,]
        signs_cand <- signs_cand[-to_remove,]
        vars_cand <- vars_cand[-to_remove,]
        beta_cand <- beta_cand[-to_remove]
      }

      # Check to see if this new model is already in the full model
      modelFull <- abind_custom(knots, signs, vars, along=3)
      duplicate_model <- any(apply(modelFull, 1,
                                   function(slice) identical(slice, model_cand)))
      # We can update the coefficients regardless of what happens
      nbasis_curr   <- length(n.int_cand)
      beta[mm,1:(nbasis_curr+1)] <- c(beta0_cand, beta_cand)


      if(duplicate_model){
        # If duplicate, just set lookup and nbasis to previous values
        lookupFull[mm] <- lookupFull[mm-1]
        nbasisFull[mm] <- nbasisFull[mm-1]
      }else{
        # If not duplicate, add everything to the full model
        nbasis_curr   <- length(n.int_cand)
        nmodelsFull   <- nmodelsFull + 1

        nbasisFull[mm] <- nbasis_curr
        lookupFull[mm] <- nmodelsFull

        # Store parameters
        beta[mm,1:(nbasis_curr+1)] <- c(beta0_cand, beta_cand)
        n.int[nmodelsFull,1:nbasis_curr] <- n.int_cand
        knots[nmodelsFull,1:nbasis_curr,] <- knots_cand
        signs[nmodelsFull,1:nbasis_curr,] <- signs_cand
        vars[nmodelsFull,1:nbasis_curr,] <- vars_cand
      }
    }else{ #If model is the same as before
      # Set things to be the same as the previous iteration
      lookupFull[mm] <- lookupFull[mm-1]
      nbasisFull[mm] <- nbasisFull[mm-1]

      # Get new coefficient vector
      for(ii in 1:nbassmodels){
        mod_curr <- mod_list[[ii]]
        tmp <- mod_curr$beta[mm,1:(1+nbasis[mm,ii])] * weights[ii]
        beta0_curr <- tmp[1]
        beta_curr <- tmp[-1]
        if(ii == 1){
          beta0_cand <- beta0_curr
          beta_cand  <- beta_curr
        }else{
          beta0_cand <- beta0_cand + beta0_curr
          beta_cand  <- c(beta_cand, beta_curr)
        }
      }

      # Remove duplicates and combine betas
      # NOTE: duplicates will be the same as it was last iteration,
      #       so no need to recompute it
      for(jj in seq_along(duplicates)){
        dup_jj <- duplicates[[jj]]
        # Add the coefficients together
        beta_cand[dup_jj[1]] <- sum(beta_cand[dup_jj])
        # Remove the excess rows
        to_remove <- dup_jj[-1]
        beta_cand <- beta_cand[-to_remove]
      }
      beta[mm,1:(nbasis_curr+1)] <- c(beta0_cand, beta_cand)
    }
  } # end loop over mcmc

  # Trim excess NA's off of data frames
  maxBasisFull <- max(nbasisFull)
  n.int <- n.int[1:nmodelsFull, 1:maxBasisFull]
  knots <- knots[1:nmodelsFull, 1:maxBasisFull,]
  signs <- signs[1:nmodelsFull, 1:maxBasisFull,]
  vars  <- vars[1:nmodelsFull, 1:maxBasisFull,]
  beta <- beta[,1:(maxBasisFull+1)]

  # Create object to return
  out <- mod_list[[1]]
  out$yhat.mean <- yhat.mean
  out$yhat <- yhat
  out$beta <- beta
  out$s2 <- s2
  out$lam <- lam
  out$nbasis <- nbasisFull
  out$degree.all <- degree
  out$degree <- degree[1]
  out$nmcmc <- nmcmc
  out$nburn <- 0
  out$thin <- 1
  out$p <- p
  out$pdes <- p
  out$beta.prec <- beta.prec
  if(is.null(yy)){
    out$y <- z[,1]
  }else{
    out$y <- yy
  }
  out$z <- z
  out$xx.des <- xx.des
  out$n.models <- nmodelsFull
  out$model.lookup <- lookupFull
  out$des <- des
  out$type <- type
  out$knotInd.des <- knots
  out$signs.des <- signs
  out$vars.des <- vars
  out$n.int.des <- n.int
  out$maxInt.des <- maxInt
  out$range.des <- range.des
  out$func <- FALSE
  out$cx <- rep("numeric", p)
  out$pfunc <- 0
  class(out) <- "bass"
  return(out)
}







