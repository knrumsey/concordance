#' Estimate the Expected Gradient with BASS
#'
#' Closed form estimator of Z = E(\nabla f)
#'
#' @param mod a fitted BASS model
#' @param prior NULL (default) (0,1])prior for each variable. See details for required structure of prior
#' @param mcmc.use set of indices telling which mcmc draws to use
#' @param scale01 logical (ignored in current version)
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details prior should be a list of length p (one object for each variable). Each element of prior should be a named list with fields.
#' See also the concordance::build_prior() function.
#'
#'    - dist - ("uniform", "normal", "beta", "gamma").
#'    - trunc - truncation bounds (a, b). These should be c(0, 1) for "beta" and c(0, Inf) for "gamma".
#'    - mean - vector of means (mixture of normals only)
#'    - sd - vector of sds (mixture of normals only)
#'    - shape1, shape2 - shape parameters for beta distribution
#'    - shape, scale - parameters for gamma distribution
#'    - weights - vector of mixture weights (currently only compatible with dist="normal")
#' @export
Z_bass <- function(mod, prior = NULL, mcmc.use=NULL, scale01=FALSE){
  if(is.null(mcmc.use)){
    mcmc.use <- length(mod$nbasis)
  }
  gbass_flag <- "gbass" %in% class(mod)
  if(gbass_flag){
    mod$knotInd.des <- mod$knots.des
  }
  # Parse the prior information
  if(is.null(prior)){
    prior <- list()
    for(i in 1:mod$pdes){
      prior[[i]] <- list(dist="uniform")
    }
  }
  #HANDLE PRIORS
  for(i in 1:(mod$pdes)){
    # 1. Scale prior truncation bounds to BASS scale
    if(is.null(prior[[i]]$trunc)){
      prior[[i]]$trunc <- mod$range.des[,i]
    }
    prior[[i]]$trunc <- scale_range(prior[[i]]$trunc, mod$range.des[,i])

    # 2. Handle each distribution type separately
    distribution = prior[[i]]$dist
    if(distribution == "normal"){
      if(is.null(prior[[i]]$weights)){
        num_mix <- length(prior[[i]]$mean)
        prior[[i]]$weights <- rep(1/num_mix, num_mix)
      }
      prior[[i]]$mean <- scale_range(prior[[i]]$mean, mod$range.des[,i])
      prior[[i]]$sd <- prior[[i]]$sd/(mod$range.des[2,i] - mod$range.des[1,i])
      prior[[i]]$z <- pnorm((prior[[i]]$trunc[2]-prior[[i]]$mean)/prior[[i]]$sd) - pnorm((prior[[i]]$trunc[1]-prior[[i]]$mean)/prior[[i]]$sd)
      cc <- sum(prior[[i]]$weights*prior[[i]]$z)
      prior[[i]]$weights <- prior[[i]]$weights/cc
      # DF: prior[[i]]$z
      # change weights with truncation
      # divide by cc instead to keep the same prior shape
      # does the truncation change the distribution shape in the non-truncated regions??

      #Check for extrapolation
      qq <- qnorm(c(0.001, 0.999), prior[[i]]$mean, prior[[i]]$sd)
      if((qq[1] < 0 - 0.1) | (qq[2] > 1 + 0.1)){
        warning('You are asking the emulator to extrapolate. This is not reccomended.')
      }
    }
    if(distribution == "uniform"){
      prior[[i]]$weights = 1
      #Check for extrapolation
      if(((prior[[i]]$trunc[1] < 0 - 0.1) | (prior[[i]]$trunc[2] > 1 + 0.1))){
        warning('You are asking the emulator to extrapolate. This is not reccomended.')
      }
    }
    if(distribution == "beta"){
      prior[[i]]$weights = 1
      if(abs(min(mod$xx.des[,i]) - 0) > 0.01 | abs(max(mod$xx.des[,i]) - 1) > 0.01){
        warning("For dist=\"beta\", data should be scaled to (0, 1) range BEFORE fitting bass model")
      }
    }
    if(distribution == "gamma"){
      prior[[i]]$weights = 1
      if(abs(min(mod$xx.des[,i]) - 0) > 0.01){
        warning("For dist=\"gamma\", data should be shifted to have a minimum of 0 BEFORE fitting bass model")
      }
      prior[[i]]$rate = prior[[i]]$rate * prior[[i]]$trunc[2]
    }
  }

  # Make Z vector
  Zf_post <- list()
  if(!gbass_flag){
    Xt <- mod$xx.des
  }
  # Get transformation matrix
  A_tform <- diag(1/apply(mod$range.des, 2, diff))
  for(r in 1:length(mcmc.use)){
    #Compute only the stuff we will need for every iteration
    rr <- mcmc.use[r]
    mod_number_new <- mod$model.lookup[rr]
    coeff      <- mod$beta[rr,]
    coeff      <- matrix(coeff[!is.na(coeff)][-1], nrow=1)
    M_new <- length(coeff)

    compute_flag <- FALSE
    if(r == 1){
      compute_flag <- TRUE
    }else{
      if(mod_number != mod_number_new){
        compute_flag <- TRUE
      }
    }

    if(compute_flag){
      mod_number <- mod_number_new
      M          <- M_new
      signs      <- mod$signs.des[mod_number, 1:M, ]
      indic      <- mod$vars.des[mod_number, 1:M, ]
      knots      <- mod$knotInd.des[mod_number, 1:M, ]
      if(M==1){
        signs <- matrix(signs, nrow=1, ncol=length(signs))
        indic <- matrix(indic, nrow=1, ncol=length(indic))
        knots <- matrix(knots, nrow=1, ncol=length(knots))
      }
      #if(gbass_flag){
      #  knots    <- mod$knots.des[mod_number, 1:M, ]
      #}else{
      #  knots    <- mod$knotInd.des[mod_number, 1:M, ]
      #}
      # Initalize arrays
      A <- B <- I4 <- I5 <- array(NA, dim=c(mod$pdes, M))
      for(i in 1:mod$pdes){
        prior_i <- prior[[i]]
        v <- apply(indic, 1, function(zz) match(i, zz))
        u <- !is.na(v)
        s <- apply(cbind(signs, v), 1, function(zz) zz[zz[mod$maxInt.des + 1]])
        if(gbass_flag){
          t <- apply(cbind(knots, v), 1, function(zz) zz[zz[mod$maxInt.des + 1]])
        }else{
          t <- Xt[apply(cbind(knots, v), 1, function(zz) zz[zz[mod$maxInt.des + 1]]), i]
        }

        # If this comes from BASS (rather than GBASS with gm2bm)
        # then we need to account for Devin's g-scaling-factors
        #if(!gbass_flag){
        # NOTE: Need to comment out, i don't remember why
        d <- 1/((s + 1)/2 - s*t)
        s <- s*d
        #}

        #Handle NA cases
        s[is.na(s)] <- 1
        t[is.na(t)] <- -Inf

        a <- b <- i4 <- i5 <- matrix(0, M)
        #browser()
        for(m in 1:M){
          um <- u[m]
          ssm <- s[m]
          sm <- sign(ssm)
          tm <- t[m]

          a[m] <- ifelse(sm==1, tm, -Inf)
          b[m] <- ifelse(sm==1, Inf, tm)

          # Compute truncated moments
          E0  <- XI_FUNC(0, a[m], b[m], prior_i)
          E1  <- XI_FUNC(1, a[m], b[m], prior_i)

          if(is.nan(E1) | is.nan(E0)){
            browser()
          }
          # Compute integrals
          i4[m] <- ssm*um*E0
          i5[m] <- ifelse(um == 0, 1, ssm*(E1 - tm*E0))
        }

          A[i,]  <- a
          B[i,]  <- b
          I4[i,] <- i4
          I5[i,] <- i5
      }
    }
    #browser()
    #Reconstruct Constantine matrix
    Zf <- matrix(NA, nrow=mod$pdes)
    for(i in 1:mod$pdes){
      zi_curr <- coeff*I4[i,]
      for(k in (1:mod$pdes)[-i]){
        zi_curr <- zi_curr * I5[k,]
      }
      Zf[i] <- sum(zi_curr)
    }
    #if(scale01 == FALSE){
    #  # Transform back to native space
    #  Zf <- t(A_tform)%*%Zf%*%A_tform
    #}
    Zf_post[[r]] <- Zf
  }
  class(Zf_post) <- "ConstantineVector"
  if(length(Zf_post) == 1){
    return(Zf_post[[1]])
  }
  return(Zf_post)
}
