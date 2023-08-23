#' C_fg matrix with Monte Carlo
#'
#' Approximates generalized C matrix with Monte Carlo for functions f and g
#'
#' @param f the function f (or gradient of f, if grad=TRUE)
#' @param g the function g (or gradient of f, if grad=TRUE)
#' @param measure the number of inputs in f. See details for more sophisticated use (for non-uniform measure)
#' @param grad if TRUE f is assumed to return the gradient of f. When FALSE, forward diff is used for gradient approximation.
#' @param nmc the number of Monte Carlo replications
#' @param names (optional) names for the functions f and g
#' @param seed optional. seed for MC draws
#' @param return_C (default FALSE). When TRUE, the object returned is a list with components Cf Cg Cfg
#' @param ... additional arguments passed to f()
#' @return the approximated C matrix
#' @details `measure` should be an argument-free function which simulates a draw x ~ p(x) where p is the prior measure. Alternatively, `measure` can be a numeric scalar, in which case the Monte Carlo draws are simulated from the standard uniform distribution as `runif(measure[1])`.
#' @export
Cfg_mc <- function(f, g, measure, grad=FALSE, nmc=1e4, names = NULL, seed=NULL, return_C=FALSE, ...){
  if(is.null(measure)){
    warning("measure not specified. Setting measure = 5, type ?C_mc for details")
    measure <- 5
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.numeric(measure)){
    nn <- measure[1]
    measure <- function() runif(nn)
  }
  n <- length(measure())
  if(grad){
    grad_f <- f
    grad_g <- g
  }else{
    grad_f <- function(x, ...) fd_grad(f, x, ...)
    grad_g <- function(x, ...) fd_grad(g, x, ...)
  }
  Cf <- Cg <- Cfg  <- matrix(0, nrow=n, ncol=n)
  for(m in 1:nmc){
    x_m <- measure()
    del_f <- matrix(grad_f(x_m), nrow=n)
    del_g <- matrix(grad_g(x_m), nrow=n)
    Cf  <- Cf  +  tcrossprod(del_f, del_f)/nmc
    Cg  <- Cg  +  tcrossprod(del_g, del_g)/nmc
    Cfg <- Cfg +  tcrossprod(del_f, del_g)/nmc
  }
  if(return_C){
    out <- list(Cfg=Cfg, Cf=Cf, Cg=Cg)
  }
  return(Cf)
}





#' Estimate Cfg with BASS
#'
#' Closed form estimator of the Cfg matrix using a BASS model
#'
#' @param mod1 a fitted BASS model for first function
#' @param mod2 a fitted BASS model for second function
#' @param prior NULL (default) `[0, 1]` prior for each variable. See details for required structure of prior
#' @param mcmc.use vector of mcmc indices to be used for both models. Otherwise, a matrix
#' @return A list representing the posterior distribution of the Co-Constantine matrix (Cfg).
#' @details prior should be a list of length p (one object for each variable). Each element of prior should be a named list with fields
#'
#'    dist - ("uniform", "normal").
#'
#'    trunc - truncation bounds (a, b)
#'
#'    mean - vector of means (mixture of normals only)
#'
#'    sigma - vector of sds (mixture of normals only)
#'
#'    weights - vector of mixture weights (mixture of normals only)
#' @export
Cfg_bass <- function(mod1, mod2, prior = NULL, mcmc.use=NULL){
  mod <- mod1
  if(is.null(mcmc.use)){
    mcmc.use <- min(length(mod$s2), length(mod2$s2))
  }
  mcmc.use <- as.matrix(mcmc.use)
  if(ncol(mcmc.use) == 1){
    mcmc.use <- cbind(mcmc.use, mcmc.use)
  }
  if(ncol(mcmc.use) > 2){
    warning("ncol(mcmc.use) should not exceed 2")
  }
  if(mod$p != mod2$p){
    stop("Detected different number of variables in mod1 and mod2")
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

    # 2. Handle ach distribution type separately
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
      prior[[i]]$weights <- prior[[i]]$weights/cc # DF: prior[[i]]$z # change weights with truncation # divide by cc instead to keep the same prior shape# does the truncation change the distribution shape in the non-truncated regions??

      #Check for extrapolation
      qq <- qnorm(c(0.0005, 0.9995), prior[[i]]$mean, prior[[i]]$sd)
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

  # Make constantine matrix

  Cfg_post <- list()
  Xt <- mod$xx.des
  for(r in 1:length(mcmc.use)){
    #Compute only the stuff we will need for every iteration
    rr1 <- mcmc.use[r,1]
    rr2 <- mcmc.use[r,2]

    mod_number_new <- mod$model.lookup[rr1]
    coeff          <- mod$beta[rr1,]
    coeff          <- matrix(coeff[!is.na(coeff)][-1], nrow=1)
    M_new          <- length(coeff)

    mod_number_new2 <- mod2$model.lookup[rr2]
    coeff2          <- mod2$beta[rr2,]
    coeff2          <- matrix(coeff2[!is.na(coeff2)][-1], nrow=1)
    M_new2          <- length(coeff2)

    compute_flag <- FALSE
    if(r == 1){
      compute_flag <- TRUE
    }else{
      if((mod_number != mod_number_new) | (mod_number2 != mod_number_new2)){
        compute_flag <- TRUE
      }
    }

    if(compute_flag){
      mod_number <- mod_number_new
      M          <- M_new
      knots      <- mod$knotInd.des[mod_number, 1:M, ]
      signs      <- mod$signs.des[mod_number, 1:M, ]
      indic      <- mod$vars.des[mod_number, 1:M, ]
      if(M==1){
        signs <- matrix(signs, nrow=1, ncol=length(signs))
        indic <- matrix(indic, nrow=1, ncol=length(indic))
        knots <- matrix(knots, nrow=1, ncol=length(knots))
      }

      mod_number2 <- mod_number_new2
      M2          <- M_new2
      knots2      <- mod2$knotInd.des[mod_number2, 1:M2, ]
      signs2      <- mod2$signs.des[mod_number2, 1:M2, ]
      indic2      <- mod2$vars.des[mod_number2, 1:M2, ]
      if(M2==1){
        signs2 <- matrix(signs2, nrow=1, ncol=length(signs2))
        indic2 <- matrix(indic2, nrow=1, ncol=length(indic2))
        knots2 <- matrix(knots2, nrow=1, ncol=length(knots2))
      }

      # Initalize arrays
      C <- A <- B <- I1 <- I2 <- I3 <- array(NA, dim=c(mod$p, M, M2))
      I1b <- array(NA, dim=c(mod$p, M2, M))
      for(i in 1:mod$p){
        prior_i <- prior[[i]]
        v <- apply(indic, 1, function(zz) match(i, zz))
        u <- !is.na(v)
        s <- apply(cbind(signs, v), 1, function(zz) zz[zz[mod$maxInt.des + 1]])
        t <- Xt[apply(cbind(knots, v), 1, function(zz) zz[zz[mod$maxInt.des + 1]]), i]

        v2 <- apply(indic2, 1, function(zz) match(i, zz))
        u2 <- !is.na(v2)
        s2 <- apply(cbind(signs2, v2), 1, function(zz) zz[zz[mod2$maxInt.des + 1]])
        t2 <- Xt[apply(cbind(knots2, v2), 1, function(zz) zz[zz[mod2$maxInt.des + 1]]), i]

        # If this comes from BASS (rather than GBASS with gm2bm)
        # then we need to account for Devin's g-scaling-factors
        if(!("gbass" %in% class(mod))){
          d <- 1/((s + 1)/2 - s*t)
          s <- s*d
        }
        if(!("gbass" %in% class(mod2))){
          d <- 1/((s2 + 1)/2 - s2*t2)
          s2 <- s2*d
        }
        #NOTE!!! DOES THE ABOVE WORK? CAN I HAVE A GBASS AND A BASS MODEL?

        #Handle NA cases
        s[is.na(s)] <- 1
        t[is.na(t)] <- -Inf

        s2[is.na(s2)] <- 1
        t2[is.na(t2)] <- -Inf

        # Get integration constants
        cc <- tcrossprod(s, s2)
        C[i,,] <- cc

        a <- bp <- i1 <- i2 <- i3 <- matrix(NA, M, M2)
        i1b <- a2 <- bp2 <- matrix(NA, M2, M)
        for(m1 in 1:M){
          for(m2 in 1:M2){
            #Main change for Cfg here (three lines)
            um <- c(u[m1], u2[m2])
            sm <- c(s[m1], s2[m2])
            tm <- c(t[m1], t2[m2])

            signm = sign(sm)

            if(signm[1] == 1 & signm[2] == 1){
              a[m1,m2]   <- max(tm)
              bp[m1,m2]  <- Inf
            }
            if(signm[1] == 1 & signm[2] == -1){
              a[m1,m2]   <- tm[1]
              bp[m1,m2]  <- tm[2]
            }
            if(signm[1] == -1 & signm[2] == 1){
              a[m1,m2]   <- tm[2]
              bp[m1,m2]  <- tm[1]
            }
            if(signm[1] == -1 & signm[2] == -1){
              a[m1,m2]  <- -Inf
              bp[m1,m2] <- min(tm)
            }
            aa <- a[m1,m2]
            bb <-  max(aa, bp[m1,m2])

            # Compute integrals
            #E0 <- Efunc(0, aa, bb, prior_i)
            #E1 <- Efunc(1, aa, bb, prior_i)
            #E2 <- Efunc(2, aa, bb, prior_i)

            E0 <- XI_FUNC(0, aa, bb, prior_i)
            E1 <- XI_FUNC(1, aa, bb, prior_i)
            E2 <- XI_FUNC(2, aa, bb, prior_i)

            # Start with I3
            if(um[1] == 0 | um[2] == 0){
              i3[m1,m2] <- 0
            }else{
              i3[m1,m2] <- E0
            }

            #Next, I1 and I2
            if(um[1] == 0){
              if(um[2] == 0){
                i1[m1,m2]  <- 0
                i1b[m2,m1] <- 0
                i2[m1,m2]  <- E0
              }else{
                i1[m1,m2]  <- 0
                i1b[m2,m1] <- E0
                i2[m1,m2]  <- E1 - tm[2]*E0
              }
            }else{
              if(um[2] == 0){
                i1[m1,m2]  <- E0
                i1b[m2,m1] <- 0
                i2[m1,m2]  <- E1 - tm[1]*E0
              }else{
                i1[m1,m2]  <- E1 - tm[2]*E0
                i1b[m2,m1] <- E1 - tm[1]*E0
                i2[m1,m2]  <- E2 - (tm[1] + tm[2])*E1 + tm[1]*tm[2]*E0
              }
            }
          }
        }
        b <- pmax(a, bp)
        A[i,,] <- a
        B[i,,] <- b
        I1[i,,] <- i1
        I2[i,,] <- i2
        I3[i,,] <- i3
        I1b[i,,] <- i1b
      }
      #Add signs to the I's
      I1 <- I1*C
      I2 <- I2*C
      I3 <- I3*C
      for(iii in 1:mod$p){
        I1b[iii,,] <- I1b[iii,,] * t(C[iii,,])
      }

    }
    #Reconstruct Constantine matrix
    Cfg <- matrix(NA, nrow=mod$p, ncol=mod$p)
    for(i in 1:mod$p){
      for(j in 1:mod$p){
        # #Naive way for now
        cij_curr <- 0
        if(i != j){
          for(m1 in 1:M){
            for(m2 in 1:M2){
              term <- coeff[m1]*coeff2[m2]*I1[i,m1,m2]*I1b[j,m2,m1]
              for(k in (1:mod$p)[-c(i,j)]){
                term <- term * I2[k,m1,m2]
              }
              cij_curr <- cij_curr + term
            }
          }
        }else{
          for(m1 in 1:M){
            for(m2 in 1:M2){
              term <- coeff[m1]*coeff2[m2]*I3[i,m1,m2]
              for(k in (1:mod$p)[-i]){
                term <- term * I2[k,m1,m2]
              }
              cij_curr <- cij_curr + term
            }
          }
        }

        Cfg[i,j] <- cij_curr


        # if(i == j){
        #   cij_curr <- crossprod(coeff)*I3[i,,]
        #   for(k in (1:mod$p)[-c(i)]){
        #     cij_curr <- cij_curr * I2[k,,]
        #   }
        #   Cf[i,j] <- sum(cij_curr)
        # }else{
        #   cij_curr <- crossprod(coeff)*I1[i,,]*t(I1[j,,])
        #   for(k in (1:mod$p)[-c(i,j)]){
        #     cij_curr <- cij_curr * I2[k,,]
        #   }
        #   Cf[i,j] <- Cf[j,i] <- sum(cij_curr)
        # }
        # Matrix-y way of doing things (kinda)
        # This might save a little time when p is large, but the way I chose to code it (for time savings) has divide by zero problems if not careful
        # if(mod$p > 30){
        #   I2_prod <- matrix(1, M, M)
        #   for(k in (1:mod$p)){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   if(i == j){
        #     Cf[i,i] <- coeff%*%(I3[i,,] * I2_prod / (I2[i,,]+ 1e-12))%*%t(coeff)
        #   }else{
        #     Cf[i,j] <- Cf[j,i] <- coeff%*%(I1[i,,] * t(I1[j,,]) * I2_prod / (I2[i,,] + 1e-12) / (I2[j,,] + 1e-12))%*%t(coeff)
        #   }
        # }else{do it another way}
        #A third way to do it
        # I2_prod <- matrix(1, M, M2)
        # if(i == j){
        #   for(k in (1:mod$p)[-i]){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   Iii <- I3[i,,] * I2_prod
        #   Cfg[i,i] <- coeff%*%Iii%*%t(coeff2)
        # }else{
        #   for(k in (1:mod$p)[-c(i,j)]){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   Iij <- I1[i,,] * t(I1b[j,,]) * I2_prod
        #   Cfg[i,j] <- coeff%*%Iij%*%t(coeff2)
        # }

      }
    }
    Cfg_post[[r]] <- Cfg
  }
  class(Cfg_post) <- "CoConstantineMatrix"
  if(length(Cfg_post) == 1){
    return(Cfg_post[[1]])
  }
  return(Cfg_post)
}
