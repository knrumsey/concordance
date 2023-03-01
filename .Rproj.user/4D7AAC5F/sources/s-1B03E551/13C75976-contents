#' C matrix with Monte Carlo
#'
#' Approximates Constantine's C with Monte Carlo for a function f
#'
#' @param f the function f (or gradient of f, if grad=TRUE)
#' @param measure the number of inputs in f. See details for more sophisticated use (for non-uniform measure)
#' @param grad if TRUE f is assumed to return the gradient of f. When FALSE, forward diff is used for gradient approximation.
#' @param nmc the number of Monte Carlo replications
#' @param seed optional. seed for MC draws
#' @param ... additional arguments passed to f()
#' @return the approximated C matrix
#' @details `measure` should be an argument-free function which simulates a draw x ~ p(x) where p is the prior measure. If `measure` is numeric, then Monte Carlo draws are simulated from the standard uniform distribution as `runif(measure[1])`.
#' @export
C_mc <- function(f, measure, grad=FALSE, nmc=1e4, seed=NULL, ...){
  if(is.null(measure)){
    warning("measure not specified. Setting measure = 5, see help files for details")
    measure <- 5
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.numeric(measure)){
    n <- measure[1]
    measure <- function() runif(n)
  }
  if(grad){
    grad_f <- f
  }else{
    grad_f <- function(x, ...) fd_grad(f, x, ...)
  }
  n <- length(measure())

  Cf <- matrix(0, nrow=n, ncol=n)
  for(m in 1:nmc){
    x_m <- measure()
    del_f <- matrix(grad_f(x_m), nrow=n)
    Cf  <- Cf  +  tcrossprod(del_f, del_f)/nmc
  }
  return(Cf)
}



#' Estimate the Constantine Matrix with BASS
#'
#' Closed form estimator of the C matrix using a BASS model
#'
#' @param mod a fitted BASS model
#' @param prior NULL (default) (0,1])prior for each variable. See details for required structure of prior
#' @param mcmc.use set of indices telling which mcmc draws to use
#' @param scaled logical (default FALSE). When TRUE, the the C matix corresponds to the (0, 1)-scaled inputs rather than the original inputs.
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
C_bass <- function(mod, prior = NULL, mcmc.use=NULL, scaled=FALSE){
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
    prior[[i]]$trunc <- BASS:::scale.range(prior[[i]]$trunc, mod$range.des[,i])

    # 2. Handle ach distribution type separately
    distribution = prior[[i]]$dist
    if(distribution == "normal"){
      if(is.null(prior[[i]]$weights)){
        num_mix <- length(prior[[i]]$mean)
        prior[[i]]$weights <- rep(1/num_mix, num_mix)
      }
      prior[[i]]$mean <- BASS:::scale.range(prior[[i]]$mean, mod$range.des[,i])
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
  Cf_post <- list()
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
      knots      <- mod$knotInd.des[mod_number, 1:M, ]
      signs      <- mod$signs.des[mod_number, 1:M, ]
      indic      <- mod$vars.des[mod_number, 1:M, ]
      # Initalize arrays
      C <- A <- B <- I1 <- I2 <- I3 <- array(NA, dim=c(mod$pdes, M, M))

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
        if(!gbass_flag){
          d <- 1/((s + 1)/2 - s*t)
          s <- s*d
        }

        #Handle NA cases
        s[is.na(s)] <- 1
        t[is.na(t)] <- -Inf

        # Get integration constants
        c1 <- tcrossprod(s)
        C[i,,] <- c1

        a <- bp <- i1 <- i2 <- i3 <- matrix(0, M, M)
        for(m1 in 1:M){
          for(m2 in m1:M){
            um  <- u[c(m1, m2)]
            ssm <- s[c(m1, m2)]
            sm  <- sign(ssm)
            tm  <- t[c(m1, m2)]

            if(sm[1] == 1 & sm[2] == 1){
              a[m1,m2]  <- a[m2,m1]  <- max(tm)
              bp[m1,m2] <- bp[m2,m1] <- Inf
            }
            if(sm[1] == 1 & sm[2] == -1){
              a[m1,m2]  <- a[m2,m1]  <- tm[1]
              bp[m1,m2] <- bp[m2,m1] <- tm[2]
            }
            if(sm[1] == -1 & sm[2] == 1){
              a[m1,m2]  <- a[m2,m1]  <- tm[2]
              bp[m1,m2] <- bp[m2,m1] <- tm[1]
            }
            if(sm[1] == -1 & sm[2] == -1){
              a[m1,m2]  <- a[m2,m1]  <- -Inf
              bp[m1,m2] <- bp[m2,m1] <- min(tm)
            }
            aa <- a[m1,m2]
            bb <-  max(a[m1,m2], bp[m1,m2])

            # Compute integrals
            #E0 <- Efunc(0, aa, bb, prior_i)
            #E1 <- Efunc(1, aa, bb, prior_i)
            #E2 <- Efunc(2, aa, bb, prior_i)

            E0 <- XI_FUNC(0, aa, bb, prior_i)
            E1 <- XI_FUNC(1, aa, bb, prior_i)
            E2 <- NULL #Only compute when needed

            if(is.nan(E1) | is.nan(E0)){
              browser()
            }

            # Start with I3
            if(um[1] == 0 | um[2] == 0){
              i3[m1,m2] <- i3[m2,m1] <- 0
            }else{
              i3[m1,m2] <- i3[m2,m1] <- E0
            }

            #Next, I1 and I2
            if(um[1] == 0){
              if(um[2] == 0){
                i1[m1,m2] <- 0
                i1[m2,m1] <- 0
                i2[m1,m2] <- i2[m2,m1] <- E0
              }else{
                i1[m1,m2] <- 0
                i1[m2,m1] <- E0
                i2[m1,m2] <- i2[m2,m1] <- E1 - tm[2]*E0
              }
            }else{
              if(um[2] == 0){
                i2[m1,m2] <- i2[m2,m1] <- E1 - tm[1]*E0
                i1[m1,m2] <- E0
                i1[m2,m1] <- 0
              }else{
                E2 <- XI_FUNC(2, aa, bb, prior_i)
                i1[m1,m2] <- E1 - tm[2]*E0
                i1[m2,m1] <- E1 - tm[1]*E0
                i2[m1,m2] <- i2[m2,m1] <- E2 - (tm[1] + tm[2])*E1 + tm[1]*tm[2]*E0
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
      }
      #Add signs to the I's
      I1 <- I1*C
      I2 <- I2*C
      I3 <- I3*C
    }
    #Reconstruct Constantine matrix
    Cf <- matrix(NA, nrow=mod$pdes, ncol=mod$pdes)
    for(i in 1:mod$pdes){
      for(j in 1:i){
        if(i == j){
          cij_curr <- crossprod(coeff)*I3[i,,]
          for(k in (1:mod$pdes)[-c(i)]){
            cij_curr <- cij_curr * I2[k,,]
          }
          Cf[i,j] <- sum(cij_curr)
        }else{
          cij_curr <- crossprod(coeff)*I1[i,,]*t(I1[j,,])
          for(k in (1:mod$pdes)[-c(i,j)]){
            cij_curr <- cij_curr * I2[k,,]
          }
          Cf[i,j] <- Cf[j,i] <- sum(cij_curr)
        }
        # Matrix-y way of doing things (kinda)
        # This might save a little time when p is large, but the way I chose to code it (for time savings) has divide by zero problems if not careful
        # if(mod$pdes > 30){
        #   I2_prod <- matrix(1, M, M)
        #   for(k in (1:mod$pdes)){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   if(i == j){
        #     Cf[i,i] <- coeff%*%(I3[i,,] * I2_prod / (I2[i,,]+ 1e-12))%*%t(coeff)
        #   }else{
        #     Cf[i,j] <- Cf[j,i] <- coeff%*%(I1[i,,] * t(I1[j,,]) * I2_prod / (I2[i,,] + 1e-12) / (I2[j,,] + 1e-12))%*%t(coeff)
        #   }
        # }else{do it another way}
        #A third way to do it
        # I2_prod <- matrix(1, M, M)
        # if(i == j){
        #   for(k in (1:mod$pdes)[-i]){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   Iii <- I3[i,,] * I2_prod
        #   Cf[i,i] <- coeff%*%Iii%*%t(coeff)
        # }else{
        #   for(k in (1:mod$pdes)[-c(i,j)]){
        #     I2_prod <- I2_prod * I2[k,,]
        #   }
        #   Iij <- I1[i,,] * t(I1[j,,]) * I2_prod
        #   Cf[i,j] <- Cf[j,i] <- coeff%*%Iij%*%t(coeff)
        # }
      }
    }
    if(scaled == FALSE){
      Cf <- A_tform%*%Cf%*%t(A_tform)
    }
    Cf_post[[r]] <- Cf
  }
  class(Cf_post) <- "ConstantineMatrix"
  if(length(Cf_post) == 1){
    return(Cf_post[[1]])
  }
  return(Cf_post)
}


# Computes the truncated moments with respect to the prior
XI_FUNC <- function(pow, a, b, prior){
  astar <- max(a, prior$trunc[1])
  bstar <- max(astar, min(b, prior$trunc[2]))

  if(bstar > astar){
    # Number of components in mixture
    L <- length(prior$weights)
    res <- 0
    for(ell in 1:L){
      distribution <- prior$dist
      if(distribution == "uniform"){
        res <- res + XI_FUNC_UNIF(pow, astar, bstar, prior)*prior$weights[ell]
      }
      if(distribution == "normal"){
        res <- res + XI_FUNC_TNORM_SVW(pow, astar, bstar, prior)*prior$weights[ell]
      }
      if(distribution == "beta"){
        res <- res + XI_FUNC_BETA(pow, astar, bstar, prior)*prior$weights[ell]
      }
      if(distribution == "gamma"){
        res <- res + XI_FUNC_GAMMA(pow, astar, bstar, prior)*prior$weights[ell]
      }
    }
    return(res)
  }else{
    return(0)
  }
}

XI_FUNC_TNORM <- function(pow, a, b, prior){
  mu    <- prior$mean
  sigma <- prior$sd
  tau0  <- prior$trunc[1]
  tau1  <- prior$trunc[2]
  AA    <- (a-mu)/sigma
  BB    <- (b-mu)/sigma
  if(is.infinite(AA)) AA <- sign(AA)*1e4*sigma
  if(is.infinite(BB)) BB <- sign(BB)*1e4*sigma
  Z0    <- (pnorm(BB) - pnorm(AA))/(pnorm((tau1 - mu)/sigma) - pnorm((tau0-mu)/sigma))
  term  <- 1
  if(pow == 1){
    Z1 <- (dnorm(AA) - dnorm(BB))/(pnorm(BB) - pnorm(AA))
    term <- mu + sigma*Z1
  }
  if(pow == 2){
    Z1 <- (dnorm(AA) - dnorm(BB))/(pnorm(BB) - pnorm(AA))
    Z2 <- 1 - (AA*dnorm(AA) - BB*dnorm(BB))/(pnorm(BB) - pnorm(AA))
    term <- mu^2 + 2*sigma*mu*Z1 + sigma^2*Z2
  }
  return(Z0*term)
}


XI_FUNC_TNORM_SVW <- function(pow, a, b, prior){
  mu    <- prior$mean
  sigma <- prior$sd
  tau0  <- prior$trunc[1]
  tau1  <- prior$trunc[2]
  AA    <- (a-mu)/sigma
  BB    <- (b-mu)/sigma
  T1    <- (tau1-mu)/sigma
  T0    <- (tau0-mu)/sigma
  if(is.infinite(AA)) AA <- sign(AA)*1e6*sigma
  if(is.infinite(BB)) BB <- sign(BB)*1e6*sigma

  Delta1 <- pnorm(BB) - pnorm(AA)
  Delta2 <- pnorm(T1) - pnorm(T0)

  term  <- 0
  if(pow == 1){
    Z1     <- dnorm(AA) - dnorm(BB)
    term   <- sigma*Z1
  }
  if(pow == 2){
    Z1     <- dnorm(AA) - dnorm(BB)
    Z2     <- Delta1 - AA*dnorm(AA) - BB*dnorm(BB)
    term   <- 2*sigma*mu*Z1 + sigma^2*Z2
  }
  res <- Delta1*mu^pow/Delta2 + term/Delta2
  return(res)
}


XI_FUNC_UNIF <- function(pow, a, b, prior){
  res <- (b^(pow+1) - a^(pow+1))/(diff(prior$trunc)*(pow+1))
  return(res)
}

XI_FUNC_BETA <- function(pow, a, b, prior){
  shape1 <- prior$shape1
  shape2 <- prior$shape2
  num    <- zipfR::Ibeta(b, shape1 + pow, shape2) - zipfR::Ibeta(a, shape1 + pow, shape2)
  den    <- beta(shape1, shape2)
  return(num/den)
}

XI_FUNC_GAMMA <- function(pow, a, b, prior){
  alpha <- prior$scale
  beta  <- prior$rate
  num   <- zipfR::Igamma(alpha+pow, beta*b) - zipfR::Igamma(alpha+pow, beta*a)
  den   <- beta^pow*gamma(alpha)
  return(num/den)
}











