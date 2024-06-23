#' Convert functional BASS model to BASS model
#'
#' The argument to this function is the output of a \code{bass()} call when a single functional variable is specified using the \code{xx.func} argument.
#' Note that the resulting model may not be a valid bass object for some applications,
#' but the resulting model can be passed to \code{concordance::C_bass()} and related functions.
#'
#' @param bfm an object of class bass, where a functional variable has been specified.
#' @export
#' @examples
#' #The following are equivalent
#' n <- 100 #Number of observations
#' p <- 4   #Number of variables (beyond p = 2, variables are inert)
#' X <- matrix(runif(n*p), nrow=n)
#' y <- apply(X, 1, ff1)
#' gm <- gbass(X, Y, nmcmc=1000, nburn=901)
#' bm <- gm2bm(gm)
#' sob <- sobol(bm)
#' plot(sob)
#' @export
bassfunc2bass <- function(bfm){
  if(bfm$pfunc == 0){
    warning("No functional variable detected")
    return(bfm)
  }
  # First manipulate X and y into long form
  nx <- nrow(bfm$y)
  n <- length(bfm$y)
  p <- bfm$pdes + bfm$pfunc
  # XX <- matrix(NA, nrow=n, ncol=p)
  # yy <- rep(NA, n)
  tt <- bfm$xx.func
  nt <- length(tt)
  # for(i in seq_along(tt)){
  #     t_curr <- tt[i]
  #     indx <- ((i-1)*nx+1):(i*nx)
  #     XX[indx,] <- cbind(bfm$xx.des, rep(t_curr, nx))
  #     yy[indx] <- bfm$y[,i]
  # }

  nm <- bfm$n.models
  maxBasis <- max(bfm$nbasis)
  maxInt <- bfm$maxInt.des + bfm$maxInt.func

  # Combine information
  n.int <- bfm$n.int.des + bfm$n.int.func
  signs <- abind_custom(bfm$signs.des, bfm$signs.func, along=3)
  knots <- abind_custom(bfm$knotInd.des, bfm$knotInd.func, along=3)
  vars  <- abind_custom(bfm$vars.des, bfm$vars.func, along=3)
  range <- cbind(bfm$range.des, bfm$range.func)

  for(i in 1:nm){
    for(j in 1:maxBasis){
      # Fix knots so that they have data values, rather than indices
      knot_curr <- knots[i,j,]
      vars_curr <- vars[i,j,]
      ind <- ind2 <- which(!is.na(knot_curr))
      if(maxInt %in% ind){
        # add time variable
        knot_curr[maxInt] <- tt[knot_curr[maxInt]]
        ind2 <- ind[-length(ind)]
      }
      if(length(ind2) > 0){
        knot_curr[ind2] <- diag(bfm$xx.des[knot_curr[ind2], vars_curr[ind2],drop=FALSE])
      }
      knots[i,j,] <- knot_curr

      # Fix vars so that functional variable is assigned correct index
      if(!is.na(vars_curr[maxInt])){
        vars_curr[maxInt] <- p
      }
      vars[i,j,] <- vars_curr

      # Shift vectors so that NAs occur at the end
      nint_curr <- n.int[i,j]
      if(is.na(knot_curr[nint_curr])){
        knots[i,j,nint_curr] <- knots[i,j,maxInt]
        knots[i,j,maxInt]    <- NA

        signs[i,j,nint_curr] <- signs[i,j,maxInt]
        signs[i,j,maxInt]    <- NA

        vars[i,j,nint_curr] <- vars[i,j,maxInt]
        vars[i,j,maxInt]    <- NA
      }
    }
  }

  out <- bfm
  out$n.int.des <- n.int
  out$signs.des <- signs
  out$vars.des  <- vars
  out$knots.des <- knots # Notice that this is not knotInd.des
  out$pdes      <- p
  #out$pfunc     <- 0
  out$maxInt.des<- maxInt
  out$range.des <- range
  out$wasfunc   <- TRUE # Use this to flag (in C_bass) that knots are not indices.
  class(out)    <- c("bass")
  return(out)
}
