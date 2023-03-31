#' Build Prior Method for C_bass and Cfg_bass
#'
#' A quick way to build priors for use in `C_bass` and `Cfg_bass`.
#' For more complicated priors, such as mixture distributions, see details in `?C_bass`.
#'
#'
#' @param dist A vector of length p. Valid entries include "uniform", "normal", "beta", "gamma".
#' @param trunc A matrix of dimension px2 (rows are recycled if nrow < p). `Inf` is a valid entry.
#' @param mean A p-vector of means (used for normal only)
#' @param sd A p-vector of sds (used for normal only)
#' @param shape1 A p-vector of shape1 parameters for beta prior
#' @param shape2 A p-vector of shape2 parameters for beta prior
#' @param shape A p-vector of shape parameters for gamma prior
#' @param scale A p-vector of scale parameters for gamma prior
#' @return a list which can be passed into C_bass or Cfg_bass as a prior.
#' @details All vectors and matrix rows are recycled for parameters. The vector `dist` cannot be recyled as it defines p.
#' @export
build_prior <- function(dist, trunc=NULL,
                        mean=NULL, sd=NULL,
                        shape1=NULL, shape2=NULL,
                        shape=NULL, scale=NULL){
  p <- length(dist)
  prior <- list()
  for(i in 1:p){
    pr <- list(dist=dist[i])
    if(!is.null(trunc)){
      q <- nrow(trunc)
      j <- ((i - 1) %% q) + 1
      pr$trunc <- trunc[j,]
    }
    if(!is.null(mean)){
      q <- length(mean)
      j <- ((i - 1) %% q) + 1
      pr$mean <- mean[j]
    }
    if(!is.null(sd)){
      q <- length(sd)
      j <- ((i - 1) %% q) + 1
      pr$sd <- sd[j]
    }

    if(!is.null(shape1)){
      q <- length(shape1)
      j <- ((i - 1) %% q) + 1
      pr$shape1 <- shape1[j]
    }
    if(!is.null(shape2)){
      q <- length(shape2)
      j <- ((i - 1) %% q) + 1
      pr$shape2 <- shape2[j]
    }

    if(!is.null(shape)){
      q <- length(shape)
      j <- ((i - 1) %% q) + 1
      pr$shape <- shape[j]
    }
    if(!is.null(scale)){
      q <- length(scale)
      j <- ((i - 1) %% q) + 1
      pr$scale <- scale[j]
    }
    prior[[i]] <- pr
  }
  return(prior)
}


