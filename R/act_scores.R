#' Activity Scores
#'
#' This function computes the activity scores for main effects of the variables
#'
#' @param C Constantines C matrix (e.g. from C_bass, C_mc, or C_gp)
#' @param k The number of columns of W to consider
#' @param plt Logical, should a plot be made?
#' @param norm Logical, should activity scores be normalized to have a maximum value of 1?
#' @return the activity scores
#' @export
act_scores <- function(C, k=1, plt=FALSE, norm=FALSE){
  eig <- eigen(C)
  W <- eig$vectors
  lam <- eig$values
  if(k == 1){
    res <- W[,k]^2*lam[k]
  }else{
    res <- apply(W[,1:k], 1, function(w, lam) sum(w^2*pmax(lam, 0)), lam=lam[1:k])
  }

  if(norm){
    res <- res/max(res)
  }

  if(plt){
    plot(res, xlab="Inputs", ylab=paste0("Activity Score (", k, ")"), pch=16, cex=2, ylim=c(0, max(res)))
  }

  return(res)
}

#' Co-Activity Scores
#'
#' This function computes the activity scores for main effects of the variables
#'
#' @param V The symmetrized co-Constantine matrix from Cfg_bass()
#' @param q The number of columns of W to consider
#' @param signed Use signed or unsigned version?
#' @param plt Logical, should a plot be made?
#' @param norm Logical, should activity scores be normalized to have a maximum value of 1?
#' @return the coactivityactivity scores
#' @export
coact_scores <- function(V, q=1, signed=TRUE, plt=FALSE, norm=FALSE){
  eig <- eigen(V)
  W <- eig$vectors
  lam <- eig$values
  ord <- rev(order(abs(lam)))
  if(signed == FALSE){
    lam <- abs(lam)
  }
  if(q == 1){
    res <- W[,ord[q]]^2*lam[ord[q]]
  }else{
    ind <- ord[1:q]
    res <- apply(W[,ind], 1, function(w, lam) sum(w^2*lam), lam=lam[1:q])
  }

  if(norm){
    res <- res/max(abs(res))
  }

  if(plt){
    plot(res, xlab="Inputs", ylab=paste0("Activity Score (", q, ")"), pch=16, cex=2, ylim=range(res))
  }

  return(res)
}


#' Active Dimension (Not validated. Might be buggy)
#'
#' This function estimates the dimensions of the active subspace using a sequential testing approach
#'
#' @param C Constantines C matrix (e.g. from C_bass, C_mc, or C_gp)
#' @param X the original input variables
#' @param y the original response variable (mod$y when using C_bass(mod))
#' @param k The maximum number of columns of W to consider
#' @param alpha significance threshold for testing procedure
#' @param all_sets should all dimension sets be returned? Or just the smallest set.
#' @param verbose should progress be printed
#' @return a list of active subspace dimensions
#' @export
act_dims <- function(C, X, y, k=ncol(C), alpha = 0.05, all_sets=TRUE, verbose=TRUE){
  EX <- eigen(C)$vectors
  T1 <- as.matrix(X)%*%EX

  #Fitting linear models
  if(verbose) cat("Creating linear models\n")
  fit <- list()
  fit[[1]] <- lm(y~1)
  for(i in 1:k){
    fit[[i+1]] <- lm(y~T1[,1:i])
  }

  res <- list()
  cnt <- 1
  tmp <- anova(fit[[1]], fit[[2]])$`Pr(>F)`[2]
  if(tmp < alpha){
    res[[cnt]] <- 1
    cnt <- cnt + 1
    kk <- 2
    if(!all_sets) return(res)
  }else{
    kk <- 1
  }

  for(i in 3:(k+1)){
    tmp <- anova(fit[[kk]], fit[[i]])$`Pr(>F)`[2]
    if(tmp < alpha){
      res[[cnt]] <- i-1
      cnt <- cnt + 1
      kk <- i
      if(!all_sets) return(res)
    }
  }
  return(res)
}








