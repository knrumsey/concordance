#' Forward diff function
#'
#' Function for approximating the gradient of a function
#'
#' @param f the function to find the gradient of
#' @param x the input values
#' @param h the tolerance
#' @param ... additional inputs to be passed to f
#' @return The approximate gradient of f at x
#' @export
fd_grad <- function(f, x, h=1e-12, ...){
  gradient <- rep(NA, length(x))
  for(i in seq_along(x)){
    x[i] <- x[i] + h*1i
    gradient[i] <- Im(f(x, ...))/h
    x[i] <- x[i] - h*1i
  }
  return(gradient)
}


#' Concordance Analysis (with BASS)
#'
#' Computes the concordance between mod1 and mod2 (BASS models representing f1 and f2)
#'
#' @param mod1 BASS model representing function 1
#' @param mod2 BASS model representing function 2
#' @param prior NULL (default) Uniform(0,1) prior for each variable. See details fr required prior structure.
#' @param mcmc.use a vector of mcmc replications to use. Can also be a 2-column matrix with indices for f1 and f2.
#' @param q order for the activity score measures
#' @param ... additional arguments passed fd_grad()
#' @return Estimates of C1, C2, C12, V12, conc(f1, f2), contributions and coactivity scores
#' @details `measure` should be an argument-free function which simulates a draw x ~ p(x) where p is the prior measure. If `measure` is numeric, then Monte Carlo draws are simulated from the standard uniform distribution as `runif(measure[1])`.
#' @export
coactive_bass <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, q=1, ...){
  if(is.null(mcmc.use)) mcmc.use <- min(length(mod1$s2), length(mod2$s2))
  mcmc.use <- as.matrix(mcmc.use)
  if(ncol(mcmc.use)==1) mcmc.use <- cbind(mcmc.use, mcmc.use)
  obj <- list()
  obj$C1 <- C_bass(mod1, prior, mcmc.use[,1])
  obj$C2 <- C_bass(mod2, prior, mcmc.use[,2])
  obj$C12 <- Cfg_bass(mod1, mod2, prior, mcmc.use)

  V12 <- t1 <- t2 <- t12 <- conc <- activity1 <- activity2 <- coactivity_signed <- coactivity_unsigned <- contrib <- list()
  nn <- nrow(mcmc.use)
  for(i in 1:nn){
    if(nn == 1){
      C1i <- obj$C1
      C2i <- obj$C2
      C12i <- obj$C12
    }else{
      C1i <- obj$C1[[i]]
      C2i <- obj$C2[[i]]
      C12i <- obj$C12[[i]]
    }

    V12[[i]] <- (C12i + t(C12i))/2
    t1[[i]] <- sum(diag(C1i))
    t2[[i]] <- sum(diag(C2i))
    t12[[i]] <- sum(diag(C12i))
    conc[[i]] <- t12[[i]]/sqrt(t1[[i]]*t2[[i]])
    activity1[[i]] <- act_scores(C1i, q)
    activity2[[i]] <- act_scores(C2i, q)
    coactivity_signed[[i]] <- coact_scores(V12[[i]], q=q, signed=TRUE)
    coactivity_unsigned[[i]] <- coact_scores(V12[[i]], q=q, signed=FALSE)
    contrib[[i]] <- eigen(V12[[i]])$values/sqrt(t1[[i]]*t2[[i]])
  }
  if(nn == 1){
    V12 <- V12[[1]]
    t1 <- t1[[1]]
    t2 <- t2[[1]]
    t12 <- t12[[1]]
    conc <- conc[[1]]
    contrib <- contrib[[1]]
    activity1 <- activity1[[1]]
    activity2 <- activity2[[1]]
    coactivity_signed <- coactivity_signed[[1]]
    coactivity_unsigned <- coactivity_unsigned[[1]]
  }
  obj$V12 <- V12
  obj$t1 <- t1
  obj$t2 <- t2
  obj$t12 <- t12
  obj$conc <- conc
  obj$contrib <- contrib
  obj$activity1 <- activity1
  obj$activity2 <- activity2
  obj$coactivity_signed <- coactivity_signed
  obj$coactivity_unsigned <- coactivity_unsigned

  class(obj) <- "CoactiveSubspace"
  return(obj)
}

CoactiveSubspace.print <- function(x, ...){
  cat("conc(f1, f2) = ", mean(unlist(x$conc)))
}


#' Concordance
#'
#' Computes the concordance between f and g
#'
#' @param f the function f (or gradient of f, if grad=TRUE)
#' @param g the function g (or gradient of g, if grad=TRUE)
#' @param measure the number of inputs in f and g. See details for more sophisticated use (for non-uniform measure)
#' @param grad if TRUE f and g are assumed to return gradients. When FALSE, forward diff is used for approximation.
#' @param nmc the number of monte carlo replications
#' @param ... additional arguments passed fd_grad()
#' @return the concordance between functions f and g
#' @details `measure` should be an argument-free function which simulates a draw x ~ p(x) where p is the prior measure. If `measure` is numeric, then Monte Carlo draws are simulated from the standard uniform distribution as `runif(measure[1])`.
#' @export
conc_mc <- function(f, g, measure, grad=FALSE, nmc=1e4, ...){
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
  tf <- tg <- tfg <- 0
  for(m in 1:nmc){
    x_m <- measure()
    del_f <- matrix(grad_f(x_m), nrow=n)
    del_g <- matrix(grad_g(x_m), nrow=n)
    tf <- tf +  crossprod(del_f, del_f)/nmc
    tg <- tg +  crossprod(del_g, del_g)/nmc
    tfg <- tfg +  crossprod(del_f, del_g)/nmc
  }
  return(as.numeric(tfg/sqrt(tf*tg)))
}

#' Concordance Analysis
#'
#' Performs a full concordance analysis between f and g
#'
#' @param f the function f (or gradient of f, if grad=TRUE)
#' @param g the function g (or gradient of g, if grad=TRUE)
#' @param measure the number of inputs in f and g. See details for more sophisticated use (for non-uniform  measure)
#' @param grad if TRUE f and g are assumed to return gradients. When FALSE, forward diff is used for approximation.
#' @param nmc the number of monte carlo replications
#' @param names names of the functions
#' @param seed optional. seed for monte carlo draws
#' @param ... additional arguments passed fd_grad()
#' @return a list with components: C (constantine matrices), principle_grads, contributions, totals, conc, dist
#' @details `measure` should be an argument-free function which simulates a draw x ~ p(x) where p is the prior measure. If `measure` is numeric, then Monte Carlo draws are simulated from the standard uniform distribution as `runif(measure[1])`.
#' @export
conc_analysis_mc <- function(f, g, measure, grad=FALSE, nmc=1e4, names=c("f", "g"), seed=NULL, ...){
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
  Cf <- Cg <- Cfg <- Cgf <- matrix(0, nrow=n, ncol=n)
  for(m in 1:nmc){
    x_m <- runif(n)
    del_f <- matrix(grad_f(x_m), nrow=n)
    del_g <- matrix(grad_g(x_m), nrow=n)
    Cf  <- Cf  +  tcrossprod(del_f, del_f)/nmc
    Cg  <- Cg  +  tcrossprod(del_g, del_g)/nmc
    Cfg <- Cfg +  tcrossprod(del_f, del_g)/nmc
    #Cgf <- Cgf +  tcrossprod(del_g, del_f)/nmc
  }
  Cgf <- t(Cfg)
  Vfg <- (Cfg + Cgf)/2

  # Get principle gradients
  Ef <- eigen(Cf)
  deltaf <- Ef$vectors%*%diag(sqrt(abs(Ef$values)))

  Eg <- eigen(Cg)
  deltag <- Eg$vectors%*%diag(sqrt(abs(Eg$values)))

  Efg <- eigen(Vfg)
  deltafg <- Efg$vectors%*%diag(sqrt(abs(Efg$values)))

  #Get contributions
  pif <- Ef$values/sum(Ef$values)
  pig <- Eg$values/sum(Eg$values)
  pifg <- Efg$values/sqrt(sum(Ef$values)*sum(Eg$values))

  #Return object
  out <- list(C = list(Cf=Cf, Cg=Cg, Cfg=Cfg, Vfg=Vfg),
              principle_grads = list(delta_f=deltaf, delta_g=deltag, delta_fg = deltafg),
              contributions = list(pi_f = pif, pi_g = pig, pi_fg = pifg),
              totals = list(t_f = sum(Ef$values), t_g = sum(Eg$values), t_fg = sum(Efg$values)),
              conc = sum(Efg$values)/sqrt(sum(Ef$values)*sum(Eg$values)),
              names=names)
  out$dist = sqrt(1 - out$conc)
  out$names <- names
  class(out) <- "ConcordanceAnalysis"
  return(out)
}

#' Summary and Print functions
#'
#' Prints a summary for an object of class "ConcordanceAnalysis"
#'
#' @param x object of class "ConcordanceAnalysis"
#' @param ... Additional arguments (ignored)
#' @export
print.ConcordanceAnalysis <- function(x, ...){
  obj <- x
  names <- obj$names

  sf <-  apply(abs(obj$principle_grads$delta_f%*%diag(obj$contributions$pi_f)), 1, sum)
  sg <-  apply(abs(obj$principle_grads$delta_g%*%diag(obj$contributions$pi_g)), 1, sum)
  sfg <- apply(abs(obj$principle_grads$delta_fg%*%diag(obj$contributions$pi_fg)), 1, sum)

  cat("Concordance analysis for functions", names[1], "and", names[2],
      "\n   Concordance:  ", round(obj$conc, 6),
      "\n   Discord:      ", round(obj$dist, 6),
      "\n\n   ", names[1],
      "\n   Root Total:   ", round(sqrt(obj$totals$t_f), 3),
      "\n   Contributions:", round(obj$contributions$pi_f, 3),
      "\n   Sensitivity:  ", round(sf/max(abs(sf)), 3),
      "\n\n   ", names[2],
      "\n   Root Total:   ", round(sqrt(obj$totals$t_g), 3),
      "\n   Contributions:", round(obj$contributions$pi_g, 3),
      "\n   Sensitivity:  ", round(sg/max(abs(sg)), 3),
      "\n\n   ", names[1], "*", names[2],
      "\n   Root Co-Total:", round(sqrt(obj$totals$t_fg), 3),
      "\n   Contributions:", round(obj$contributions$pi_fg, 3),
      "\n   Sensitivity:  ", round(sfg/max(abs(sfg)), 3),
      sep=" ")
}

#' Summary and Print functions
#'
#' Prints a summary for an object of class "ConcordanceAnalysis"
#'
#' @param object object of class "ConcordanceAnalysis"
#' @param ... Ignored
#' @export
summary.ConcordanceAnalysis <- function(object, ...){
  print.ConcordanceAnalysis(object)
}

#' Plot contributions
#'
#' Plots the contributions pi_f, pi_g, and pi_fg
#'
#' @param obj object of class "ConcordanceAnalysis"
#' @param ... additional arguments passed to barplot
#' @export
plot_contributions <- function(obj, ...){
  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  plot(NULL, ylim=c(min(c(0, min(unlist(obj$contributions)))), max(unlist(obj$contributions))), xlim=c(1,length(obj$contributions$pi_f)),
                    xlab="Active Direction Index", ylab="Contribution", main="Contribution of Active Gradients", ...)
  points(obj$contributions$pi_f, pch=1, col=pal[1], cex=2)
  points(obj$contributions$pi_g, pch=2, col=pal[2], cex=1.8)
  points(obj$contributions$pi_fg, pch=3, col=pal[3], cex=1.9)
  abline(h=sum(obj$contributions$pi_fg), col=pal[4], lty=3)
  legend("topright", legend=c(obj$names, paste0(obj$names[1], "*", obj$names[2])), pch=1:3, col=pal, cex=1.5, bty='n')
}

#' Plot sensitivities
#'
#' The sensitivity of variable j (wrt to f) is defined as
#' sum_{i=1}^n pi_f(i) * delta_f(j,i)
#' The definition is similar for g or for fg
#'
#' @param obj object of class "ConcordanceAnalysis"
#' @param vnames optional vector of variable names
#' @param ... additional arguments passed to barplot
#' @export
plot_sensitivities <- function(obj, vnames=NULL, ...){
  sf <-  apply(abs(obj$principle_grads$delta_f%*%diag(obj$contributions$pi_f)), 1, sum)
  sg <-  apply(abs(obj$principle_grads$delta_g%*%diag(obj$contributions$pi_g)), 1, sum)
  sfg <- apply(abs(obj$principle_grads$delta_fg%*%diag(obj$contributions$pi_fg)), 1, sum)
  S <- rbind(sf, sg, sfg)
  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  if(is.null(vnames))
    vnames <- 1:ncol(S)
  barplot(S, beside=TRUE, names.arg=vnames, col=pal[1:3], main="Variable Sensitivity", ...)
  #legend("topright", legend=c(obj$names, paste0(obj$names[1], "*", obj$names[2])), fill=pal, cex=1.5, bty='y')
}

#' Plot components of the kth Principle Gradient
#'
#' @param obj object of class "ConcordanceAnalysis"
#' @param k which principle gradient is desired?
#' @param vnames optional vector of variable names
#' @param ... additional arguments passed to barplot
#' @export
plot_active_grad_k <- function(obj, k = 1, vnames=NULL, ...){
  grad <- rbind(obj$principle_grads$delta_f[,k], obj$principle_grads$delta_g[,k], obj$principle_grads$delta_fg[,k])
  if(is.null(vnames))
    vnames <- 1:ncol(grad)

  if(k == 1){
    ord_suffix <- "st"
  }
  if(k == 2){
    ord_suffix <- "nd"
  }
  if(k == 3){
    ord_suffix <- "rd"
  }
  if(k > 3){
    ord_suffix <- "th"
  }
  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  barplot(grad, beside=TRUE, names.arg=vnames, col=pal[1:3], main=paste0(k, ord_suffix, " Principle Gradient"), ...)
  #legend("topright", legend=c(obj$names, paste0(obj$names[1], "*", obj$names[2])), fill=pal, cex=1.5, bty='y')
}

#' Plotting Function for object of class ConcordanceAnalysis
#'
#' @param x object of class "ConcordanceAnalysis"
#' @param ... arguments to be passed to individual plot functions
#' @export
plot.ConcordanceAnalysis <- function(x, ...){
  obj <- x
  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
  par(mfrow=c(2,2))
  plot_contributions(obj)
  plot_sensitivities(obj, ...)
  plot_active_grad_k(obj, 1, ...)
  plot_active_grad_k(obj, 2, ...)
  par(mfrow=c(1,1))
}








