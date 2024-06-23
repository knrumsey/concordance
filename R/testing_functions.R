#' Trace of a matrix
#'
#' Shortcut for sum(diag(A))
#'
#' @param A a matrix
#' @return The trace of a matrix
#' @export
tr <- function(A) sum(diag(A))

#' The Borehole Function
#'
#' This function models the flow of water through a borehole.
#'
#' @param xx the 7 inputs, restricted to the unit interval
#' @param design the radius of the borehole (typically rw)
#' @details PARAMETER RANGES
#' rw in `[0.05, 0.15]` 	radius of borehole (m)
#' r in `[100, 50000]` 	radius of influence (m)
#' Tu in `[63070, 115600]`    	transmissivity of upper aquifer (m2/yr)
#' Hu in `[990, 1110]` 	potentiometric head of upper aquifer (m)
#' Tl in `[63.1, 116]` 	transmissivity of lower aquifer (m2/yr)
#' Hl in `[700, 820]` 	potentiometric head of lower aquifer (m)
#' L in `[1120, 1680]` 	length of borehole (m)
#' Kw in `[9855, 12045]` 	hydraulic conductivity of borehole (m/yr)
#' @return The output of the borehole function
#' @export
f_borehole <- function(xx, design=0.5){
  rr<-matrix(c(
    0.01, 0.9,
    100, 50000,
    63070, 115600,
    990, 1110,
    63.1, 116,
    700, 820,
    1120, 1680,
    9855, 12045
  ),nrow=2)

  xx <- c(design, xx)

  rw <- xx[1]*diff(rr[,1]) + rr[1,1]
  r  <- xx[2]*diff(rr[,2]) + rr[1,2]
  Tu <- xx[3]*diff(rr[,3]) + rr[1,3]
  Hu <- xx[4]*diff(rr[,4]) + rr[1,4]
  Tl <- xx[5]*diff(rr[,5]) + rr[1,5]
  Hl <- xx[6]*diff(rr[,6]) + rr[1,6]
  L  <- xx[7]*diff(rr[,7]) + rr[1,7]
  Kw <- xx[8]*diff(rr[,8]) + rr[1,8]

  frac1 <- 2 * pi * Tu * (Hu-Hl)
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  y <- frac1 / frac2
  y <- y / (pi * rw^2)
  return(y)
}

#' The Gradient of Borehole Function
#'
#' This function returns the gradient of the borehole function.
#'
#' @param xx the 7 inputs, restricted to the unit interval
#' @param design the radius of the borehole (typically rw)
#' @param adjust logical. adjustment for scaling needed?
#' @details PARAMETER RANGES
#' rw in `[0.05, 0.15]` 	radius of borehole (m)
#' r in `[100, 50000]` 	radius of influence (m)
#' Tu in `[63070, 115600]`    	transmissivity of upper aquifer (m2/yr)
#' Hu in `[990, 1110]` 	potentiometric head of upper aquifer (m)
#' Tl in `[63.1, 116]` 	transmissivity of lower aquifer (m2/yr)
#' Hl in `[700, 820]` 	potentiometric head of lower aquifer (m)
#' L in `[1120, 1680]` 	length of borehole (m)
#' Kw in `[9855, 12045]` 	hydraulic conductivity of borehole (m/yr)
#' @return The output of the borehole function
#' @export
borehole_grad <- function(xx, design=0.5, adjust=TRUE){
  rr<-matrix(c(
    0.01, 0.9,
    100, 50000,
    63070, 115600,
    990, 1110,
    63.1, 116,
    700, 820,
    1120, 1680,
    9855, 12045
  ),nrow=2)

  xx <- c(design, xx)
  rw <- xx[1]*diff(rr[,1]) + rr[1,1];r  <- xx[2]*diff(rr[,2]) + rr[1,2];Tu <- xx[3]*diff(rr[,3]) + rr[1,3];Hu <- xx[4]*diff(rr[,4]) + rr[1,4]
  Tl <- xx[5]*diff(rr[,5]) + rr[1,5];Hl <- xx[6]*diff(rr[,6]) + rr[1,6];L  <- xx[7]*diff(rr[,7]) + rr[1,7];Kw <- xx[8]*diff(rr[,8]) + rr[1,8]


  a <- r
  b <- Tu
  c <- Hu
  d <- Tl
  e <- Hl
  f <- L
  g <- Kw
  r <- rw

  #SWAP BECAUSE OF TYPO IN WOLFRAM
  tmp <- d
  d <- b
  b <- tmp

  grad <- rep(NA, 7)
  grad[1] <- -(2 *pi* b^2 *g^2 *r^4 *(b + d) *(c - e))/(a *(g* r^2 *(b + d)* log(a/r) + 2 *b *d *f)^2)
  grad[2] <-  (2 *pi* (c - e))/(log(a/r)* ((2* d *f)/(g *r^2 *log(a/r)) + d/b + 1)) + (2 *pi* d *(c - e))/(b *log(a/r) *((2 *d *f)/(g* r^2* log(a/r)) + d/b + 1)^2)
  grad[3] <- (2 *pi* b^2 *g *r^2)/(b* g* r^2 *log(a/r) + d* g *r^2 *log(a/r) + 2* b* d* f)
  grad[4] <- -(2 *pi* b^2* g *r^2 *(c - e)* (g *r^2 *log(a/r) + 2 *b *f))/(g* r^2 *(b + d)* log(a/r) + 2* b* d* f)^2
  grad[5] <- -(2 *pi* b^2* g *r^2)/(b* g* r^2* log(a/r) + d *g *r^2* log(a/r) + 2 *b* d* f)
  grad[6] <- (4 *pi* b^3 *d* g* r^2* (e - c))/(g* r^2* (b + d)* log(a/r) + 2* b* d* f)^2
  grad[7] <- (4 *pi* b^3 *d *f *r^2* (c - e))/(b *g* r^2* log(a/r) + d *g* r^2* log(a/r) + 2* b* d* f)^2

  #CORRECTION BECAUSE OF TYPO
  grad <- grad * d/b

  #REDO GRADIENT FOR T VARS
  grad[2] <- (2 *pi* b^2* g^2* r^4* (c - e) *log(a/r))/(b *g *r^2 *log(a/r) + d* g *r^2* log(a/r) + 2 *b* d* f)^2
  grad[4] <- (2 *pi* d^2 *g^2* r^4 *(c - e)* log(a/r))/(b *g *r^2 *log(a/r) + d *g *r^2* log(a/r) + 2 *b *d *f)^2

  #Adjust for scaling (if needed)
  if(adjust){
    grad <- grad * apply(rr, 2, diff)[2:8]
  }

  grad <- grad / (pi * rw^2)

  return(grad)
}


#' A Modified Borehole Function
#'
#' This function is for testing, it is a modified borehole function designed to have a more interesting active subspace
#'
#' @param xx 5 inputs, restricted to the unit interval. More inputs can be used but they are completely inert.
#' @param design the radius of the borehole (typically rw)
#' @details PARAMETER RANGES
#' rw in `[0.05, 0.15]` 	radius of borehole (m)
#' r in `[100, 50000]` 	radius of influence (m)
#' Tu in `[63070, 115600]`    	transmissivity of upper aquifer (m2/yr)
#' Hu in `[990, 1110]` 	potentiometric head of upper aquifer (m)
#' Tl in `[63.1, 116]` 	transmissivity of lower aquifer (m2/yr)
#' Hl in `[700, 820]` 	potentiometric head of lower aquifer (m)
#' L in `[1120, 1680]` 	length of borehole (m)
#' Kw in `[9855, 12045]` 	hydraulic conductivity of borehole (m/yr)
#' @return The output of the borehole function
modified_borehole <- function(xx, design=0.5){
  z1 <- (3*xx[1] - 2*xx[2] + 0.25*xx[3] + 2)/5.25
  z2 <- (4*xx[2] -5*xx[3] + 7*xx[4] + 5)/16
  borehole(c(z1, z2, xx), design=design)
}



#' Piston Function
#'
#' Piston function studied by Constantine in global sensitivity metrics paper
#'
#' @param x 7 inputs. See Constantine paper for details
#' @details PARAMETER RANGES:
#'measure <- function(){
#'res <- c(
#'  runif(1, 30, 60),
#'  runif(1, .005, .02),
#'  runif(1, .002, .01),
#'  runif(1, 1000, 5000),
#'  runif(1, 90000, 110000),
#'  runif(1, 290, 296),
#'  runif(1, 340, 360)
#')
#'return(res)
#'}
#'
#' @return Time to fire for piston
f_piston <- function(x){
  M  <- x[1]
  S  <- x[2]
  V0 <- x[3]
  k  <- x[4]
  P0 <- x[5]
  Ta <- x[6]
  T0 <- x[7]

  A <- P0*S + 19.62*M - k*V0/S
  V <- S/(2*k)*(sqrt(A^2+4*k*P0*V0*Ta/T0) - A)
  res <- 2*pi*sqrt(M/(k+S^2*P0*V0/T0*Ta/V^2))
  return(res)
}













