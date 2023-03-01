library(concordance)
library(BASS)
library(activegp)
library(hetGP)
library(lhs)
set.seed(1111)

# A simple function of two variables
f <- function(x){
  x[1]^2 + x[1]*x[2] + x[2]^3/9
}

# Get truth with monte carlo
measure <- function(){
  flag <- TRUE
  while(flag){
    xx = runif(2)
    if(xx[1] > xx[2])
      flag = FALSE
  }
  return(xx)
}
C0 = C_mc(f, measure)


# Generate Data
N <- 500
X <- maximinLHS(N-2, 2)
X <- rbind(X, rep(0, 2))
X <- rbind(X, rep(1, 2))
Yf <- apply(X, 1, f)

#Fit BASS model
mod <- bass(X, Yf)

#Construct prior list
L <- 3 # Number of mixture components
rho <- list()
for(i in 1:L){
  tmp <- list()
  tmp[[1]] <- list(dist="uniform", trunc=c(i/(L+1), 1))
  tmp[[2]] <- list(dist="uniform", trunc=c((i-1), i)/(L+1))
  rho[[i]] <- tmp
}


C_ell <- list()
for(i in 1:L){
  A <- diag(c((1 - i/(L+1)), 1/(L+1)))
  C_ell[[i]] <- A%*%C_bass(mod, prior=rho[[i]])[[1]]%*%t(A)
}

C1 <- matrix(0, nrow=2, ncol=2)
for(i in 1:L){
  C1 <- C1 + C_ell[[i]]/L
}
C1





