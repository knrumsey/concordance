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

# Sim Study Parameters
N <- 100
p <- 5

X <- lhs::maximinLHS(N, p)
y <- apply(X, 1, f)
mod <- bass(X, y)

# Set Gaussian prior
pr <- build_prior(rep("normal", 5), mean=0.5, sd=0.1)
C_gaus <- C_bass(mod, prior=pr)

