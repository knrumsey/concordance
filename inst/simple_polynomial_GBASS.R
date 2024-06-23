library(GBASS)
library(BASS)
library(concordance)
library(lhs)
set.seed(1111)

# A simple function of two variables
f <- function(x){
  x[1]^2 + x[1]*x[2] + x[2]^3/9
}

n <- 500
n_corrupted <- 5
p <- 2
X <- lhs::randomLHS(n, p)
y <- apply(X, 1, f) + rnorm(n, 0, 0.01)
y[1:n_corrupted] <- y[1:n_corrupted] + rnorm(n_corrupted, 0.5)

#Fit BASS and GBASS Models
mod0 <- BASS::bass(X, y)
mod1 <- GBASS::tbass(X, y, df=5)
mod1 <- GBASS::gm2bm(mod1)

# Get C matrices
C0 <- C_bass(mod0)
C1 <- C_bass(mod1)
Ctrue <- matrix(1/45*c(120, 50, 50, 21), nrow=2, byrow=TRUE)

# Get error
D0 <- C0 - Ctrue
D1 <- C1 - Ctrue
sum(abs(D0))
sum(abs(D1))

#Print matrices
C0
C1
Ctrue
