library(concordance)
library(BASS)
set.seed(1000809)

t0 <- Sys.time()
f1 <- function(x){
  x[1]^2 + x[1]*x[2]
}

f3 <- function(x){
  f1(x) + x[3]
}

X <- matrix(runif(300), ncol=3)
y1 <- apply(X, 1, f1)
y3 <- apply(X, 1, f3)

# FIT BASS MODELS
mod1 <- bass(X, y1)
mod3 <- bass(X, y3)

# CONCORDANCE
C1 <- C_bass(mod1)
C3 <- C_bass(mod3)
C13 <- Cfg_bass(mod1, mod3)
tr <- function(A) sum(diag(A))
tr(C13)/sqrt(tr(C1)*tr(C3))

tf <- Sys.time()


