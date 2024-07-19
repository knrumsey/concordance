# Load libraries
library(concordance)
library(BASS)
library(activegp)
library(hetGP)
library(lhs)
set.seed(1111)

# Define functions
f <- function(x){
  x[1]^2 + x[1]*x[2]
}
g <- function(x){
  x[1]^2 + x[1]*x[2] + x[2]^3/9
}

# True C matrix
Ctrue <- matrix(c(8/3, 10/9, 11/12, 7/18), nrow=2, byrow=2) #Worked this out on paper

# Monte Carlo (using true functions)
C0 <- Cfg_mc(f, g, measure=2)

# Create synthetic data
X <- maximinLHS(500, 3)
X <- rbind(X, rep(0, 3))
X <- rbind(X, rep(1, 3))
Yf <- apply(X, 1, f)
Yg <- apply(X, 1, g)

# Fit BASS models
mod1 <- bass(X, Yf)
mod2 <- bass(X, Yg)

# Estimate Cfg posterior (and posterior mean)
C1_list <- Cfg_bass(mod1, mod2, mcmc.use=seq(1, 1000, by=10))
C1 <- matrix(0, nrow=3, ncol=3)
for(i in 1:100){
  C1 <- C1 + C1_list[[i]]/100
}
C1

#Plot Cfg Matrices
par(mfrow=c(1,2))
image(Ctrue, main="True Cfg Matrix")
image(C1[1:2, 1:2], main="Estimated Cfg Matrix")

# Estimate Concordance
Cfg  <- Cfg_bass(mod1, mod2)
Cf   <- C_bass(mod1)
Cg   <- C_bass(mod2)
conc <- sum(diag(Cfg))/sqrt(sum(diag(Cf))*sum(diag(Cg)))
conc









