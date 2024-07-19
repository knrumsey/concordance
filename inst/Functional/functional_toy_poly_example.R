library(BASS)
library(lhs)

f1 <- function(x, tt){
  res <- rep(x[1]^2 + x[1]*x[2], length(tt))
}

f2 <- function(x, tt, a=0, b=15){
  x[1]^2 + x[1]*x[2] + (tt*(b-a) + a)*x[2]^3
}

# Generate toy data
#aa = -20
#bb = 5
aa = 0
bb = 15
n <- 300
p <- 3
X <- lhs::maximinLHS(n, p)
nfunc <- 20
tt <- seq(0, 1, length.out=nfunc)
y1 <- apply(X, 1, f1, tt=tt) + rnorm(n, 0, 0.001)
y2 <- apply(X, 1, f2, tt=tt) + rnorm(n, 0, 0.001)

# Augmentation
fit1 <- bass(X, y1, xx.func=tt)
fit2 <- bass(X, y2, xx.func=tt)
conc_method1 <- conc_bass(fit1, fit2)

fit1_modified <- bassfunc2bass(fit1)
fit2_modified <- bassfunc2bass(fit2)
conc_method1b <- conc_bass(fit1_modified, fit2_modified)

# As a function of t
np <- max(3, BASS:::bassPCAsetup(X, y1)$n.pc)
fitPCA1 <- bassPCA(X, y1, n.pc=np, n.cores=3)
fitPCA2 <- bassPCA(X, y2, n.pc=np, n.cores=3)

conc_vec <- rep(NA, nfunc)
K1 <- K2 <- K12 <- matrix(0, nrow=p, ncol=p)
for(t in 1:nfunc){
  print(t)
  mod_list <- fitPCA1$mod.list
  weights <- fitPCA1$dat$basis[t,]
  fit1_lc <- lcbass2bass(mod_list, weights)
  mod_list <- fitPCA2$mod.list
  weights <- fitPCA2$dat$basis[t,]
  fit2_lc <- lcbass2bass(mod_list, weights)

  C1 <- C_bass(fit1_lc)
  C2 <- C_bass(fit2_lc)
  C12 <- Cfg_bass(fit1_lc, fit2_lc)
  conc_vec[t] <- tr(C12)/sqrt(tr(C1)*tr(C2))

  K1 <- K1 + C1/nfunc
  K2 <- K2 + C2/nfunc
  K12 <- K12 + C12/nfunc
}
conc1 <- tr(K12)/sqrt(tr(K1)*tr(K2))
conc2 <- mean(conc_vec)

# Make figure
plot(tt, conc_vec, type='l', pch=16, lwd=2, xlab="t", ylab="concordance(t)", ylim=c(-0.2, 1))
abline(h=conc_method1b$conc, lwd=2, lty=2, col='dodgerblue')
abline(h=conc1, lwd=2, lty=3, col='orange')
abline(h=conc2, lwd=2, lty=4, col='firebrick')
lines(tt, conc_vec, lwd=3)

conc_fg_func <- function(t, a, b){
  beta <- (b-a)*t + a
  term1 <- 3 + beta/2
  term2 <- 9 + 3*beta*(9/5*beta + 1)
  term1/sqrt(term2)
}
curve(conc_fg_func(x, a=aa, b=bb), add=TRUE, lty=3, col="red")
