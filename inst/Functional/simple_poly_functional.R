library(BASS)
library(lhs)

f1 <- function(x, tt){
  res <- rep(x[1]^2 + x[1]*x[2], length(tt))
}

f2 <- function(x, tt, a=-10, b=10){
  x[1]^2 + x[1]*x[2] + (tt*(b-a) + a)*x[2]^3
}

grlee1 <- function(x, scale01=FALSE){
  if(scale01){
    x <- 0.5 + 2*x
  }
  term1 <- sin(10*pi*x) / (2*x)
  term2 <- (x-1)^4

  y <- term1 + term2
  return(y)
}

f3 <- function(x, tt, a, b, phi){
  #f1(x, tt)*(1+plogis(tt, 0.5, 0.05))
  f2(x, tt, a ,b)*phi(tt)
}

phi1 <- function(tt) 1
phi2 <- function(tt) dnorm(tt, 0.5, 1/8)^2 + 0.1

n <- 300
p <- 2
X <- lhs::maximinLHS(n, p)
nfunc <- 20
tt <- seq(0, 1, length.out=nfunc)
y1 <- apply(X, 1, f1, tt=tt) + rnorm(n, 0, 0.001)
y2 <- apply(X, 1, f3, tt=tt, a=-10, b=10, phi1) + rnorm(n, 0, 0.001)

# Augmentation
fit1 <- bass(X, y1, xx.func=tt)
fit2 <- bass(X, y2, xx.func=tt)
conc_method1 <- conc_bass(fit1, fit2)

# same thing but slightly more rigorous
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


plot(tt, conc_vec, type='l', pch=16, lwd=2, xlab="t", ylab="concordance(t)",
     ylim=c(-1, 1), main="phi1")
abline(h=conc_method1b$conc, lwd=2, lty=2, col='dodgerblue')
abline(h=conc1, lwd=2, lty=3, col='orange')
abline(h=conc2, lwd=2, lty=4, col='firebrick')
lines(tt, conc_vec, lwd=3)



#### SECOND VERSION
n <- 300
p <- 2
X <- lhs::maximinLHS(n, p)
nfunc <- 20
tt <- seq(0, 1, length.out=nfunc)
y1 <- apply(X, 1, f1, tt=tt) + rnorm(n, 0, 0.001)
y2 <- apply(X, 1, f3, tt=tt, a=-10, b=10, phi2) + rnorm(n, 0, 0.001)

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


# COMPARE RESULTS WITH ANALYTIC RESULTS

par(mfrow=c(1,2))
plot(tt, conc_vec, type='l', pch=16, lwd=2, xlab="t", ylab="concordance(t)",
     ylim=c(-.25, 1), main="phi2")
abline(h=conc_method1b$conc, lwd=2, lty=2, col='dodgerblue')
abline(h=conc1, lwd=2, lty=3, col='orange')
abline(h=conc2, lwd=2, lty=4, col='firebrick')
lines(tt, conc_vec, lwd=3)



a <- 50
mu <- 0.5
sig <- 1/6
pow <- 3
s = function(tt) 0.1+ a*dnorm(tt, mu, sig)^pow
sp = function(tt) dnorm_pow_prime(tt, mu, sig, pow, a)
m = function(tt) 1

a <- 50
s = function(tt) 0.1 +  gamma(2*a)/gamma(a)^2*(tt)^a*(1-tt)^a
sp = function(tt) gamma(2*a)/gamma(a)^2*a*tt^(a-1)*(1-tt)^(a-1)*(1-2*tt)


a = 5000
eps = 0.15
s = function(tt) 0 + a * as.numeric(tt < 0.5 + eps ) * as.numeric(tt > 0.5 - eps)
sp = function(tt) 0

curve(conc_func_poly(x), ylim=c(-0.25, 1), xlab="t", ylab="concordance(t)", main="Analytic Results")
abline(h=conc_jit_poly, col="firebrick", lty=4, lwd=1)
abline(h=conc_agg_poly(s, sp, m), col="orange", lty=4, lwd=1)
abline(h=conc_aug_poly(s, sp, m), col="dodgerblue", lty=4, lwd=1)


