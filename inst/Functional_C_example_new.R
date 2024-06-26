library(BASS)
library(lhs)

pollutant_modified <- function(X, xi, nfunc){
  Z <- rep(NA, 4)
  Z[1] <- X[1]
  Z[2] <- X[2]
  Z[3] <- (X[3]*2 + X[4])/3
  Z[4] <- X[5]
  duqling::pollutant(Z, space=xi, time=seq(0.3, 60, length.out=nfunc), scale01=TRUE)
}

X <- lhs::maximinLHS(300, 5)
nfunc <- 20
y1 <- apply(X, 1, pollutant_modified, xi=1.5, nfunc=nfunc) #+ rnorm(200, 0, 1)
y2 <- apply(X, 1, pollutant_modified, xi=2.0, nfunc=nfunc) #+ rnorm(200, 0, 1)
tt <- seq(0, 1, length.out=nfunc)


# METHOD ONE (treat t as input) - Augmentation
fit1 <- bass(X, y1, xx.func=tt)
fit2 <- bass(X, y2, xx.func=tt)

conc_method1 <- conc_bass(fit1, fit2)

# method one, but surpress warnings.
fit1_modified <- bassfunc2bass(fit1)
fit2_modified <- bassfunc2bass(fit2)
conc_method1b <- conc_bass(fit1_modified, fit2_modified)

# Note that we can do other things with the modified models
# sob1 <- sobol(fit1_modified) # need to fix this, i think it has to do with knots instead of knot locations.

# METHOD TWO (Analyze as a function of t)
np <- BASS:::bassPCAsetup(X, y1)$n.pc
fitPCA1 <- bassPCA(X, y1, n.pc=np, n.cores=3)
fitPCA2 <- bassPCA(X, y2, n.pc=np, n.cores=3)
#conc_method2 <- conc_bass(fitPCA1, fitPCA2)


conc_vec <- rep(NA, nfunc)
K1 <- K2 <- K12 <- matrix(0, nrow=5, ncol=5)
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

#par(mfrow=c(1,2))
#matplot(y1, type='l', lty=1, col=adjustcolor("dodgerblue", alpha.f=0.25))
#matplot(y2, type='l', lty=1, col=adjustcolor("grey", alpha.f=0.25), add=TRUE)


plot(conc_vec, type='o', pch=16, lwd=2, xlab="t")
abline(h=conc_method1b$conc, lwd=2, lty=2, col='dodgerblue')
abline(h=conc1, lwd=2, lty=2, col='orange')
abline(h=conc2, lwd=2, lty=2, col='firebrick')



