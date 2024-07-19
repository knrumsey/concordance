library(BASS)
library(lhs)

# Modified from duqling package. I also changed th priors from their default values
my_pollutant <- function(x, scale01=FALSE,
                      space=c(0.5, 1, 1.5, 2, 2.5),
                      time=seq(from=0.3, to=60, by=0.3)){
  if(scale01){
    RR <- cbind(c(5,   5, 0.02, 0.01, 0.003, 28),
                c(15, 15, 0.12,    3, 0.03,  32))
    x[1:6] <- x[1:6]*(RR[,2] - RR[,1]) + RR[,1]
  }


  M1  <- x[1]
  M2  <- x[2]
  D   <- x[3]
  L0  <- x[4]
  Lt  <- x[5]
  tau <- x[6]
  L   <- L0 + Lt*tau
  s   <- space
  t   <- time


  ds <- length(s)
  dt <- length(t)
  dY <- ds * dt
  Y <- matrix(0, ds, dt)

  # Create matrix Y, where each row corresponds to si and each column
  # corresponds to tj.
  for (ii in 1:ds) {
    si <- s[ii]
    for (jj in 1:dt) {
      tj <- t[jj]

      term1a <- M1 / sqrt(4*pi*D*tj)
      term1b <- exp(-si^2 / (4*D*tj))
      term1 <- term1a * term1b

      term2 <- 0
      if (tau < tj) {
        term2a <- M2 / sqrt(4*pi*D*(tj-tau))
        term2b <- exp(-(si-L)^2 / (4*D*(tj-tau)))
        term2 <- term2a * term2b
      }

      C <- term1 + term2
      Y[ii, jj] <- sqrt(4*pi) * C
    }
  }

  # Convert the matrix into a vector (by rows).
  Yrow <- t(Y)
  y <- t(as.vector(Yrow))
  return(y)
}

pollutant_modified <- function(X, xi, nfunc){
  my_pollutant(X, space=xi, time=seq(0.3, 60, length.out=nfunc), scale01=TRUE)
}

n <- 300
p <- 6
nfunc <- 30
X <- lhs::maximinLHS(300, p)
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
#np <- BASS:::bassPCAsetup(X, y1)$n.pc
fitPCA1 <- bassPCA(X, y1, n.pc=np, n.cores=3)
fitPCA2 <- bassPCA(X, y2, n.pc=np, n.cores=3)
#conc_method2 <- conc_bass(fitPCA1, fitPCA2)

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

#par(mfrow=c(1,2))
#matplot(y1, type='l', lty=1, col=adjustcolor("dodgerblue", alpha.f=0.25))
#matplot(y2, type='l', lty=1, col=adjustcolor("grey", alpha.f=0.25), add=TRUE)


plot(tt, conc_vec, type='l', pch=16, lwd=2, xlab="t", ylab="concordance",
     ylim=c(min(c(conc_vec,0)), 1.1), yaxt='n')
axis(2, seq(-1, 1, by=0.2), seq(-1, 1, by=0.2))
abline(h=conc_method1b$conc, lwd=2, lty=2, col='dodgerblue')
abline(h=conc1, lwd=2, lty=3, col='orange')
abline(h=conc2, lwd=2, lty=4, col='firebrick')
lines(tt, conc_vec, lwd=3)
legend("top", c(expression(paste(kappa["fg"],"(t),")),
                expression(paste(kappa["fg"]^"0",",")),
                expression(paste(kappa["fg"]^"1", ",")),
                expression(paste(kappa["fg"]^"2"))),
       lwd=2, lty=1:4, col=c("black", "dodgerblue", "orange", "firebrick"),
       horiz=TRUE, bty='n')


