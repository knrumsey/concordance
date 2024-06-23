pollutant_modified <- function(X, xi, nfunc){
  Z <- rep(NA, 4)
  Z[1] <- X[1]
  Z[2] <- X[2]
  Z[3] <- (X[3]*2 + X[4])/3
  Z[4] <- X[5]
  duqling::pollutant(Z, space=xi, time=seq(0.3, 60, length.out=nfunc), scale01=TRUE)
}

X <- lhs::maximinLHS(200, 5)
y1 <- apply(X, 1, pollutant_modified, xi=1.5, nfunc=10) + rnorm(200, 0, 0.5)
y2 <- apply(X, 1, pollutant_modified, xi=2.0, nfunc=10) + rnorm(200, 0, 0.5)
tt <- seq(0, 1, length.out=10)


# METHOD ONE (treat t as input)
fit1 <- bass(X, y1, xx.func=tt)
fit2 <- bass(X, y2, xx.func=tt)
conc_method1 <- conc_bass(fit1, fit2)

# METHOD TWO (Analyze as a function of t)
fitPCA1 <- bassPCA(X, y1)
fitPCA2 <- bassPCA(X, y2)
conc_method1 <- conc_bass(fitPCA1, fitPCA2)


