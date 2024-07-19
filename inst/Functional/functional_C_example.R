n  <- 200
p  <- 4
nt <- 50
X <- matrix(runif(n*p), ncol=p)
y <- apply(X, 1, duqling::pollutant, space=2,
           time = seq(0.3, 60, length.out=nt), scale01=TRUE)
t <- seq(0, 1, length.out=nt)
fit_func <- BASS::bass(X, y, xx.func=t)
fit <- bassfunc2bass(fit_func)

#sob <- sobol(fit) # Commented out, but verified that this works.
C <- C_bass(fit)
plot(eigen(C)$values)
act_scores(C, k=2, plt=TRUE, norm=TRUE)
