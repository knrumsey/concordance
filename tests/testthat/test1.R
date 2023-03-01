test_that("simple polynomial example BASS-based AS matches truth", {
  cat('simple polynomial example test')

  f <- function(x){
    x[1]^2 + x[1]*x[2] + x[2]^3/9
  }

  # Sim Study Parameters
  N <- 1000
  p <- 3

  # Get true value of C matrix (analytically)
  Ctrue <- matrix(0, nrow=p, ncol=p)
  Ctrue[1:2, 1:2] <- matrix(1/45*c(120, 50, 50, 21), nrow=2, byrow=TRUE)

  X <- matrix(runif((N-2) * p), nrow=N-2, ncol=p)
  Yf <- apply(X, 1, f)

  mod_bass <- BASS::bass(X, Yf, verbose=FALSE)
  Cbass <- C_bass(mod_bass)

  d1 <- sum(abs(Cbass - Ctrue))
  expect_that(d1, is_less_than(0.04))
})
