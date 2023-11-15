test_that("simple polynomial example with Gaussian measure", {
  cat('simple polynomial example test with Gaussian measure')

  f <- function(x){
    x[1]^2 + x[1]*x[2] + x[2]^3/9
  }

  # Sim Study Parameters
  N <- 5000
  p <- 3
  sd0 <- 0.05

  # Get true value of C matrix (with monte carlo)
  measure <- function() rnorm(p, 0.5, sd0)
  Cmc <- C_mc(f, measure, nmc=1e5)

  X <- matrix(rnorm(N*p, 0.5, sd0), nrow=N, ncol=p)
  X <- lhs::randomLHS(N, p)
  Yf <- apply(X, 1, f)

  mod_bass <- BASS::bass(X, Yf, verbose=FALSE)
  pr <- list()
  pr[[1]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)
  pr[[2]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)
  pr[[3]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)

  Cbass <- C_bass(mod_bass, prior=pr, scale01=TRUE)

  d1 <- sum(abs(Cbass - Cmc))
  expect_that(d1, is_less_than(0.04))
})
