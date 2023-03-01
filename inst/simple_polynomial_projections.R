library(concordance)
library(BASS)
library(activegp)
library(hetGP)
library(lhs)
library(rgl)
set.seed(1111)s

# A simple function of two variables
f <- function(x){
  x[1]^2 + x[1]*x[2] + x[2]^3/9
}

# Simulation Parameters
set.seed(1111)
N <- 500
p <- 6
X <- maximinLHS(N, p)
A <- matrix(c(4, 3, 2, 1/3, 1/3, 1/3,
              1/3, 1/3, 1/3, 5, 2, 2), nrow=p, ncol=2)
Z <- X%*%A
y <- apply(Z, 1, f)

# Fit BMARS model and get active subspace
mod <- bass(X, y)
C <- C_bass(mod)
W <- eigen(C)$vectors

mycol <- function(yy, nn){
  yy <- (yy-min(yy))/diff(range(yy))
  1 + floor(yy*(nn-1))
}


# Make Figure
#' options(rgl.printRglwidget = TRUE)
#' png("../figures/poly_2dplot.png", height=5, width=5, units="in", res=300)
#' bob <- RColorBrewer::brewer.pal(11, "RdBu")
#' plot(X%*%W[,1], y, xlab="First active direction", ylab="Model output",
#'      pch=21, bg=bob[mycol(y, length(bob))], col="black")
#' dev.off()

#' png("../figures/poly_3dplot.png", height=5, width=5, units="in", res=300)
#' plot3d(X%*%W[,1], X%*%W[,2], y,
#'        xlab="First active direction",
#'        ylab="Second active direction",
#'        zlab="Model output",
#'        col=bob[mycol(y, length(bob))], size=0.8, type=c("s"))
#' dev.off()

rmspe <- rep(NA, 6)
for(i in 1:6){
  tmp <- bass(X%*%W[,1:i], y)
  rmspe[i] <- sqrt(mean((tmp$yhat.mean - y)^2))
}

tab <- matrix(NA, nrow=2, ncol=6)
rownames(tab) <- c("Square Root Eigenvalue", "Pct of SD Explained")
tab[1,] <- sqrt_eval/max(sqrt_eval)
tab[2,] <- (1 - rmspe/sd(y))*100
stargazer::stargazer(tab)






