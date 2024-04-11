library(tictoc)
library(lhs)
library(BASS)
f1 <- function(x) x[1]^2 + x[1]*x[2]
f2 <- function(x, beta=1/9) x[1]^2 + x[1]*x[2] + beta*x[2]^3

X <- lhs::maximinLHS(5000, 2)
y1 <- apply(X, 1, f1)
y2 <- apply(X, 1, f2, beta=3)

tic()
m1 <- bass(X, y1)
t1 = toc()
tic()
m2 <- bass(X, y2)
t2 = toc()

tic()
C1  <- C_bass(m1)
C2  <- C_bass(m2)
C12 <- Cfg_bass(m1, m2)
t3 = toc()

sum(diag(C12))/sqrt(sum(diag(C1)*sum(diag(C2))))


# Analytic stuff
conc <- function(beta) (3+beta/2)/sqrt(9 + 3*beta*(1.8*beta + 1))
C1_f <- function(beta){
  matrix(c(480, 165, 165, 60), nrow=2, byrow=TRUE)/180
}
C2_f <- function(beta){
  matrix(c(480, 165 + 315*beta, 165 + 315*beta, 60 + beta*(324*beta + 180)), nrow=2, byrow=TRUE)/180
}
C12_f <- function(beta){
  matrix(c(480, 165 + 315*beta, 165, 60 + beta*90), nrow=2, byrow=TRUE)/180
}
V12_f <- function(beta){
  matrix(c(480, 165 + 315*beta/2, 165 + 315*beta/2, 60 + beta*90), nrow=2, byrow=TRUE)/180
}
coactivity <- function(beta,q=1){
  C12 <- C12_f(beta)
  V12 <- (C12+ t(C12))/2
  Eig <- eigen(V12)
  a <- rep(NA, nrow(Eig$vectors))
  for(i in 1:nrow(Eig$vectors)){
    a[i] <- sum(Eig$values[1:q]*Eig$vectors[i,1:q]^2)
  }
  return(a)
}

coactivity2 <- function(beta,q=1){
  C12 <- C12_f(beta)
  V12 <- (C12+ t(C12))/2
  Eig <- eigen(V12)
  a <- rep(NA, nrow(Eig$vectors))
  for(i in 1:nrow(Eig$vectors)){
    a[i] <- sum(abs(Eig$values[1:q])*Eig$vectors[i,1:q]^2)
  }
  return(a)
}

# 1. Signed Activity Score
tmpf1 <- function(x) coactivity(x,1)[1]
tmpf2 <- function(x) coactivity(x,1)[2]
tmpf3 <- function(x) coactivity(x,2)[1]
tmpf4 <- function(x) coactivity(x,2)[2]

beta_vec <- seq(-15, 5, by=0.25)
zz1 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf1(beta_vec[zz]))))
zz2 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf2(beta_vec[zz]))))
zz3 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf3(beta_vec[zz]))))
zz4 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf4(beta_vec[zz]))))

#png(filename="figs/simple_ex2.png", width=5, height=5, units="in", res=300)
bob <- RColorBrewer::brewer.pal(2, "Set1")
plot(NULL, xlim=c(-15, 5), ylim=range(c(zz1, zz2, zz3, zz4)),
     xlab=expression(paste(beta)), ylab="Signed co-activity score", bty='n')
points(beta_vec, zz1, col=bob[1], pch=15, cex=0.5, lwd=1)
points(beta_vec, zz3, col=bob[1], pch=0, cex=0.5, lwd=1)
points(beta_vec, zz2, col=bob[2], pch=16, cex=0.5, lwd=1)
points(beta_vec, zz4, col=bob[2], pch=1, cex=0.5, lwd=1)
legend("bottomright",
       c(expression(paste(alpha[1],", q=1")),
         expression(paste(alpha[1],", q=2")),
         expression(paste(alpha[2],", q=1")),
         expression(paste(alpha[2],", q=2"))),
       bty='n', col=bob[c(1,1,2,2)], pch=c(15,0,16,1), cex=1)
#dev.off()

# 2. Unsigned Activity Score
tmpf1 <- function(x) coactivity2(x,1)[1]
tmpf2 <- function(x) coactivity2(x,1)[2]
tmpf3 <- function(x) coactivity2(x,2)[1]
tmpf4 <- function(x) coactivity2(x,2)[2]

beta_vec <- seq(-15, 5, by=0.25)
zz1 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf1(beta_vec[zz]))))
zz2 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf2(beta_vec[zz]))))
zz3 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf3(beta_vec[zz]))))
zz4 <- Re(unlist(lapply(seq_along(beta_vec), function(zz) tmpf4(beta_vec[zz]))))

#png(filename="figs/simple_ex3.png", width=5, height=5, units="in", res=300)
bob <- RColorBrewer::brewer.pal(5, "Set1")
plot(NULL, xlim=c(-15, 5), ylim=range(c(zz1, zz2, zz3, zz4)),
     xlab=expression(paste(beta)), ylab="Unsigned co-activity score", bty='n')
points(beta_vec, zz1, col=bob[1], pch=15, cex=0.7, lwd=1)
points(beta_vec, zz3, col=bob[1], pch=0, cex=0.7, lwd=1)
points(beta_vec, zz2, col=bob[2], pch=16, cex=0.7, lwd=1)
points(beta_vec, zz4, col=bob[2], pch=1, cex=0.7, lwd=1)
legend("top",
       c(expression(paste(tilde(alpha)[1],", q=1")),
         expression(paste(tilde(alpha)[1],", q=2")),
         expression(paste(tilde(alpha)[2],", q=1")),
         expression(paste(tilde(alpha)[2],", q=2"))),
       bty='n', col=bob[c(1,1,2,2)], pch=c(15, 0, 16,1), cex=1)
##dev.off()

# 3. Concordance Plot
#png(filename="figs/simple_ex1.png", width=5, height=5, units="in", res=300)
curve(conc(x), from=-15, to=5, xlab=expression(paste(beta)),ylim=c(-1, 1), bty='n',
      lwd=2, col=bob[3], ylab=expression(paste("conc(",f[1], ", ", f[2], ")")))
abline(h=c(-1,0,1), lwd=2, lty=3)
#dev.off()

# 4. Relative Contribution to Concordance

contrib <- function(beta){
  C1 = C1_f(beta)
  C2 = C2_f(beta)
  C12 = C12_f(beta)
  V12 = (C12 + t(C12))/2
  t1 <- tr(C1)
  t2 <- tr(C2)
  t12 <- tr(V12)
  p1 <- eigen(C1)$values/t1
  p2 <- eigen(C2)$values/t2
  p12a <- eigen(V12)$values/t12
  p12b <- eigen(V12)$values/sqrt(t1*t2)
  return(p12b)
}

#png(filename="figs/simple_ex1b.png", width=5, height=5, units="in", res=300)
zz0 <- matrix(unlist(lapply(seq_along(beta_vec), function(zz) contrib(beta_vec[zz]))), ncol=2, byrow=TRUE)
plot(NULL, xlim=c(-15, 5), ylim=range(c(zz0[,1], zz0[,2])),
     xlab=expression(paste(beta)), ylab="Contribution to concordance", bty='n')
lines(beta_vec, zz0[,1], lty=1, lwd=2, col=bob[4])
lines(beta_vec, zz0[,2], lty=2, lwd=2, col=bob[5])
legend("topleft",
       c(expression(paste(pi[1])),
         expression(paste(pi[2]))),
       bty='n', col=bob[4:5], lwd=2, lty=1:2, cex=1)
abline(h=0, lty=3, lwd=2)
#dev.off()


