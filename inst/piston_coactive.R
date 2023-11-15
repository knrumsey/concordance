fname <- "piston"
p <- quack(fname)$input_dim
X <- lhs::maximinLHS(1000, p)
y <- duq(X, fname, scale01=TRUE)
mod <- bass(X, y)
C <- C_bass(mod)
plot(sqrt(abs(eigen(C)$values)), pch=16)
abline(h=0)

act_dims(C, X, y)


piston1 <- function(xx, ...){
  xx <- c(xx[1:4], 0, -1, xx[-(1:4)])
  duqling::piston(xx, scale01=TRUE, ...)
}

piston2 <- function(xx, ...){
  xx <- c(xx[1:4], 1, 2, xx[-(1:4)])
  duqling::piston(xx, scale01=TRUE, ...)
}

X <- lhs::maximinLHS(1000, 5)
y1 <- apply(X, 1, piston1)
y2 <- apply(X, 1, piston2)

# Analysis with concordance package
mod1 <- bass(X, y1)
mod2 <- bass(X, y2)

C1_post <- C_bass(mod1, mcmc.use=seq(1, 1000, by=5))
C2_post <- C_bass(mod2, mcmc.use=seq(1, 1000, by=5))
C12_post <- Cfg_bass(mod1, mod2, mcmc.use=seq(1, 1000, by=5))
C1 <- C2 <- C12 <- matrix(0, nrow=5, ncol=5)
for(i in 1:length(C1_post)){
  C1  <- C1 +  C1_post[[i]]/length(C1_post)
  C2  <- C2 +  C2_post[[i]]/length(C2_post)
  C12 <- C12 + C12_post[[i]]/length(C12_post)
}
V12 <- (C12+t(C12))/2

# Analysis with MC
B <- 100000
C1_mc <- C2_mc <- C12_mc <- matrix(0, nrow=5, ncol=5)
for(b in 1:B){
  xx <- runif(5)
  grad1 <- fd_grad(piston1, xx)
  grad2 <- fd_grad(piston2, xx)

  C1_mc  <- C1_mc +  tcrossprod(grad1, grad1)/B
  C2_mc  <- C2_mc +  tcrossprod(grad2, grad2)/B
  C12_mc <- C12_mc + tcrossprod(grad1, grad2)/B
}

png("figs/C12_mc.png", units="in", height=5, width=4, res=300)
image.plot(C12_mc, xaxt='n', yaxt='n', main="Monte Carlo Estimate", col=hcl.colors(50, "YlOrRd", rev=TRUE))
axis(1, seq(0, 1, length=5), 1:5)
axis(2, seq(0, 1, length=5), 1:5)
dev.off()

png("figs/C12_conc.png", units="in", height=5, width=4, res=300)
image.plot(C12, xaxt='n', yaxt='n', main="Concordance Estimate", col=hcl.colors(50, "YlOrRd", rev=TRUE))
axis(1, seq(0, 1, length=5), 1:5)
axis(2, seq(0, 1, length=5), 1:5)
dev.off()

png("figs/C12_comp.png", units="in", height=5, width=4, res=300)
plot(as.numeric(C12_mc), as.numeric(C12),
     pch=21, bg="dodgerblue")
abline(0, 1, lwd=2, lty=3)
dev.off()




