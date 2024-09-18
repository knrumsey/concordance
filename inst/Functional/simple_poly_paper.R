library(concordance)
library(BASS)
library(lhs)

# Define functions
f <- function(x, nT=30){
  res <- (x[1]^2 + x[1]*x[2])*rep(1, nT)
  return(res)
}

g <- function(x, nT=30){
  tt <- seq(0, 1, length=nT+2)[-c(1, nT+2)]
  res <- (x[1]^2 + x[1]*x[2] + 10*(2*tt-1)*x[2]^3)*tt^2*(1-tt)
  return(res)
}

# Define true values
Cf_true <- function(t){
  res <- matrix(c(32,11,11,4)/12,
                nrow=2, byrow=TRUE)
  return(res)
}
Cg_true <- function(t){
  u <- 2*t-1
  a <- t^4*(1-t)^2/12
  B <- matrix(c(32, 11+210*u, 11+210*u, 4 + 120*u + 2160*u^2),
              nrow=2, byrow=TRUE)
  return(a*B)
}
Cfg_true <- function(t){
  u <- 2*t-1
  a <- t^2*(1-t)/12
  B <- matrix(c(32, 11+210*u, 11, 4 + 60*u), nrow=2, byrow=TRUE)
  return(a*B)
}
conc_true <- function(t){
  u <- 2*t-1
  (3+5*u)/sqrt(9+30*u + 540*u^2)
}
Kf_true <- function(){
  res <- matrix(c(32, 11, 11, 4), nrow=2, byrow=TRUE)/12
  return(res)
}
Kg_true <- function(){
  res <- matrix(c(64, 127, 127, 788), nrow=2, byrow=TRUE)/2520
  return(res)
}
Kfg_true <- function(){
  res <- matrix(c(32, 53, 11, 16), nrow=2, byrow=TRUE)/144
  return(res)
}

# Generate training data
p <- 2
n <- 800
nT <- 25
tt <- seq(0, 1, length=nT+2)[-c(1, nT+2)]
X <- lhs::randomLHS(n, p)
y1 <- apply(X, 1, f, nT=nT)
y2 <- apply(X, 1, g, nT=nT)

# Fit BASS models
np1 <- BASS:::bassPCAsetup(X, y1)$n.pc
fit1 <- bassPCA(X, y1, n.pc=np1, n.cores=min(np1, 3))

np2 <- BASS:::bassPCAsetup(X, y2)$n.pc
fit2 <- bassPCA(X, y2, n.pc=np2, n.cores=min(np2, 3))

# Estimate C models
tic()
C1 <- C_bassPCA(fit1)
toc()
tic()
C2 <- C_bassPCA(fit2)
toc()
tic()
C12 <- Cfg_bassPCA(fit1, fit2)
toc()

cat("Alternate version")
tic()
C1_alt <- C_bassPCA_v2(fit1, prior=pr)
toc()
tic()
C2_alt <- C_bassPCA_v2(fit2, prior=pr)
toc()
tic()
C12_alt <- Cfg_bassPCA_v2(fit1, fit2, prior=pr)
toc()

# Calculate concordance
tr1 <- unlist(lapply(C1_alt, function(A) sum(diag(A))))
tr2 <- unlist(lapply(C2_alt, function(A) sum(diag(A))))
tr12 <- unlist(lapply(C12_alt, function(A) sum(diag(A))))
conc_est <- tr12/sqrt(tr1*tr2)
curve(conc_true(x), lwd=2)
points(tt, conc_est, pch=16, col='dodgerblue')

# Check errors
for(t in 1:nT){
  print(sum(abs(Cf_true(tt[t]) - C1_alt[[t]]))/sum(abs(Cf_true(tt[t]))))
  print(sum(abs(Cg_true(tt[t]) - C2_alt[[t]]))/sum(abs(Cg_true(tt[t]))))
  print(sum(abs(Cfg_true(tt[t]) - C12_alt[[t]]))/sum(abs(Cfg_true(tt[t]))))
}

# Get K matrices
Kf <- Kg <- Kfg <- matrix(0, nrow=2, ncol=2)
for(t in 1:nT){
  Kf <- Kf + C1[[t]]/nT
  Kg <- Kg + C2[[t]]/nT
  Kfg <- Kfg + C12[[t]]/nT
}

print(sum(abs(Kf - Kf_true())))
print(sum(abs(Kg - Kg_true())))
print(sum(abs(Kfg - Kfg_true())))

# Try augmented case
y1n <- t(apply(y1, 1, function(row) row + rnorm(length(row), mean = 0, sd = 0.01 * sd(row))))
y2n <- t(apply(y2, 1, function(row) row + rnorm(length(row), mean = 0, sd = 0.01 * sd(row))))

fit1_aug <- bass(X, y1n, xx.func=tt)
fit2_aug <- bass(X, y2n, xx.func=tt)

fit1_aug_mod <- bassfunc2bass(fit1_aug)
fit2_aug_mod <- bassfunc2bass(fit2_aug)
conc_aug <- conc_bass(fit1_aug_mod, fit2_aug_mod)

C1_aug <- C_bass(fit1_aug_mod, prior=pr)
C2_aug <- C_bass(fit2_aug_mod, prior=pr)
C12_aug <- Cfg_bass(fit1_aug_mod, fit2_aug_mod, prior=pr)
tr(C12_aug)/sqrt(tr(C1_aug)*tr(C2_aug))


library(RColorBrewer)
png("figs/simple_poly_conc.png", width=8, height=5, units="in", res=300)
par(mfrow=c(1,1))
curve(conc_true(x), lwd=2, lty=3, xlab="t", ylab="Concordance")
abline(h=c(0.1459, 0.331, 0.347), col=brewer.pal(5, "Set2"), lwd=1)
points(c(0,0,0), c(0.1459, 0.331, 0.347), pch=c(1,2,5), cex=1.4, lwd=2, col=brewer.pal(5, "Set2"))
points(tt, conc_est, pch=3, lwd=1, col=brewer.pal(5, "Set2")[4])
legend('topleft', c("Concordance(t)", "Estimated values", "Augmentation", "Direct aggregation", "JIT aggregation"),
       lty=c(3,-1,1,1,1), lwd=c(2,1,1,1,1), col=c("black", brewer.pal(5, "Set2")[c(4,1:3)]), pch=c(-1,3,1,2,5))
dev.off()
