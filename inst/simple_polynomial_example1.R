library(concordance)
library(BASS)
library(activegp)
library(hetGP)
library(lhs)
set.seed(1111)

# A simple function of two variables
f <- function(x){
  x[1]^2 + x[1]*x[2] + x[2]^3/9
}

# Sim Study Parameters
N <- 100
p_vec <- c(2, 3, 4, 5, 10, 20, 50, 100)
pp <- max(p_vec)

# Get true value of C matrix (analytically)
Ctrue <- matrix(0, nrow=pp, ncol=pp)
Ctrue[1:2, 1:2] <- matrix(c(8/3, 10/9, 10/9, 21/45), nrow=2, byrow=TRUE)
Ctrue[1:2, 1:2] <- matrix(1/45*c(120, 50, 50, 21), nrow=2, byrow=TRUE)

# Allocate storage matrices
metrics1 <- matrix(0, nrow=length(p_vec), ncol=6)
metrics2 <- matrix(0, nrow=length(p_vec), ncol=6)

t00 <- Sys.time()
REPS <- 1
for(rep in 1:REPS){
  # Generate Data
  X <- maximinLHS(N-2, pp)
  X <- rbind(X, rep(0, pp))
  X <- rbind(X, rep(1, pp))
  Yf <- apply(X, 1, f)

  for(i in 1:length(p_vec)){
    p <- p_vec[i]

    # activegp method
    t0 <- Sys.time()
    mod_gp <- mleHomGP(X[,1:p], Yf)
    tgp1 <- Sys.time() - t0
    Cgp = C_GP(mod_gp)
    tgp2 <- Sys.time() - t0

    # BASS method
    t0 <- Sys.time()
    mod_bass <- bass(X[,1:p], Yf, verbose=FALSE, temp.ladder=1.1^(0:5))
    tbass1 <- Sys.time() - t0
    Cbass = C_bass(mod_bass)
    tbass2 <- Sys.time() - t0

    # Save time metrics
    metrics1[i,1] <- metrics1[i,1] + tgp1/REPS
    metrics2[i,1] <- metrics2[i,1] + tbass1/REPS

    metrics1[i,2] <- metrics1[i,2] + tgp2/REPS
    metrics2[i,2] <- metrics2[i,2] + tbass2/REPS

    # Save matrix norm metrics
    D1 <- Cgp[[1]] - Ctrue[1:p,1:p]
    D2 <- Cbass - Ctrue[1:p,1:p]

    #L1
    metrics1[i,3] <- metrics1[i,3] + sum(abs(D1))/REPS
    metrics2[i,3] <- metrics2[i,3] + sum(abs(D2))/REPS

    #L2
    metrics1[i,4] <- metrics1[i,4] + sqrt(sum(D1^2))/REPS
    metrics2[i,4] <- metrics2[i,4] + sqrt(sum(D2^2))/REPS

    #Linf
    metrics1[i,5] <- metrics1[i,5] + max(abs(D1))/REPS
    metrics2[i,5] <- metrics2[i,5] + max(abs(D2))/REPS

    #Frobenius
    metrics1[i,6] <- metrics1[i,6] + sqrt(sum(diag(D1%*%t(D1))))/REPS
    metrics2[i,6] <- metrics2[i,6] + sqrt(sum(diag(D2%*%t(D2))))/REPS

    print(p)
  }
  print("REP")
  print(rep)
}
print(Sys.time() - t00)

# Make plot
library(RColorBrewer)
bob <- brewer.pal(7, "Set2")
par(mfrow=c(1,2))




# Time plot
plot(p_vec, t1, ylim=c(0, max(rbind(metrics1[,1:2], metrics2[,1:2]))),
     xlab="p", ylab="time (s)", col=bob[1], lwd=2, lty=1, type='o', pch=16, main='Runtime')
lines(p_vec, t1, lwd=2, lty=2, col=bob[1])
points(p_vec, t1, pch=16, col=bob[1])

lines(p_vec, t2, lwd=2, lty=1, col=bob[2])
points(p_vec, t2, pch=15, col=bob[2])
#lines(p_vec, metrics2[,1], lwd=2, lty=2, col=bob[2])
#points(p_vec, metrics2[,1], pch=15, col=bob[2])

# Error Plot
k <- 6
ind <- 1:8
plot(p_vec[ind], metrics1[ind,k], ylim=c(range(rbind(metrics1[ind,k], metrics2[ind,k]))),
     xlab="p", ylab="Error", col=bob[1], lwd=2, lty=1, type='o', pch=16, main='Error')
lines(p_vec[ind], metrics2[ind,k], lwd=2, lty=1, col=bob[2])
points(p_vec[ind], metrics2[ind,k], pch=15, col=bob[2])


colnames(metrics1) <- c("model time", "full time", "L1", "L2", "Linf", "Lfrob")
rownames(metrics1) <- p_vec
colnames(metrics2) <- c("model time", "full time", "L1", "L2", "Linf", "Lfrob")
rownames(metrics2) <- p_vec

save(metrics1, metrics2, file=paste0("metrics_N", N, "_", cntr, ".rda"))
cntr <- cntr + 1



