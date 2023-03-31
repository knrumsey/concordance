library(BASS)
library(concordance)
library(parallel)
library(lhs)

tr <- function(xx) sum(diag(xx))
conc_C <- function(C1, C2, C12){
  tr(C12)/sqrt(tr(C1)*tr(C2))
}
# Full model should come first
concordance <- function(mod1, mod2){
  C1 <- C_bass(mod1)
  C2 <- C_bass(mod2)
  C12 <- Cfg_bass(mod1, mod2)
  conc_C(C1, C2, C12)
}


X <- matrix(lhs::randomLHS(10000, 2), ncol=2)
y <- apply(X, 1, duqling::dms_additive)
mod<-bass(X,y) # full model
nfolds<-20 # I’ll break the training data into 20 parts
scramble<-sample(nfolds,size=length(y),replace=T) # randomly allocate training data into folds
fit<-function(i){
  ind<-which(scramble<=i) # if i=5, train with the first 5 folds
  bass(X[ind,],y[ind])
}
mods<-parallel::mclapply(1:nfolds,fit,mc.cores = 5)
#qoi<-unlist(lapply(mods, function(m) mean(sqrt(m$s2/mod$s2)))) # for each model, get the relative error (relative to the full model)
# would rather have concordance as the qoi:
qoi <- unlist(lapply(mods, function(m) concordance(mod, m)))
plot(cumsum(table(scramble)),qoi,xlab='sample size')
