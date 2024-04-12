#library(duqling) #devtools::install_github("knrumsey/duqling")
library(BASS)
library(lhs)

# TRUE ANSWER IS (1.5, beta+0.5, 0)
f <- function(xx, beta=1/3) xx[1]^2 + xx[1]*xx[2] + beta*xx[2]^3

# WITH CONCORDANCE
X <- lhs::maximinLHS(100, 3)
y <- apply(X, 1, f)
mod <- bass(X, y)
Z = Z_bass(mod)
Z

# WITH MONTE CARLO
Zmc <- rep(0, 3)
for(i in 1:nrow(X)){
  Zmc <- Zmc + fd_grad(f, X[i,])/nrow(X)
}
Zmc
