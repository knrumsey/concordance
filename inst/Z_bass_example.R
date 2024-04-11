library(duqling) #devtools::install_github("knrumsey/duqling")
library(BASS)
library(lhs)

X <- lhs::maximinLHS(1000, 3)
y <- apply(X, 1, duqling::detpep_curve)
mod <- bass(X, y)

Z = Z_bass(mod)
