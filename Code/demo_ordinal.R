# Simulating the "effective" correlation for Ordinal datasets

library(mvtnorm)
library(here)

set.seed(2023)

params <- list(Sigma.X=1, Sigma.Y=2, Rho=0.5, N=10000, Delta=0,
               Prop.matched=1.0, Distribution="Ordinal")

S <- diag(2)
diag(S) <- c(params$Sigma.X, params$Sigma.Y)
S[upper.tri(S) | lower.tri(S)] <- rep(params$Rho) * params$Sigma.X * params$Sigma.Y
sim_mat <- rmvnorm(
  n=params$N,
  mean=c(0, params$Delta),
  sigma=S,
  checkSymmetry=FALSE)
X <- sim_mat[, 1]
Y <- sim_mat[, 2]

n_matched <- floor(params$Prop.matched * params$N)

# if ordinal: cut to integers from 1 to 7
if (params$Distribution=="Ordinal") {
  cut_seq <- c(-999, c(-2, 1, 0, 0.7, 1.3, 2.0), 999)
  X <- as.numeric(cut(X, cut_seq))
  Y <- as.numeric(cut(Y, cut_seq))
}

cor(X, Y)

