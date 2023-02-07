# Simulating the "effective" correlation for Ordinal datasets

library(mvtnorm)
library(here)

set.seed(2023)

N.ITER <- 1000

rho_init <- seq(-0.9, 0.9, 0.1)
eff_rho_avg <- rep(NA, length(rho_init))

for (rho in rho_init) {

  rho_temp <- rep(NA, N.ITER)
  
  params <- list(Sigma.X=1, Sigma.Y=1, Rho=rho, N=1000, Delta=0,
                 Prop.matched=1.0, Distribution="Ordinal")
  
  for (i in 1:N.ITER) {
    S <- diag(2)
    diag(S) <- c(params$Sigma.X, params$Sigma.Y)
    S[upper.tri(S) | lower.tri(S)] <- rep(params$Rho) * params$Sigma.X * params$Sigma.Y
    sim_mat <- rmvnorm(
      n=params$N,
      mean=c(0, params$Delta),
      sigma=S,
      checkSymmetry=TRUE)
    X <- sim_mat[, 1]
    Y <- sim_mat[, 2]

    n_matched <- floor(params$Prop.matched * params$N)

    # if ordinal: cut to integers from 1 to 7
    if (params$Distribution=="Ordinal") {
      cut_seq <- c(-999, c(-2, 1, 0, 0.7, 1.3, 2.0), 999)
      X <- as.numeric(cut(X, cut_seq))
      Y <- as.numeric(cut(Y, cut_seq))
    }

    rho_temp[i] <- cor(X, Y)
  }

  eff_rho_avg[which(rho_init==rho)] <- mean(rho_temp)
}


plot(rho_init, eff_rho_avg, type='l', xlim=c(-1, 1), ylim=c(-1, 1))
lines(rho_init, rho_init, col='blue')
title('Initial values of rho versus "effective correlation"')

