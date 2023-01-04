################################################### -
## Title: Script to run simulation
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-04
################################################### -

library(doParallel)
library(doRNG)
library(mvtnorm)

# HARDCODED PARAMETERS
DISTRIB <- c("Normal")
RHO <- c(-0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9)
DELTA <- c(0, 0.25, 0.5)
N <- c(10, 20, 50, 100, 200)
PROP.MATCHED <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
SIGMA.X <- c(1)
SIGMA.Y <- c(1)
REP <- 1:1000

# register parallel backend
cl <- detectCores() - 2
registerDoParallel(cl)

# ensures reproducibility with `%dorng%`
set.seed(2023)

# expand grid of parameter choices
df_grid <- expand.grid(DISTRIB, RHO, DELTA, N, PROP.MATCHED, SIGMA.X, SIGMA.Y, REP)
colnames(df_grid) <- c("Distribution", "Rho", "Delta", "N", "Prop.matched",
                       "Sigma.X", "Sigma.Y", "Repetition")
n_datasets <- nrow(df_grid)
cat("Simulating", n_datasets, "datasets using", cl, "cores...")

# iterate over all datasets in parallel
res <- foreach(i=1:n_datasets, .combine=rbind) %dorng%{
  
  # generate a dataset given a choice of parameters
  params <- df_grid[i, ]
  S <- diag(2)
  diag(S) <- c(params$Sigma.X, params$Sigma.Y)
  S[upper.tri(S) | lower.tri(S)] <- rep(params$Rho) * prod(params$Sigma.X, params$Sigma.Y)
  sim_mat <- rmvnorm(
    n=params$N,
    mean=c(0, params$Delta),
    sigma=S,
    checkSymmetry=FALSE)
  X <- sim_mat[, 1]
  Y <- sim_mat[, 2]
  
  # setup
  n_matched <- floor(params$Prop.matched * params$N)
  
  # estimate correlation using various methods
  cor(X[1:n_matched], Y[1:n_matched])
  
  # aggregate results
  c(mean(X) - mean(Y))
}



