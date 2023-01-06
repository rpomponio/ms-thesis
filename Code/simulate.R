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
REP <- 1:2

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

# iterate over all datasets in parallel, compute estimates
results <- foreach(i=1:1000, .combine=rbind) %dorng%{
  
  # generate a dataset given a parameter set
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
  boot_res <- cor.boot(X, Y, n_matched)
  
  # estimate correlation using various methods
  rho_conserv <- cor.conserv(X, Y)
  rho_matched <- cor.matched(X, Y, n_matched)
  rho_boot_mean <- boot_res$mean
  rho_boot_p5 <- boot_res$p5
  rho_boot_p20 <- boot_res$p20
  rho_emalg <- cor.emalg(X, Y, n_matched)
  
  # aggregate results
  list(
    "Distribution"=params$Distribution,
    "Rho"=params$Rho,
    "Delta"=params$Delta,
    "N"=params$N,
    "Prop.matched"=params$Prop.matched,
    "Sigma.X"=params$Sigma.X,
    "Sigma.Y"=params$Sigma.Y,
    "Repetition"=params$Repetition,
    "Max.conserv"=rho_conserv,
    "Matched"=rho_matched, "Boot.mean"=rho_boot_mean,
    "Boot.p5"=rho_boot_p5, "Boot.p20"=rho_boot_p20,
    "EM.alg"=rho_emalg)

}



