################################################### -
## Title: Script to run simulation (in parallel)
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-04
################################################### -

library(doParallel)
library(doRNG)
library(mvtnorm)
library(here)

source(here("Code/estimators.R"))

# HARDCODED PARAMETERS
DISTRIB <- c("Normal")
RHO <- c(-0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9)
DELTA <- c(0) ###c(0, 0.25, 0.5)
N <- c(10) ###c(10, 20, 50, 100, 200)
PROP.MATCHED <- c(0.4, 0.8, 1) ###c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
SIGMA.X <- c(1)
SIGMA.Y <- c(1)
REP <- 1:10000

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
time_start <- proc.time()
results <- foreach(i=1:n_datasets, .combine=rbind) %dorng%{
  
  # generate a dataset given a parameter set
  params <- df_grid[i, ]
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
  
  # setup
  n_matched <- floor(params$Prop.matched * params$N)
  # boot_res <- cor.boot(X, Y, n_matched)
  
  # estimate correlation using various methods and aggregate results
  list(
    "Repetition"=params$Repetition,
    "Distribution"=as.numeric(params$Distribution),
    "Rho"=params$Rho,
    "Delta"=params$Delta,
    "Sigma.X"=params$Sigma.X,
    "Sigma.Y"=params$Sigma.Y,
    "N"=params$N,
    "M"=n_matched,
    "Prop.matched"=params$Prop.matched,
    # include all estimators below this line:
    "Max.conserv"        = est.cor.conserv(X, Y),
    "Matched"            = est.cor.matched(X, Y, n_matched),
    "Shrunken"           = est.cor.shrunken(X, Y, n_matched),
    "Unbiased"           = est.cor.unbiased(X, Y, n_matched),
    "Freq.20th.quantile" = est.cor.quantile(X, Y, n_matched, q=0.2),
    # "Boot.mean"          = boot_res$mean,
    # "Boot.5th.quantile"  = boot_res$q.05,
    # "Boot.20th.quantile" = boot_res$q.20,
    "EM.alg"             = est.cor.emalg(X, Y, n_matched),
    "Bayes.unif"         = est.cor.bayesian.unif(X, Y, n_matched)
    )
}
time_elapsed <- proc.time() - time_start

# save to file
sys_date <- gsub(" \\d+:\\d+:\\d+", "", Sys.time())
saveRDS(
  data.matrix(as.data.frame(results, row.names=F)),
  here("DataRaw", paste0("simulation_results_", sys_date, ".rds"))
  )

# print to signal completion
cat("Completed in", time_elapsed[3], "seconds.")
