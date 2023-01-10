################################################### -
## Title: Script containing functions for estimation
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-04
################################################### -

source(here::here("Code/helpers.R"))

#
# All estimators must take vectors X, Y as arguments
# Optionally, an estimator may take n.matched as an argument
#
# All estimators must return a scalar (or list of scalers) or NA values
#

# compute maximally-conservative estimate of correlation
est.cor.conserv <- function(X, Y) {
  cor(sort(X), sort(Y, decreasing=TRUE))
}

# compute sample correlation of matched samples
est.cor.matched <- function(X, Y, n.matched) {
  if (n.matched==0) {
    return(NA)
  } else {
    return(cor(X[1:n.matched], Y[1:n.matched]))
  }
}

# compute shrunken correlation of matched samples
est.cor.shrunken <- function(X, Y, n.matched) {
  
}

# compute "unbiased" correlation of matched samples
est.cor.unbiased <- function(X, Y, n.matched) {
  
}

# compute 20th quantile estimate of matched samples
est.cor.q.20 <- function(X, Y, n.matched) {
  
}


# compute bootstrap-based correlation estimate(s)
cor.boot <- function(X, Y, n.matched, N.BOOT=9999, SEED=123) {
  
  set.seed(SEED)
  
  if (n.matched==0) {
    return(list("mean"=NA, "p5"=NA, "p20"=NA, "p50"=NA))
  } else {
    rho.estimates <- rep(NA, N.BOOT)
    for (i in 1:N.BOOT){
      boot.indices <- sample(1:n.matched, replace=TRUE)
      X.boot <- X[boot.indices]
      Y.boot <- Y[boot.indices]
      if (is.na(sd(X.boot)) | sd(X.boot)==0) {
        rho.estimates[i] <- NA
      } else {
        rho.estimates[i] <- cor(X.boot, Y.boot)        
      }
    }
    
    # @Ryan: how to handle cases where correlation can't be estimated?
    if (mean(is.na(rho.estimates)) > 0.1){
      rho.p5 <- NA
      rho.p20 <- NA
      rho.p50 <- NA
      rho.mean <- NA
    } else {
      rho.p5 <- quantile(rho.estimates, 0.05, na.rm=TRUE)
      rho.p20 <- quantile(rho.estimates, 0.20, na.rm=TRUE)
      rho.p50 <- quantile(rho.estimates, 0.50, na.rm=TRUE)
      rho.mean <- mean(rho.estimates, na.rm=TRUE)
    }
    
    return(list("mean"=rho.mean, "p5"=unname(rho.p5),
                "p20"=unname(rho.p20), "p50"=unname(rho.p50)))
  }
}

# compute EM algorithm-based estimate of correlation
est.cor.emalg <- function(X, Y, n.matched, RHO.INIT=0, MAX.ITER=500) {
  n <- length(X)
  n.unmatched <- n - n.matched
  matched.indices <- 1:n.matched
  
  X.unmatched <- X[-matched.indices]
  Y.unmatched <- Y[-matched.indices]
  
  if (n.matched > 0) {
    X.matched <- X[matched.indices]
    Y.matched <- Y[matched.indices]
    crossprod.obs <- crossprod(X.matched, Y.matched)[1, 1]
  } else {
    # @Ryan: is it correct to assume the cross product is zero?
    crossprod.obs <- 0
  }
  
  # compute MLES of mu, sigma
  mu.X <- sum(X) / n
  mu.Y <- sum(Y) / n
  sigma.X <- sqrt( sum((X - mu.X) ^2) / n )
  sigma.Y <- sqrt( sum((Y - mu.Y) ^2) / n )
  
  # MLEs of mu and sigma found on unmatched samples, only
  mu.X.unmatched <- sum(X.unmatched) / n.unmatched
  mu.Y.unmatched <- sum(Y.unmatched) / n.unmatched
  sigma.X.unmatched <- sqrt( sum((X.unmatched - mu.X.unmatched) ^2) / n.unmatched )
  sigma.Y.unmatched <- sqrt( sum((Y.unmatched - mu.Y.unmatched) ^2) / n.unmatched )
  
  # initialize starting value(s) of rho
  rho.p <- -2
  rho.p2 <- RHO.INIT
  
  # define convergence criterion as difference in estimates of rho
  n.iter <- 0
  eps <- 1e-5
  while(abs(rho.p - rho.p2) >= eps & n.iter <= MAX.ITER) {
    n.iter <- n.iter + 1
    rho.p <- rho.p2
    # e-step (unobserved cross product of X, Y)
    crossprod.unobs <- sum(
      rho.p * sigma.X.unmatched * sigma.Y.unmatched,
      mu.X.unmatched * mu.Y.unmatched) * n.unmatched
    # crossprod augmented includes both observed and unobserved cases
    crossprod.aug <- crossprod.obs + crossprod.unobs
    # m-step (conditional MLE of rho given partially observed data)
    rho.p2 <- (crossprod.aug - n * mu.X * mu.Y) / (n * sigma.X * sigma.Y)
    if(rho.p2 > 1 | rho.p2 < (-1)){
      warning(paste0("Estimate of rho=", rho.p2, " in iteration:", n.iter))
    }
    
  }
  
  if(n.iter >= MAX.ITER){
    warning("Reached maximum number of iterations.")
  }
  return(rho.p2)
}

# wrapper for computing bayesian estimate with uniform prior
est.cor.bayesian.unif <- function(X, Y, n.matched) {
  SS.X <- sum(X[1:n.matched]^2)
  SS.Y <- sum(Y[1:n.matched]^2)
  SS.XY <- sum(X[1:n.matched] * Y[1:n.matched])
  
  if (n.matched==0) {
    return(NA)
  } else {
    return(unif.cor.est(SS.XY, SS.X, SS.Y, n=n.matched))
  }
}





