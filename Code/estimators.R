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
# All estimators must return a scalar (or list of scalars), or NA value(s)
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

# compute "shrunken" correlation of matched samples
est.cor.shrunken <- function(X, Y, n.matched) {
  if (n.matched==0) {
    return(NA)
  } else {
    r <- cor(X[1:n.matched], Y[1:n.matched])
    return(sign(r) * sqrt(abs((1 - (1-r^2) * ((n.matched - 1)/(n.matched - 2))))))
  }
}

# compute "unbiased" correlation of matched samples
est.cor.unbiased <- function(X, Y, n.matched) {
  if (n.matched==0) {
    return(NA)
  } else {
    r <- cor(X[1:n.matched], Y[1:n.matched])
    return(r * (1  + (1 - r^2) / (2 * (n.matched - 3))))
  }
}

# compute 20th quantile estimate of matched samples (requires 3 samples)
est.cor.quantile <- function(X, Y, n.matched, q=0.2) {
  if (n.matched < 3) {
    return(NA)
  } else {
    test_res <- cor.test(X[1:n.matched], Y[1:n.matched],
                         alternative="greater", conf.level=(1 - q))
    return(test_res$conf.int[1])
  }
}

# compute bootstrap-based correlation estimate(s)
cor.boot <- function(X, Y, n.matched, N.BOOT=9999, SEED=123) {
  
  set.seed(SEED)
  
  if (n.matched==0) {
    return(list("mean"=NA, "q.05"=NA, "q.20"=NA, "q.50"=NA))
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
      rho.q.05 <- NA
      rho.q.20 <- NA
      rho.q.50 <- NA
      rho.mean <- NA
      warning("Ten or more percent of bootstrap iterations failed.")
    } else {
      rho.q.05 <- quantile(rho.estimates, 0.05, na.rm=TRUE)
      rho.q.20 <- quantile(rho.estimates, 0.20, na.rm=TRUE)
      rho.q.50 <- quantile(rho.estimates, 0.50, na.rm=TRUE)
      rho.mean <- mean(rho.estimates, na.rm=TRUE)
    }
    
    return(list("mean"=rho.mean, "q.05"=unname(rho.q.05),
                "q.20"=unname(rho.q.20), "q.50"=unname(rho.q.50)))
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
  if (n.matched < n) {
    mu.X.unmatched <- sum(X.unmatched) / n.unmatched
    mu.Y.unmatched <- sum(Y.unmatched) / n.unmatched
    sigma.X.unmatched <- sqrt( sum((X.unmatched - mu.X.unmatched) ^2) / n.unmatched )
    sigma.Y.unmatched <- sqrt( sum((Y.unmatched - mu.Y.unmatched) ^2) / n.unmatched )
  }

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
    if (n.matched < n) {
      crossprod.unobs <- sum(
        rho.p * sigma.X.unmatched * sigma.Y.unmatched,
        mu.X.unmatched * mu.Y.unmatched) * n.unmatched
    } else {
      crossprod.unobs <- 0
    }
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





