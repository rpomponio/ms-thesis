################################################### -
## Title: Script containing functions for estimation
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-04
################################################### -

# function to compute maximally-conservative estimate of correlation
cor.conserv <- function(X, Y) {
  cor(sort(X), sort(Y, decreasing=TRUE))
}

# function to compute sample correlation of matched samples
cor.matched <- function(X, Y, n.matched) {
  if (n.matched==0) {
    return(NA)
  } else {
    return(cor(X[1:n.matched], Y[1:n.matched]))
  }
}

# function to compute bootstrap-based correlation estimates
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

# function to compute EM algorithm-based estimate of correlation
cor.emalg <- function(X, Y, n.matched, RHO.INIT=0, MAX.ITER=500) {
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

# function to compute Bayesian estimate of correlation
unif.cor.est <- function(S12,S1,S2,n)
{
  if(abs(S1+S2-2*S12)<(2e-6)*n || abs(S1+S2+2*S12)<(2e-6)*n)
  {
    print("CHEATED!!!! -- Unif")
    return(sign(S12)*1)
  }
  
  exp.rho.bayes <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (rho/2)*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  
  normalizing.const <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (1/2)*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  
  finite=FALSE
  iter = 0
  lim.up = 0.9
  while(finite==FALSE && iter < 8)
  {
    int.num1 <- integrate(exp.rho.bayes,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num2 <- integrate(exp.rho.bayes,lower=-lim.up,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num3 <- integrate(exp.rho.bayes,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    
    finite  <- (int.num1$message=="OK" & int.num2$message=="OK" & int.num3$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Unif num failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  }
  
  finite=FALSE
  iter = 0
  lim.up = 0.9
  
  while(finite==FALSE && iter < 8)
  {
    int.denom1 <- integrate(normalizing.const,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom2 <- integrate(normalizing.const,lower=-lim.up,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom3 <- integrate(normalizing.const,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)       
    
    finite  <- (int.denom1$message=="OK" & int.denom2$message=="OK" & int.denom3$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Unif denom failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  } 
  
  int.numerator <- int.num1$value + int.num2$value + int.num3$value
  int.denominator <- int.denom1$value + int.denom2$value + int.denom3$value
  
  rho.bayes <- int.numerator/int.denominator
  
  return(rho.bayes)
}

# wrapper
cor.bayesian.unif <- function(X, Y, n.matched) {
  SS.X <- sum(X[1:n.matched]^2)
  SS.Y <- sum(Y[1:n.matched]^2)
  SS.XY <- sum(X[1:n.matched] * Y[1:n.matched])
  
  if (n.matched==0) {
    return(NA)
  } else {
    return(unif.cor.est(SS.XY, SS.X, SS.Y, n=n.matched))
  }
}





