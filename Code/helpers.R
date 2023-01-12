################################################### -
## Title: Functions for Bayesian Estimators
## Author: Bailey Fosdick
################################################### -


########## BAYESIAN METHODS ############

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


jeff.cor.est <- function(S12,S1,S2,n)
{
  
  if(abs(S1+S2-2*S12)<(2e-6)*n || abs(S1+S2+2*S12)<(2e-6)*n)
  {
    print("CHEATED!!!! -- Jeff")
    return(sign(S12)*1)
  }    
  
  exp.rho.jeff <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (rho*sqrt(1+rho^2))/(1-rho^2)*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  
  normalizing.const.jeff <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (sqrt(1+rho^2))/(1-rho^2)*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  
  finite=FALSE
  iter = 0
  lim.up = 0.9
  while(finite==FALSE && iter < 8)
  {
    int.num1 <- integrate(exp.rho.jeff,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num2 <- integrate(exp.rho.jeff,lower=-lim.up,upper=0,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num3 <- integrate(exp.rho.jeff,lower=0,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num4 <- integrate(exp.rho.jeff,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    
    finite  <- (int.num1$message=="OK" & int.num2$message=="OK" & int.num3$message=="OK" & int.num4$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Jeff num failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  } 
  
  
  finite=FALSE
  iter = 0
  lim.up = 0.9
  
  while(finite==FALSE && iter < 8)
  {
    int.denom1 <- integrate(normalizing.const.jeff,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom2 <- integrate(normalizing.const.jeff,lower=-lim.up,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom3 <- integrate(normalizing.const.jeff,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)       
    
    finite  <- (int.denom1$message=="OK" & int.denom2$message=="OK" & int.denom3$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Jeff denom failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  } 
  
  int.numerator <- int.num1$value + int.num2$value + int.num3$value + int.num4$value
  int.denominator <- int.denom1$value + int.denom2$value + int.denom3$value
  
  rho.jeff <- int.numerator/int.denominator
  
  return(rho.jeff)
}

arcsine.cor.est <- function(S12,S1,S2,n)
{
  if(abs(S1+S2-2*S12)<(2e-6)*n || abs(S1+S2+2*S12)<(2e-6)*n)
  {
    print("CHEATED!!!! -- Arcsine")
    return(sign(S12)*1)
  }
  
  exp.rho.arcsine <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (rho/pi)*(1/sqrt(1-rho^2))*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  
  normalizing.const.arcsine <- function(rho)
  {
    exp.term <- -(1/(2*(1-rho^2)))*(S1-2*rho*S12+S2)
    result <-  (1/pi)*(1/sqrt(1-rho^2))*((2*pi*sqrt(1-rho^2))^(-n))*exp(exp.term)
    return(result)
  }
  finite=FALSE
  iter = 0
  lim.up = 0.9
  while(finite==FALSE && iter < 8)
  {
    int.num1 <- integrate(exp.rho.arcsine,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num2 <- integrate(exp.rho.arcsine,lower=-lim.up,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.num3 <- integrate(exp.rho.arcsine,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    
    finite  <- (int.num1$message=="OK" & int.num2$message=="OK" & int.num3$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Arcsine num failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  } 
  
  finite=FALSE
  iter = 0
  lim.up = 0.9
  
  while(finite==FALSE && iter < 8)
  {
    int.denom1 <- integrate(normalizing.const.arcsine,lower=-1,upper=-lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom2 <- integrate(normalizing.const.arcsine,lower=-lim.up,upper=lim.up,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)
    int.denom3 <- integrate(normalizing.const.arcsine,lower=lim.up,upper=1,subdivisions=50000,abs.tol=10e-50,stop.on.error=FALSE)       
    
    finite  <- (int.denom1$message=="OK" & int.denom2$message=="OK" & int.denom3$message=="OK")
    lim.up <- lim.up*1.1
    iter <- iter+1 
  }
  
  if(finite==FALSE)
  {
    res <- min(abs(S1+S2-2*S12),abs(S1+S2+2*S12))
    print(paste("Arcsine denom failed: diff= ",res))
    print(paste("S1=",S1,", S2=",S2,", S12=",S12))
    return(sign(S12)*1)
  } 
  
  int.numerator <- int.num1$value + int.num2$value + int.num3$value
  int.denominator <- int.denom1$value + int.denom2$value + int.denom3$value
  
  rho.ArcSine <- int.numerator/int.denominator
  
  return(rho.ArcSine)
}