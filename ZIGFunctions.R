# Useful functions, including Zero Inflated Gamma (ZIG) density, 
# probability, quantile, and random generation functions.
# Also included is the negative log-likelihood.
# The gamma distribution functions are parametrized in terms of
#  mean and beta instead of the usual alpha and beta.

# expit and logit functions
alogit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

# ZIG functions
dsgamma <- function(y,mu,beta,epsilon) {
  (y>=epsilon)*dgamma(y-epsilon,shape=mu/beta,scale=beta)
}
psgamma <- function(y,mu,beta,epsilon) {
  (y>=epsilon)*pgamma(y-epsilon,shape=mu/beta,scale=beta)
}
rsgamma <- function(y,mu,beta,epsilon) {
  (y>=epsilon)*rgamma(1,shape=mu/beta,scale=beta)
}
qsgamma <- function(y,mu,beta,epsilon) {
  qgamma(y,shape=mu/beta,scale=beta) + epsilon
}
pzig <- function(y,mu,beta,epsilon,Pi) {
  Pi*punif(y,0,epsilon) + (1-Pi)*psgamma(y,mu,beta,epsilon)
}
dzig <- function(y,mu,beta,epsilon,Pi) {
  Pi*dunif(y,0,epsilon) + (1-Pi)*dsgamma(y,mu,beta,epsilon)
}

qzig <- function(u,mu,beta,epsilon,Pi) {
  ix <- which(u<=Pi)
  if (length(ix) == 0) {
      # no standard uniforms are less than Pi
      result <- qsgamma((u - Pi)/(1 - Pi), mu, beta, epsilon)
      return(result)
  } else {
      result <- numeric(length=length(u))
      result[ix] <- qunif(u[ix]/Pi[ix],0,epsilon) 
      result[-ix] <- qsgamma((u[-ix]-Pi[-ix])/(1-Pi[-ix]),mu[-ix],beta,epsilon)
      return(result)
      #ifelse(u<=Pi,qunif(u/Pi,0,epsilon),qsgamma((u-Pi)/(1-Pi),mu,beta,epsilon)) # gives warnings.
  }
}

# qzig <- function(u,mu,beta,epsilon,Pi) {
#   ix <- which(u<=Pi)
#   result <- numeric(length=length(u))
#   result[ix] <- qunif(u[ix]/Pi[ix],0,epsilon) 
#   result[-ix] <- qsgamma((u[-ix]-Pi[-ix])/(1-Pi[-ix]),mu[-ix],beta,epsilon)
#   return(result)
#   #ifelse(u<=Pi,qunif(u/Pi,0,epsilon),qsgamma((u-Pi)/(1-Pi),mu,beta,epsilon)) # gives warnings.
# }

rzig <- function(n,mu,beta,epsilon,Pi) {
  p <- rbinom(n,1,Pi)
  p*runif(n,0,epsilon) + (1-p)*rsgamma(n,shape=mu/beta,scale=beta,epsilon) 
}

# Negative Log Likelihood.
# theta is the vector of parameters, Y is the response vector, epsilon is the upper
# limit of the tiny uniforms, ZX is the design matrix for the Bernoulli part, 
# and NZX is the design matrix for the continuous part.
NL <- function(theta,Y,epsilon,H,X=NULL) {
  #print(theta)
  alphaN <- alogit(theta[1]) # spatial correlation "nugget"
  alphaR <- exp(theta[2])  # spatial correlation range
  lambda <- exp(theta[3])  # gamma *rate* parameter 1/scale
  beta <-  1/lambda # gamma scale parameter
  xi_Z <- matrix(theta[4:(3+ncol(X))],ncol=1) # Zero-inflation parameters
  Pi <- alogit(X%*%xi_Z) # probability of a zero.
  xi_NZ <- matrix(theta[(4+ncol(X)):(3+2*ncol(X))],ncol=1) # Gamma part parameters
  mu <- exp(X%*%xi_NZ) # gamma means
  alpha <- mu/beta # gamma shape parameters
  
  Id <- diag(length(Y))
  Sigma<-alphaN*exp(-H/alphaR)
  diag(Sigma)<-1
  logdetSigma <- determinant(Sigma)$modulus[1]
  Sig_inv <- solve(Sigma)
  
  FF <- pzig(Y,alpha,beta,epsilon,Pi)
  # Avoid numerical problems:
  FF[which(FF==1)] <- 1-.Machine$double.eps
  FF[which(FF==0)] <- .Machine$double.eps
  Z <- qnorm(FF)
  0.5*(logdetSigma
       + t(Z)%*%(Sig_inv-Id)%*%Z)
  -sum(log(dzig(Y,alpha,beta,epsilon,Pi)))
}
