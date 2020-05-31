# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for calculating spatial parameters for use in the spatial Gaussian copula and using them for prediction

source("../ZIGFunctions.R")
source("kriging_func.R")
source("preprocess.R")
library(gstat)
library(sp)
library(fields)

calculate_spatial_params <- function(simdata, testdata, zeros) {
    
    #' Calculates Method of Moments estimates for Gamma model and the spatial covariance parameters for ordinary kriging
    #' For use with the simulated datasets - no covariate information used here
    #'
    #' @param simdata dataframe to calculate spatial parameters on
    #' @param test dataframe of values of unobserved locations to predict values on
    #' @param zeros boolean vector of zero values in response column
    #' 
    #' @return list of training data, testing data, correlation matrices required for prediction, correlogram parameters, and MOM estimators for gamma distribution
    
    H <- as.matrix(dist(simdata[,c("x", "y")]))
    
    # Use method of moment estimator for gamma distribution 
    if (length(zeros) > 0) {
      mu <- mean(simdata$resp[-zeros])
      beta <- var(simdata$resp[-zeros])/mu
      mu <- rep(mu, times=nrow(simdata))
    } else {
      mu <- mean(simdata$resp)
      beta <- var(simdata$resp)/mu
      mu <- rep(mu, times=nrow(simdata))
    }
    
    simdata.cpy <- simdata
    coordinates(simdata.cpy) <- ~x+y
    #Calculate covariance parameters using autofitVariogram
    v_fit <- autofitVariogram(resp~1,
                              input_data = simdata.cpy)
    
    # v_fit <- tryCatch({
    #     autofitVariogram(resp~1,
    #                      input_data = simdata)
    # }, warning = function(war) {
    #     warning("Variogram did not converge.")
    #     stop("Variogram did not converge.")
    # }, message = function(mes) {
    #     message("Hmm something weird.")
    # })
    
    alphaN<-v_fit$var_model$psill[2]/sum(v_fit$var_model$psill) # nugget parameter
    alphaR<-v_fit$var_model$range[2] # decay parameter
   
    #Calculate spatial correlation matrix w. exponential form for observed points!
    Sigma_obs<-alphaN*exp(-(H/alphaR)^2)
    diag(Sigma_obs)<-1
    
    #Calculate spatial correlation matrix for observed/unobserved points! 
    H_obs_noobs <- as.matrix(rdist(simdata[,c("x", "y")], testdata[,c("x", "y")]))
    Sigma_obs_noobs <- alphaN*exp(-(H_obs_noobs/alphaR)^2)
    
    output <- list(
      Sigma_obs = Sigma_obs,
      Sigma_obs_noobs = Sigma_obs_noobs,
      alphaN = alphaN,
      alphaR = alphaR,
      mu = mu,
      beta = beta
      # epsilon = tmp$epsilon,
      # Pi = tmp$Pi,
    )
    
    return(output)
}

calculate_zhat_gauscop <- function(resp, mu, beta, epsilon, Pi, Sigma_obs, Sigma_obs_noobs) {
    
  #' Calculate predictions for unobserved locations using spatial gaussian copula and matrix algebra
  #' 
  #' @param resp vector of transformed observed values
  #' @param mu vector of computed row-wise values
  #' @param beta scalar beta value for Gamma model
  #' @param epsilon smallest nonzero value
  #' @param Pi vector of proportions of zeros in responses
  #' @param Sigma_obs spatial correlation matrix of observed values
  #' @param Sigma_obs_noobs spatial correlation matrix of observed values and unobserved values
  #' 
  #' @return dataframe with columns of predicted standard normal values (zhat) and predicted ZIG values (zigs)
  
  z <- qnorm(pzig(y = resp,
                  mu = mu,
                  beta = beta,
                  epsilon = epsilon,
                  Pi = Pi))
  
  S_obs_inv <- solve(Sigma_obs, tol = 1e-22)
  S_noobs_tr <- t(Sigma_obs_noobs)
  zhat <- S_noobs_tr %*% S_obs_inv %*% z
  
  pz <- pnorm(zhat)
  pz[which(pz==1)] <- 1-.Machine$double.eps
  pz[which(pz==0)] <- .Machine$double.eps
  
  zigs <- qzig(u = pz,
               mu = rep(mu[1], length(pz)),
               beta = beta,
               epsilon = epsilon,
               Pi = rep(Pi[1], length(pz)))
  
  output <- data.frame(zhat = zhat,
                       zigs = zigs)
  return(output)
}

krige_zhat_gauscop <- function(training_data, test_data, resp, mu, beta, epsilon, Pi) {
    
  #' Calculate predictions for unobserved locations using spatial gaussian copula and ordinary kriging
  #' 
  #' @param resp vector of transformed observed values
  #' @param mu vector of computed row-wise values
  #' @param beta scalar beta value for Gamma model
  #' @param epsilon smallest nonzero value
  #' @param Pi vector of proportions of zeros in responses
  #' @param Sigma_obs spatial correlation matrix of observed values
  #' @param Sigma_obs_noobs spatial correlation matrix of observed values and unobserved values
  #' 
  #' @return dataframe with columns of predicted standard normal values (zhat) and predicted ZIG values (zigs)
  
  z <- qnorm(pzig(y = resp,
                  mu = mu,
                  beta = beta,
                  epsilon = epsilon,
                  Pi = Pi))
  if (sum(duplicated(z)) == (length(z)-1)) {
    zhat <- rep(z[1], nrow(test_data))
  } else {
    zhat <- ord_kriging_z(training_data = training_data, test_points = test_data, test_z = z)$krige_output$var1.pred
  }
  
  pz <- pnorm(zhat)
  pz[which(pz==1)] <- 1-.Machine$double.eps
  pz[which(pz==0)] <- .Machine$double.eps
  
  zigs <- qzig(u = pz,
               mu = rep(mu[1], length(pz)),
               beta = beta,
               epsilon = epsilon,
               Pi = rep(Pi[1], length(pz)))
  
  output <- data.frame(zhat = zhat,
                       zigs = zigs)
  return(output)
}
