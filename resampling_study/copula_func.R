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
library(automap)


krige_zhat_gauscop <- function(training_data, test_data, resp, mu, beta, epsilon, Pi, zero.glm.coe, nonzero.glm.coe, s2) {
    
  #' Calculate predictions for unobserved locations using spatial gaussian copula and ordinary kriging
  #' 
  #' @param training_data dataframe of observed values, including annual precipitation
  #' @param test_data dataframe of unobserved values, including annual precipitation
  #' @param resp vector of transformed observed values
  #' @param mu vector of computed row-wise values
  #' @param beta scalar beta value for Gamma model
  #' @param epsilon smallest nonzero value
  #' @param Pi vector of proportions of zeros in responses
  #' @param zero.glm.coe coefficients for binomial model
  #' @param nonzero.glm.coe coefficients for gamma generalized linear model
  #' @param s2 estimated sample variance of observed responses
  #' 
  #' @return dataframe with columns of predicted standard normal values (zhat) and predicted ZIG values (zigs)
  
  z <- qnorm(pzig(y = resp,
                  mu = mu,
                  beta = beta,
                  epsilon = epsilon,
                  Pi = Pi))
  zhat <- ord_kriging_z(training_data = training_data, test_points = test_data, train_z = z)$krige_output$var1.pred
  
  pz <- pnorm(zhat)
  pz[which(pz==1)] <- 1-.Machine$double.eps
  pz[which(pz==0)] <- .Machine$double.eps
  
  X_test <- matrix(cbind(rep(1, nrow(test_data)), 
                         test_data$annpre), 
                   nrow = nrow(test_data)
                   )
  Pi_test <- alogit(X_test %*% zero.glm.coe)
  mu_test <- exp(X_test %*% nonzero.glm.coe)
  beta_test <- mean(s2/mu_test)
  
  zigs <- qzig(u = pz,
               mu = mu_test,
               beta = beta_test,
               epsilon = epsilon,
               Pi = Pi_test)
  
  output <- data.frame(zhat = zhat,
                       zigs = zigs)
  return(output)
}

calculate_spatial_params <- function(simdata, testdata) {
    
    #' Calculates estimates for the spatial covariance parameters using automap
    #' the spatial covariance parameters for regression kriging.
    #' For use with the simulated datasets - no covariate information used here
    #'
    #' @param simdata dataframe to calculate spatial parameters on
    #' @param test dataframe of values of unobserved locations to predict values on
    #' @param zeros boolean vector of zero values in response column
    #' 
    #' @return list of training data, testing data, correlation matrices required for prediction, correlogram parameters, and MOM estimators for gamma distribution
   
    H <- as.matrix(dist(simdata[,c("x", "y")]))
    if (class(simdata) == "data.frame") {
      coordinates(simdata) <- ~x+y
    }
    #Calculate covariance parameters using autofitVariogram
    suppressWarnings({
      v_fit <- autofitVariogram(resp~annpre,
                                input_data = simdata)
    })
    
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
      alphaR = alphaR
      # mu = mu,
      # beta = beta
      # epsilon = tmp$epsilon,
      # Pi = tmp$Pi,
    )
    
    return(output)
}
#' calculate_zhat_gauscop <- function(resp, mu, beta, epsilon, Pi, Sigma_obs, Sigma_obs_noobs) {
#'     
#'   #' Calculate predictions for unobserved locations using spatial gaussian copula and matrix algebra
#'   #' 
#'   #' @param resp vector of transformed observed values
#'   #' @param mu vector of computed row-wise values
#'   #' @param beta scalar beta value for Gamma model
#'   #' @param epsilon smallest nonzero value
#'   #' @param Pi vector of proportions of zeros in responses
#'   #' @param Sigma_obs spatial correlation matrix of observed values
#'   #' @param Sigma_obs_noobs spatial correlation matrix of observed values and unobserved values
#'   #' 
#'   #' @return dataframe with columns of predicted standard normal values (zhat) and predicted ZIG values (zigs)
#'   
#'   z <- qnorm(pzig(y = resp,
#'                   mu = mu,
#'                   beta = beta,
#'                   epsilon = epsilon,
#'                   Pi = Pi))
#'   
#'   S_obs_inv <- solve(Sigma_obs, tol = 1e-22)
#'   S_noobs_tr <- t(Sigma_obs_noobs)
#'   zhat <- S_noobs_tr %*% S_obs_inv %*% z
#'   
#'   pz <- pnorm(zhat)
#'   pz[which(pz==1)] <- 1-.Machine$double.eps
#'   pz[which(pz==0)] <- .Machine$double.eps
#'   
#'   zigs <- qzig(u = pz,
#'                mu = rep(mu[1], length(pz)),
#'                beta = beta,
#'                epsilon = epsilon,
#'                Pi = rep(Pi[1], length(pz)))
#'   
#'   output <- data.frame(zhat = zhat,
#'                        zigs = zigs)
#'   return(output)
#' }