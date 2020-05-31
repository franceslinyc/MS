# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Degenerate case where Gaussian copula produces all 0s
# This is just a sandbox to figure out how to handle it

setwd("~/MS_Project/simulation_study/")
simdata <- readRDS("../simulated_Y.rds")
dat <- read.csv("../ForestDataFuzzed.csv")

source("../ZIGFunctions.R")
source("preprocess.R")
source("copula_func.R")
source("kriging_func.R")
source("rfsp_func.R")
source("metrics.R")

library(gstat)
library(tidyverse)
library(rgdal)
library(sp)
library(automap)

alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
                  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
colnames(alb.xy) <- c("x","y")
dat <- cbind(alb.xy,dat)

full_simdata <- data.frame(
    x = dat$x,
    y = dat$y,
    resp = simdata[404,]
)
    
# Create spatial Gaussian copula model
set.seed(182)
tmp <- split_data(full_simdata, n = 300)
hist(tmp$train$resp)
hist(tmp$test$resp)
theta <- NULL

transformed.training.set <- cube_root_fix_zeros(tmp$train)
while(is.null(theta)) {
    try(
        theta <- calculate_spatial_params(simdata = transformed.training.set$transformed_data,
                                          testdata = tmp$test,
                                          zeros = transformed.training.set$zeros)
    )
}
krige_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                   mu = theta$mu,
                   beta = theta$beta,
                   epsilon = transformed.training.set$epsilon,
                   Pi = transformed.training.set$Pi,
                   training_data = tmp$train,
                   test_data = tmp$test) -> zhat.zigs.df.krige

calculate_zhat_gauscop(resp = transformed.training.set$transformed_data$resp,
                       mu = theta$mu,
                       beta = theta$beta,
                       epsilon = transformed.training.set$epsilon,
                       Pi = transformed.training.set$Pi,
                       Sigma_obs = theta$Sigma_obs,
                       Sigma_obs_noobs = theta$Sigma_obs_noobs) -> zhat.zigs.df.alg

yhat_gauscop.krige <- backtransform_gauscop(zhat.zigs.df.krige$zigs, transformed.training.set$epsilon)
yhat_gauscop.alg <- backtransform_gauscop(zhat.zigs.df.alg$zigs, transformed.training.set$epsilon)

hist(yhat_gauscop.krige)
hist(yhat_gauscop.alg)

transformed.training.set$epsilon
transformed.training.set$Pi

z <- qnorm(pzig(y = transformed.training.set$transformed_data$resp,
              mu = theta$mu,
              beta = theta$beta,
              epsilon = transformed.training.set$epsilon,
              Pi = transformed.training.set$Pi))
hist(z)

# # Below seems alright -- issue seems to lie with qzig() function
# 
# mu = theta$mu
# mu
# theta$beta
# theta$alphaN
# theta$alphaR
# 
# training_data <- tmp$train 
# coordinates(training_data) <- ~x+y 
# emp_var<-variogram(resp~1,loc=~x+y,training_data)
# plot(emp_var)
# v_fit<-fit.variogram(emp_var,vgm("Gau"))
# plot(emp_var,v_fit)
#   
# H <- as.matrix(dist(tmp$train[,c("x", "y")]))
# hist(H[upper.tri(H)])
# 
# H_obs_noobs <- as.matrix(rdist(tmp$train[,c("x", "y")], tmp$test[,c("x", "y")]))
# hist(H_obs_noobs[upper.tri(H_obs_noobs)])

S_obs_inv <- solve(theta$Sigma_obs, tol = 1e-22)
det(theta$Sigma_obs)
S_noobs_tr <- t(theta$Sigma_obs_noobs)
zhat <- S_noobs_tr %*% S_obs_inv %*% z
hist(zhat)

pz <- pnorm(zhat)
hist(pz)
# pz[which(pz==1)] <- 1-.Machine$double.eps
# pz[which(pz==0)] <- .Machine$double.eps
pz2 <- (pz - rep(transformed.training.set$Pi[1], length(pz)))/(1 - rep(transformed.training.set$Pi[1], length(pz)))
hist(pz2)

# Need to redefine qzig() function to handle case where none of the standard uniforms are less than the value of Pi

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

zigs <- qzig(u = pz,
           mu = rep(theta$mu[1], length(pz)),
           beta = theta$beta,
           epsilon = transformed.training.set$epsilon,
           Pi = rep(transformed.training.set$Pi[1], length(pz)))
hist(zigs)

yhat_gauscop.alg2 <- backtransform_gauscop(zigs, transformed.training.set$epsilon)
hist(yhat_gauscop.alg2)
hist(yhat_gauscop.krige)

rmspe(yhat_gauscop.alg2, tmp$test$resp)
rmspe(yhat_gauscop.krige, tmp$test$resp)


ix <- which(pz <= transformed.training.set$Pi[1])
result <- numeric(length=length(pz))
result[ix] <- qunif(pz[ix]/transformed.training.set$Pi[ix],0,transformed.training.set$epsilon) 
result[-ix] <- qsgamma((pz[-ix]-transformed.training.set$Pi[-ix])/(1-Pi[-ix]),mu[-ix],beta,epsilon)
return(result)
