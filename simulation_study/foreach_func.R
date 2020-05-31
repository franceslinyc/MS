# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Simulation function for use in foreach

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

# alb.xy <- project(with(dat,cbind(lon_fuzzed,lat_fuzzed)),
#                   "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# colnames(alb.xy) <- c("x","y")
# dat <- cbind(alb.xy,dat)

run.sim <- function(index, dat, train.n = 300) {
    
    set.seed(182)
    
    full_simdata <- data.frame(
        x = dat$x,
        y = dat$y,
        resp = simdata[index,]
    )
    
    theta <- NULL
    rfsp_model <- NULL
    
    # Create spatial random forest model
    while(is.null(rfsp_model)){
        tmp <- split_data(full_simdata, n = train.n)
        try(
            calculate_rfsp_predictions(training_data = tmp$train,
                                       test_points = tmp$test) -> rfsp_model
        )
    }
    
    # Create spatial Gaussian copula model
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
                       test_data = tmp$test) -> zhat.zigs.df
    yhat_gauscop <- backtransform_gauscop(zhat.zigs.df$zigs, transformed.training.set$epsilon)
    
    # Create ordinary kriging model
    training_data.cpy <- tmp$train
    test_data.cpy <- tmp$test
    coordinates(training_data.cpy) <- ~x+y
    coordinates(test_data.cpy) <- ~x+y
    
    ak_model <- autoKrige(resp ~ 1, training_data.cpy, test_data.cpy)
    
    run.sim.output <- list(
        obs_responses = tmp$test$resp,
        gauscop_predictions = yhat_gauscop,
        rfsp_predictions = rfsp_model$predictions,
        ok_predictions = ak_model$krige_output$var1.pred
    ) 
    
    return(run.sim.output)
}

