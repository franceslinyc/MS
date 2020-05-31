# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for using spatial random forests to predict values

library(ranger)
library(raster)
library(GSIF)

calculate_rfsp_predictions_rasterpd <- function(training_data, test_data, num.trees = 150) {
    
    #' Function for using RFSP to calculate predictions for unobserved data on untransformed data
    #' Uses raster::pointDistances instead of GSIF::buffer.dist
    #' 
    #' @param training_set dataframe of observed locations, including observed responses
    #' @param test_points dataframe of unobserved locations
    #' 
    #' @return ranger output with the predictions for the unobserved locations
  
    if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
        stop("Need to specify coordinates using x and y in training set")
    }
    if (!("x" %in% colnames(test_data) && "y" %in% colnames(test_data))){
        stop("Need to specify coordinates using x and y in test set")
    }
    
    training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
    test_points <- SpatialPointsDataFrame(test_data[,c("x", "y")], data.frame(resp = test_data[,c("resp")]))
    train_distances <- data.frame(pointDistance(training_points, lonlat = FALSE))
    train_distances$annpre <- training_data$annpre
    test_distances <- data.frame(pointDistance(test_points, training_points, lonlat = FALSE))
    test_distances$annpre <- test_data$annpre
    colnames(train_distances) <- c(paste0("layer.", 1:nrow(training_data)), "annpre")
    colnames(test_distances) <- c(paste0("layer.", 1:nrow(training_data)), "annpre")
    
    dn0 <- paste(names(train_distances), collapse="+") 
    fm0 <- as.formula(paste("resp ~ ", dn0))   
    
    rm.resp <- cbind(data.frame(resp = training_data$resp), train_distances) 
    suppressWarnings({
      resp_ranger <- ranger(fm0,
                            rm.resp,
                            quantreg=TRUE,
                            num.trees=num.trees) 
    })
    
    resp.predict <- predict(resp_ranger, test_distances)
    return(resp.predict)
}

#' calculate_rfsp_predictions <- function(training_data, test_points, num.trees = 150, tolerance = .7) {
#'     
#'     #' Function for using RFSP to calculate predictions for unobserved data on untransformed data
#'     #' 
#'     #' @param training_set dataframe of observed locations, including observed responses
#'     #' @param test_points dataframe of unobserved locations
#'     #' 
#'     #' @return ranger output with the predictions for the unobserved locations
#'   
#'     if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
#'         stop("Need to specify coordinates using x and y in training set")
#'     }
#'     if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
#'         stop("Need to specify coordinates using x and y in test set")
#'     }
#'     
#'     training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
#'     suppressWarnings(
#'       training_points_spdf <- SpatialPixelsDataFrame(points = training_points, data = training_points@data, tolerance = tolerance)
#'     )
#'     
#'     test_points <- SpatialPointsDataFrame(test_points[,c("x", "y")], data.frame(resp = test_points[,c("resp")]))
#'     suppressWarnings(
#'       test_points_spdf <- SpatialPixelsDataFrame(points = test_points, data = test_points@data, tolerance = tolerance)
#'     )
#'     
#'     suppressWarnings(
#'       grid.dist0 <- GSIF::buffer.dist(training_points, training_points_spdf, as.factor(1:nrow(training_data)))
#'     )
#'     
#'     dn0 <- paste(names(grid.dist0), collapse="+") 
#'     fm0 <- as.formula(paste("resp ~ ", dn0))   
#'     ov.resp <- over(training_points, grid.dist0)
#'     
#'     suppressWarnings(
#'       grid.dist1 <- GSIF::buffer.dist(training_points, test_points_spdf, as.factor(1:nrow(training_data)))
#'     )
#'     
#'     rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
#'     suppressWarnings({
#'       resp_ranger <- ranger(fm0,
#'                             rm.resp,
#'                             quantreg=TRUE,
#'                             num.trees=num.trees) 
#'     })
#'     
#'     resp.predict <- predict(resp_ranger, grid.dist1@data)
#'     return(resp.predict)
#' }