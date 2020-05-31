# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for using spatial random forests to predict values

library(ranger)
library(raster)
library(GSIF)

calculate_rfsp_predictions <- function(training_data, test_points, num.trees = 150, tolerance = .7) {
    
    #' Function for using RFSP to calculate predictions for unobserved data on untransformed data
    #' 
    #' @param training_set dataframe of observed locations, including observed responses
    #' @param test_points dataframe of unobserved locations
    #' 
    #' @return ranger output with the predictions for the unobserved locations
  
    if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
        stop("Need to specify coordinates using x and y in training set")
    }
    if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
        stop("Need to specify coordinates using x and y in test set")
    }
    
    training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
    suppressWarnings(
      training_points_spdf <- SpatialPixelsDataFrame(points = training_points, data = training_points@data, tolerance = tolerance)
    )
    
    test_points <- SpatialPointsDataFrame(test_points[,c("x", "y")], data.frame(resp = test_points[,c("resp")]))
    suppressWarnings(
      test_points_spdf <- SpatialPixelsDataFrame(points = test_points, data = test_points@data, tolerance = tolerance)
    )
    
    suppressWarnings(
      grid.dist0 <- GSIF::buffer.dist(training_points, training_points_spdf, as.factor(1:nrow(training_data)))
    )
    
    dn0 <- paste(names(grid.dist0), collapse="+") 
    fm0 <- as.formula(paste("resp ~ ", dn0))   
    ov.resp <- over(training_points, grid.dist0)
    
    suppressWarnings(
      grid.dist1 <- GSIF::buffer.dist(training_points, test_points_spdf, as.factor(1:nrow(training_data)))
    )
    
    rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
    suppressWarnings(
    resp_ranger <- ranger(fm0,
                          rm.resp,
                          quantreg=TRUE,
                          num.trees=num.trees) 
    )
    
    resp.predict <- predict(resp_ranger, grid.dist1@data)
    return(resp.predict)
}

calculate_rfsp_predictions_rasterpd <- function(training_data, test_points, num.trees = 150) {
    
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
    if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
        stop("Need to specify coordinates using x and y in test set")
    }
    
    training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
    test_points <- SpatialPointsDataFrame(test_points[,c("x", "y")], data.frame(resp = test_points[,c("resp")]))
    train_distances <- data.frame(pointDistance(training_points, lonlat = FALSE))
    test_distances <- data.frame(pointDistance(test_points, training_points, lonlat = FALSE))
    colnames(train_distances) <- paste0("layer.", 1:nrow(training_data))
    colnames(test_distances) <- paste0("layer.", 1:nrow(training_data))
    
    dn0 <- paste(names(train_distances), collapse="+") 
    fm0 <- as.formula(paste("resp ~ ", dn0))   
    
    rm.resp <- cbind(data.frame(resp = training_data$resp), train_distances) 
    suppressWarnings(
    resp_ranger <- ranger(fm0,
                          rm.resp,
                          quantreg=TRUE,
                          num.trees=num.trees) 
    )
    
    resp.predict <- predict(resp_ranger, test_distances)
    return(resp.predict)
}

# nonrasterized.rfsp <- calculate_rfsp_predictions(tmp$train, tmp$test)
# rasterized.rfsp <- calculate_rfsp_predictions_rasterpd(tmp$train, tmp$test)
# 
# hist(nonrasterized.rfsp$predictions)
# hist(rasterized.rfsp$predictions)
# 
# rmspe(nonrasterized.rfsp$predictions, tmp$test$resp)
# rmspe(rasterized.rfsp$predictions, tmp$test$resp)

#' calculate_rfsp_predictions_rasterized <- function(training_data, test_points, num.trees = 150, tolerance = .7) {
#'     
#'     #' Function for using RFSP to calculate predictions for unobserved data on untransformed data
#'     #' This function rasterizes the points first to create a pixel layer.
#'     #' This is in effort to avoid issues with regular grid spacing in the other function.
#'     #' DUMB IDEA DOESNT WORK
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
#'     # coordinates(training_data) <- ~x+y
#'     # coordinates(test_points) <- ~x+y
#'  
#'     training_points <- SpatialPointsDataFrame(training_data[,c("x", "y")], as.data.frame(training_data[,c("resp")]))
#'     training_rast <- raster()
#'     extent(training_rast) <- extent(training_points)
#'     rast2 <- rasterize(training_points, training_rast, training_points@data$`training_data[, c("resp")]`, fun = mean)
#'     training_points_spdf <- as(rast2, "SpatialPixelsDataFrame")
#'     
#'     suppressWarnings(
#'       grid.dist0 <- GSIF::buffer.dist(training_points, training_points_spdf, as.factor(1:nrow(training_data)))
#'     )
#'     grid.dist0
#'     
#'     dn0 <- paste(names(grid.dist0), collapse="+") 
#'     fm0 <- as.formula(paste("resp ~ ", dn0))   
#'     grid.dist0@proj4string <- training_points@proj4string
#'     training_points@proj4string <- grid.dist0@proj4string
#'     ov.resp <- over(training_points, grid.dist0)
#'     
#'     # test_points <- tmp$test
#'     test_points <- SpatialPointsDataFrame(test_points[,c("x", "y")], data.frame(resp = test_points[,c("resp")]))
#'     test_rast <- raster()
#'     extent(test_rast) <- extent(test_points)
#'     rast3 <- rasterize(test_points, test_rast, test_points@data$resp, fun = mean)
#'     test_points_spdf <- as(rast3, "SpatialPixelsDataFrame")
#'     
#'     suppressWarnings(
#'       grid.dist1 <- GSIF::buffer.dist(training_points, test_points_spdf, as.factor(1:nrow(training_data)))
#'     )
#'     
#'     rm.resp <- cbind(data.frame(resp = training_data$resp), ov.resp) 
#'     suppressWarnings(
#'     resp_ranger <- ranger(fm0,
#'                           rm.resp,
#'                           quantreg=TRUE,
#'                           num.trees=num.trees) 
#'     )
#'     
#'     resp.predict <- predict(resp_ranger, grid.dist1@data)
#'     return(resp.predict)
#' }

# start <- Sys.time()
# train.n <- 300
# sim.sets.n <- 10
# rfsp_predictions_raster <- vector(mode = "list",
#                           length = sim.sets.n)
# rfsp_predictions_nonraster <- vector(mode = "list",
#                           length = sim.sets.n)
# obs_responses <- vector(mode = "list",
#                         length = sim.sets.n)
# for (i in 1:sim.sets.n) {
# 
#     set.seed(182)
#     full_simdata <- data.frame(
#         x = dat$x,
#         y = dat$y,
#         resp = simdata[i,]
#     )
# 
#     calculate_rfsp_predictions_rasterpd(training_data = tmp$train, test_points = tmp$test) -> rfsp_model_raster
#     calculate_rfsp_predictions(training_data = tmp$train, test_points = tmp$test) -> rfsp_model_nonraster
# 
#     rfsp_predictions_raster[[i]] <- rfsp_model_raster$predictions
#     rfsp_predictions_nonraster[[i]] <- rfsp_model_nonraster$predictions
#     obs_responses[[i]] <- tmp$test$resp
# }
# end <- Sys.time()
# end - start
# 
# rmspe(unlist(rfsp_predictions_raster), unlist(obs_responses))
# rmspe(unlist(rfsp_predictions_nonraster), unlist(obs_responses))
