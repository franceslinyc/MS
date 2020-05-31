# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for using ordinary kriging to predict values, either using transformed ZIG or original values

library(geoR)
library(gstat)

ord_kriging_z <- function(training_data, test_points, test_z) {
    
    #' Perform ordinary kriging on z-hat standard normal values to get predictions for unobserved values.
    #' Makes use of function in the automap package.
    #' 
    #' @param training_data dataframe of observed locations
    #' @param test_points dataframe of unobserved locations
    #' @param test_z vector of transformed standard normal values
    #' 
    #' @return Ordinary kriging model
  
    if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
    stop("Need to specify coordinates using x and y")
    }
    if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
    stop("Need to specify coordinates using x and y")
    }
    
    training_data$z <- test_z
    test_points <- test_points[,c("x", "y")]
    coordinates(training_data) <- ~x+y
    coordinates(test_points) <- ~x+y
    
    ok_model_z <- autoKrige(z ~ 1,
                            input_data = training_data,
                            new_data = test_points)
    # z.vgm <- variogram(z ~ 1, training_data)
    # tryCatch({
    #     z.fit <- fit.variogram(z.vgm, vgm("Gau"))
    # }, warning = function(war) {
    #     stop("Variogram did not converge.")
    # })
    
    # tryCatch({
    #     invisible(ok_model_z <- krige(z ~ 1, training_data, test_points, model = z.fit))
    # }, warning = function(war) {
    #     stop("Ordinary Kriging not successful.")
    # })
    
    return(ok_model_z)
}

ord_kriging <- function(training_data, test_points) {
    
    #' Perform ordinary kriging on untransformed to get predictions for unobserved values.
    #' Makes use of function from gstat package.
    #' 
    #' @param training_data dataframe of observed locations
    #' @param test_points dataframe of unobserved locations
    #' 
    #' @return Ordinary kriging model
  
  if (!("x" %in% colnames(training_data) && "y" %in% colnames(training_data))){
    stop("Need to specify coordinates using x and y in training set")
  }
  if (!("x" %in% colnames(test_points) && "y" %in% colnames(test_points))){
    stop("Need to specify coordinates using x and y in test set")
  }
  
  test_points <- test_points[,c("x", "y")]
  coordinates(training_data) <- ~x+y
  coordinates(test_points) <- ~x+y
  
  resp.vgm <- variogram(resp ~ 1, training_data)
  tryCatch({
    resp.fit <- fit.variogram(resp.vgm, vgm("Gau"))
  }, warning = function(war) {
    stop("Variogram did not converge.")
  })
  
  tryCatch({
    invisible(ok_model <- krige(resp ~ 1, training_data, test_points, model = resp.fit))
  }, warning = function(war) {
    stop("Ordinary Kriging not successful.")
  })
  
  return(ok_model)
}

