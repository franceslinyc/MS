# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for splitting the data and performing the cube root transformation / replace 0s with runifs values
# as well as the code for backtransforming the transformation

split_data <- function(full_data, n = 300) {
    
    #' Function to split data frame into a training set and a test set
    #' 
    #' @param full_data dataframe to split
    #' @param n integer 
    #' @return output list of two datafram
  
    colnames(full_data) <- c("x", "y", "resp")
    training_idx <- sample(1:nrow(full_data), n, replace=FALSE)
    training_data <- full_data[training_idx,]
    test_data <- full_data[-training_idx,]
    output <- list(train = training_data, 
                 test = test_data)
    return(output)
}

cube_root_fix_zeros <- function(data, cube.root.transform = TRUE, correct.zeros = TRUE) {
   
    #' Function to perform zero correction and cube root transform 
    #' and output the proportion of zeros and smallest non-zero value 
    #' 
    #' @param data dataframe to perform transformation on
    #' @param cube.root.transform boolean value, defaults to TRUE
    #' @param correct.zeros boolean value, defaults to TRUE
    #' 
    #' @return list containing transformed data, logical vector of zero values in responses, smallest nonzero, and vector of proportions of zero
    
    if (!("resp" %in% colnames(data))) {
    stop("Needs response column in dataset")
    }
    
    if (cube.root.transform) {
    data$resp <- data$resp^(1/3)
    }
    zeros <- which(data$resp == 0)
    if (length(zeros) > 0){
      epsilon <- min(data$resp[-zeros])
      Pi <- rep(length(zeros)/nrow(data), nrow(data))
    } else {
      # Needs a continuity correction
      epsilon <- min(data$resp)
      Pi <- rep(.001, nrow(data))
    }
    
    if (correct.zeros) {
      data$resp[zeros] <- runif(length(zeros), 0, epsilon)
    }
    
    output <- list(transformed_data = data,
                 zeros = zeros,
                 epsilon = epsilon,
                 Pi = Pi)
    return(output)
}

backtransform_gauscop <- function(preds, epsilon) {
    
    #' Function to back transform cube root and replace small values with 0
    #' 
    #' @param pred vector of predicted values
    #' @param epsilon smallest nonzero value in original data
    #' 
    #' @return vector of backtransformed predicted values

    preds[preds < epsilon] <- 0
    return(preds^3)
}

small_preds_to_zeros <- function(preds, epsilon, backtransform=TRUE) {
  
  #' Function to transform small predicted values (below a given epsilon) to 0
  #' 
  #' @param preds vector of predicted values
  #' @param epsilon smallest nonzero value in original data
  #' 
  #' @return vector of new predicted values
  
  if (backtransform) {
    epsilon <- epsilon^3
    preds[preds < epsilon] <- 0
  } else {
    preds[preds < epsilon] <- 0
  }
  
  return(preds)
}
