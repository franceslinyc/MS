# MSc Project 2020: Comparison of Spatial Gaussian Copula and Random Forests in Spatial Prediction
# Sun, N
#
# Code for calculating prediction metrics

rmspe <- function(predicted_values, observed_values) {
    
    #' Computes root mean squared prediction error 
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return RMSPE single value
    
    if(length(predicted_values) != length(observed_values)) {
        stop("Vectors are of unequal length.")
    }
    
    sq.diffs <- (predicted_values - observed_values)^2
    return(sqrt(mean(sq.diffs)))
     
}

mspe <- function(predicted_values, observed_values) {
    #' Computes sum of squared residuals
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return RMSPE single value
    
    if(length(predicted_values) != length(observed_values)) {
        stop("Vectors are of unequal length.")
    }
    
    sum.sq.diffs <- sum((predicted_values - observed_values)^2)
    return((sum.sq.diffs))
}

srb <- function(predicted_values, observed_values) {
  
    #' Computes Signed Relative Bias for predictions 
    #' using formula from ver Hoef (2013)
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return signed relative bias
  
    if(length(predicted_values) != length(observed_values)) {
        stop("Vectors are of unequal length.")
    }
    tau <- mean((predicted_values - observed_values))
    srb <- sign(tau)*sqrt(tau^2/(rmspe(predicted_values, observed_values)^2 - tau^2))
    return(srb) 
}

pic90 <- function(predicted_values, observed_values) {
  
    #' Computes 90% prediction interval coverage
    #' using formula from ver Hoef (2013).
    #' 
    #' The value produced here should be close to 90%.
    #' 
    #' @param predicted_values vector of predicted values
    #' @param observed_values vector of observed values
    #' 
    #' @return signed relative bias
    
    se <- sd(predicted_values)
    upper.bound <- predicted_values + 1.645*se
    lower.bound <- predicted_values - 1.645*se
    coverage <- as.integer((observed_values >= lower.bound) & (observed_values <= upper.bound))
    return(mean(coverage))
}

make_plot <- function(data) {
    
    ggplot(data, aes(x=x/1000,y=y/1000,color=resp)) +
      geom_point(size=2) +
      scale_color_gradient(low = "white", high = "red") +
      scale_x_continuous(breaks = round(seq(-2150, -2000, by = 50),2)) +
      theme(axis.text.x=element_text(angle=90,hjust=1))
    
}
