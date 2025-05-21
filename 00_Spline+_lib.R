# TODO: make this a more general usage for samples?

# Find Cubic Smooth Spline of 1D Data
# INPUT: x_new is a set of 
# INPUT: uses a df with a sample,holding data in sample_df$mean as the x variable
find_spline <- function(sample_name, samplewise_df, x_new) {
  print(sample_name)
  sample_df <- samplewise_df[samplewise_df$sample == sample_name,]
  
  # remove NaN's:
  sample_df <- sample_df[!is.nan(sample_df$mean),]
  
  # Create a smoothing spline for Sample:
  smooth_spline <- smooth.spline(sample_df$distance, sample_df$mean, cv = FALSE) # spar controls the amount of 
  
  # Estimating Function with Spline:
  y_hat <- predict(smooth_spline, x_new)$y # Evaluate the spline at the new x-values
  
  # Note that predict may extrapolate beyond last measurement. Set such points to NaN
  y_hat[x_new > max(sample_df$distance)] <- NaN
  
  y_hat
}

find_derivative <- function(sample_name, samplewise_df, x_new) {
  print(sample_name)
  sample_df <- samplewise_df[samplewise_df$sample == sample_name,]
  
  # remove NaN's:
  sample_df <- sample_df[!is.nan(sample_df$mean),]
  
  # Create a smoothing spline for Sample:
  smooth_spline <- smooth.spline(sample_df$distance, sample_df$mean, cv = FALSE) # spar controls the amount of 
  
  # Estimating Derivative from Spline:
  y_hat_prime <- predict(smooth_spline, x_new, deriv = 1)$y # deriv = 1 for first derivative
  
  # Note that predict may extrapolate beyond last measurement. Set such points to NaN
  y_hat_prime[x_new > max(sample_df$distance)] <- NaN
  
  y_hat_prime
}

find_nth_derivative <- function(sample_name, samplewise_df, x_new, deriv) {
  print(sample_name)
  sample_df <- samplewise_df[samplewise_df$sample == sample_name,]
  
  # remove NaN's:
  sample_df <- sample_df[!is.nan(sample_df$mean),]
  
  # Create a smoothing spline for Sample:
  smooth_spline <- smooth.spline(sample_df$distance, sample_df$mean, cv = FALSE) # spar controls the amount of 
  
  # Estimating Derivative from Spline:
  y_hat_prime <- predict(smooth_spline, x_new, deriv = deriv)$y # nth derivative
  
  # Note that predict may extrapolate beyond last measurement. Set such points to NaN
  y_hat_prime[x_new > max(sample_df$distance)] <- NaN
  
  y_hat_prime
}

# Find finite difference in samplewise_df$mean
find_relative_change <- function(sample_name, samplewise_df) {
  print(sample_name)
  sample_df <- samplewise_df[samplewise_df$sample == sample_name,]
  
  y <- sample_df$mean
  x <- sample_df$distance
  # Estimating Derivative using first difference:
  y_delta <- y[2:length(y)] -  y[1:length(y)-1]
  x_delta <- x[2:length(y)] - x[1:length(y)-1]
  y_hat_prime <- y_delta / x_delta
  x_prime <- (x[1:length(y)-1] + x[2:length(y)]) / 2
  
  names(y_hat_prime) <- x_prime * 10
  
  y_hat_prime
}