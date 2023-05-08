calc_dev_threshold <- function(data_deviance_bias, max_bias, tau){
  
  # add a root term for the deviance to the bias vs. deviance data
  data_deviance_bias$dev_root <- sqrt(data_deviance_bias$deviance)
  colnames(data_deviance_bias) <- c("deviance", "bias", "dev_root")
  
  # fit polynomial quantile regression model for absolute bias
  quanfit <- rq(abs(bias) ~ deviance + dev_root -1, data = data_deviance_bias, tau = tau)
  
  # extract coefficient (assume intercept to be 0)
  b1 <- as.numeric(quanfit$coefficients[1])
  b2 <- as.numeric(quanfit$coefficients[2])
  
  # function to solve for bias
  f <- function(x) {
    b1*x + b2*sqrt(x) - max_bias
  }
  
  # calculate deviance threshold by solving function
  dev_threshold <- uniroot(f, c(0, max(data_deviance_bias$deviance)))$root
  
  return(dev_threshold)
  
}

