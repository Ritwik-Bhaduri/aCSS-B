library(mixtools)
generate_parameters = function(){
  parameters = list()
  
  # example specific parameters
  parameters$example = list()
  parameters$example$m = 10
  parameters$example$n = 10
  
  parameters$sigma = 0.5 
  # the sampling distribution for \tilde X is written for sigma = 0.5
  # If sigma is changed one would need to change functions for log_posterior and its (double) derivatives
  parameters$U = rnorm(parameters$example$m, 0, 1)
  parameters$V = rnorm(parameters$example$n, 0, 1)
  
  parameters$signal = seq(0,1,length.out = 8)
  
  parameters$non_null$U = rnorm(parameters$example$m, 0, 1)
  parameters$non_null$V = rnorm(parameters$example$n, 0, 1)
  parameters$non_null$mat = parameters$non_null$U %*% t(parameters$non_null$V)
  return(parameters)
}

generate_experiment = function(isignal,parameters){
  experiment = list() 
  noise = matrix(rnorm(parameters$example$m*parameters$example$n, 0, parameters$sigma), nrow = parameters$example$m)
  X = parameters$U %*% t(parameters$V) + 
    parameters$signal[isignal]*parameters$non_null$mat +  
    noise
  experiment$X = X  
  return(experiment)
}

generate_X_null = function(parameters){
  experiment = list() 
  X = parameters$U %*% t(parameters$V)
  noise = matrix(rnorm(parameters$example$m*parameters$example$n, 0, parameters$sigma), nrow = parameters$example$m)
  experiment$X = X + noise
  return(experiment)
}

Tfun = function(X,parameters,experiment){ 
  svd_X = svd(X)
  return(svd_X$d[2]^2)
}

