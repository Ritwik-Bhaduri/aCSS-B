generate_knots = function(range_Z = c(-5,5)){
  knot_vec = seq(range_Z[1], range_Z[2], length.out = 4)
  knot_vec = knot_vec[-c(1,length(knot_vec))]
  return(knot_vec)
}

generate_beta = function(signal, intercept = 0){
  target_beta_vec = c(-3, 1, 1 - signal) # denotes the final slopes of the linear spline
  
  # adjust beta so that cumulative beta matches above beta_vec
  beta_vec = c(target_beta_vec[1], target_beta_vec[2] - target_beta_vec[1], target_beta_vec[3] - target_beta_vec[2])
  beta_vec = c(intercept, beta_vec)
  return(beta_vec)
}

generate_parameters = function(){
  parameters = list()
  parameters$range_Z = c(-5,5)
  parameters$k = 2        # denotes number of segments under the null, current implementation only works for one knot, i.e., k=2
  parameters$signal = seq(from=0,to=1.8, length.out = 8)
  parameters$knots = generate_knots(parameters$range_Z)
  parameters$beta = lapply(parameters$signal, function(sd) generate_beta(sd))
  parameters$sigma = 0.5
  
  # example specific parameters
  parameters$n = 50 
  parameters
}

generate_experiment = function(isignal,parameters){
  experiment = list()
  Z = runif(parameters$n, min = parameters$range_Z[1], max = parameters$range_Z[2])
  h_t <- cbind(1,Z,pmax(outer(Z, parameters$knots, "-"),0))
  beta = parameters$beta[[isignal]]
  X = h_t %*% beta + rnorm(parameters$n,0,parameters$sigma)
  experiment$X = X
  experiment$example = list()
  experiment$example$Z = Z
  experiment$example$beta = beta
  experiment$example$t = parameters$knots
  return(experiment)
}

generate_X_null = function(parameters,experiment){
  Z = experiment$example$Z
  t = parameters$knots
  h_t <- cbind(1,Z,pmax(outer(Z,t,"-"),0))
  beta = parameters$beta[[1]]
  X = h_t %*% beta + rnorm(parameters$n,0,parameters$sigma)
  experiment = list() 
  experiment$X = X  
  experiment$example$Z = Z
  return(experiment)
}

Tfun = function(x, parameters, experiment){
  Z = experiment$example$Z
  initial_model <- lm(x ~ Z)
  segmented_model <- segmented(initial_model, seg.Z = ~Z, npsi = 1)
  x_pred <- predict(segmented_model)
  return(sum((x-x_pred)^2))
}