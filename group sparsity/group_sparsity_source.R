library(gglasso)

get_grp_sizes = function(p, G){
  equal_sizes = rep(floor(p/G), G)
  remainder = p - G * floor(p/G)
  extra_sizes = rep(0, G)
  extra_sizes[sample(1:G, size = remainder, replace = FALSE)] = 1
  sizes = equal_sizes + extra_sizes
  return(sizes)
}

generate_parameters = function(){
  parameters = list()
  
  parameters$n = 100
  parameters$p = 50
  parameters$G = 10
  parameters$sigma = 0.5
  parameters$n_signal = 7 # total number of signals including null(signal = 0)
  parameters$signal = seq(0, 0.6, length.out = parameters$n_signal)
  
  parameters$grp_sizes = get_grp_sizes(parameters$p, parameters$G)
  parameters$grp_indices = split(1:parameters$p, unlist(sapply(1:parameters$G, function(g) rep(g, parameters$grp_sizes[g]))))
  parameters$beta = c(rnorm(parameters$grp_sizes[1], 0, 1), rep(0, sum(parameters$grp_sizes[-1]))) 
  
  parameters$non_null = list()
  parameters$non_null$beta = matrix(NA, nrow = parameters$p, ncol = parameters$n_signal - 1)
  
  second_group_coef = rnorm(parameters$grp_sizes[2], 0, 1)
  for(j in 1:(parameters$n_signal - 1)){
    if(j == 1){
      parameters$non_null$beta[,1] = parameters$beta
    }else{
      parameters$non_null$beta[,j] = parameters$non_null$beta[,j-1]
    }
    parameters$non_null$beta[parameters$grp_indices[[2]],j] = parameters$signal[j+1] * second_group_coef
  }
  
  return(parameters)
}

generate_experiment = function(isignal,parameters){
  if(isignal == 1){
    beta = parameters$beta
  }else{
    beta = parameters$non_null$beta[,isignal-1]
  }
  
  Z = matrix(rnorm(parameters$n*parameters$p),parameters$n,parameters$p)
  X = Z %*% beta + rnorm(parameters$n, 0, parameters$sigma)
    
  experiment = list() 
  experiment$X = X  
  experiment$example$Z = Z
  return(experiment)
}

generate_X_null = function(parameters,experiment){
  Z = experiment$example$Z 
  X = Z %*% parameters$beta + rnorm(parameters$n, 0, parameters$sigma)
  
  experiment = list() 
  experiment$X = X  
  experiment$example$Z = Z
  return(experiment)
}

Tfun = function(X,parameters,experiment){
  grp = unlist(sapply(1:parameters$G, function(g) rep(g, parameters$grp_sizes[g])))
  lambda_best = cv.gglasso(experiment$example$Z, X, group = grp, lambda = NULL, pred.loss = "L1", nfolds = 5)$lambda.min
  beta_hat = gglasso(experiment$example$Z, X, group = grp, lambda = lambda_best)$beta
  
  max_abs_beta = abs(unlist(lapply(split(beta_hat, grp), max)))
  max_beta_index = which.max(max_abs_beta)
  t = as.numeric(sum(max_abs_beta[-max_beta_index]) / max_abs_beta[max_beta_index])
  
  return(t)
}

