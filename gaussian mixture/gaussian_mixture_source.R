library(mixtools)
generate_parameters = function(){
  parameters = list()
  parameters$J = 2
  parameters$d = 3 * parameters$J - 1
  parameters$theta0 = c(0.5, 0.4, 0.1, -0.4, 0.1)
  parameters$third_center = 0
  parameters$third_sd = 0.1
  parameters$signal = seq(from=0,to=0.21,by=0.03) #seq(from=0,to=0.21,by=0.03)
  
  # example specific parameters
  parameters$example = list()
  parameters$example$n = 200
  parameters
}

generate_experiment = function(isignal,parameters){
  experiment = list() 
  X = rnormmix(parameters$example$n, 
               c(parameters$signal[isignal], (1 - parameters$signal[isignal])/2, (1 - parameters$signal[isignal])/2), 
               c(parameters$third_center, parameters$theta0[2], parameters$theta0[4]), 
               c(parameters$third_sd, parameters$theta0[3], parameters$theta0[5]))
  experiment$X = X  
  return(experiment)
}

generate_X_null = function(parameters){
  experiment = generate_experiment(1, parameters)
  return(experiment)
}

Tfun = function(x, parameters, experiment){ 
  kmeans_2 = kmeans(x,2,nstart = 30)
  kmeans_3 = kmeans(x,3,nstart = 30)
  (kmeans_2$tot.withinss - kmeans_3$tot.withinss)/kmeans_2$tot.withinss
}
