


generate_parameters = function(){
	parameters = list()
	parameters$d = 3
	parameters$theta0 = c(0,1,2)
	parameters$sigma = 1
	parameters$tol = 1e-6
	parameters$signal = seq(from=0,to=1,by=0.1)
	parameters$try_mcmc_params = lapply(seq(from=5,to=25,by=5),list) 
	parameters$try_mcmc_params_ntrial = 50
	
	# example specific parameters
	parameters$example = list()
	parameters$example$n0 = 50
	parameters$example$n1 = 50
	parameters$example$numstep_mult = 2
	parameters$example$numstep_max = 500
	parameters$example$acceptance_prob_min = 0.2

	parameters
}

generate_experiment = function(isignal,parameters){
	experiment = list()
	
	mu0 = parameters$theta0[1]
	mu1 = parameters$signal[isignal]
	sig2_0 =parameters$theta0[2]
	sig2_1 = parameters$theta0[3]
	
	experiment$X = c(mu0 + sqrt(sig2_0)*rnorm(parameters$example$n0), 
					mu1 + sqrt(sig2_1)*rnorm(parameters$example$n1))
	
	# example specific random variables
	experiment$example = list()

	experiment
}

Tfun = function(x,parameters,experiment){
	abs(mean(x[1:parameters$example$n0])-mean(x[parameters$example$n0 + (1:parameters$example$n1)]))
}
generate_X_null = function(theta,parameters,experiment){
	mu = theta[1]
	sig2_0 =theta[2]
	sig2_1 = theta[3]

	c(mu + sqrt(sig2_0)*rnorm(parameters$example$n0), 
	mu + sqrt(sig2_1)*rnorm(parameters$example$n1))
}
negloglik = function(x,theta,parameters,experiment){
	mu = theta[1]
	sig2_0 = theta[2]
	sig2_1 = theta[3]
	parameters$example$n0/2 *log(2*pi*sig2_0) + 
		1/2/sig2_0 * sum((x[1:parameters$example$n0]-mu)^2) +
		parameters$example$n1/2 *log(2*pi*sig2_1) + 
		1/2/sig2_1 * sum((x[parameters$example$n0 + (1:parameters$example$n1)]-mu)^2)
}
gradient = function(x,theta,parameters,experiment){
	mu = theta[1]
	sig2_0 = theta[2]
	sig2_1 = theta[3]
	c(-1/sig2_0*sum(x[1:parameters$example$n0]-mu)
		-1/sig2_1*sum(x[parameters$example$n0 + (1:parameters$example$n1)]-mu),
		parameters$example$n0/2/sig2_0-1/2/sig2_0^2*sum((x[1:parameters$example$n0]-mu)^2),
		parameters$example$n1/2/sig2_1-1/2/sig2_1^2*sum((x[parameters$example$n0 + (1:parameters$example$n1)]-mu)^2))
}
hessian = function(x,theta,parameters,experiment){
	mu = theta[1]
	sig2_0 = theta[2]
	sig2_1 = theta[3]
	matrix(c(parameters$example$n0/sig2_0+parameters$example$n1/sig2_1,
		1/sig2_0^2*sum(x[1:parameters$example$n0]-mu),
		1/sig2_1^2*sum(x[parameters$example$n0 + (1:parameters$example$n1)]-mu),
		1/sig2_0^2*sum(x[1:parameters$example$n0]-mu),
		-parameters$example$n0/2/sig2_0^2+1/sig2_0^3*sum((x[1:parameters$example$n0]-mu)^2),
		0,
		1/sig2_1^2*sum(x[parameters$example$n0 + (1:parameters$example$n1)]-mu),
		0,
		-parameters$example$n1/2/sig2_1^2+1/sig2_1^3*sum((x[parameters$example$n0 + (1:parameters$example$n1)]-mu)^2)),3,3)
}
thetahat = function(x,w,parameters,experiment){
	mu_init = mean(x)
	sig2_0_init = mean((x[1:parameters$example$n0]-mu_init)^2)
	sig2_1_init = mean((x[parameters$example$n0 + (1:parameters$example$n1)]-mu_init)^2)
	theta_init = c(mu_init,sig2_0_init,sig2_1_init)
	Lw = function(theta){negloglik(x,theta,parameters,experiment)+parameters$sigma*sum(w*theta)}
	gradLw = function(theta){gradient(x,theta,parameters,experiment)+parameters$sigma*w}
	optim(theta_init,Lw,gradLw,control=list(reltol=1e-12))$par
}

mcmc_proposal_draw = function(x,theta,mcmc_params,parameters,experiment){
	S = sample(parameters$example$n0+parameters$example$n1,mcmc_params[[1]])
	x_new = x
	x_new[S] = generate_X_null(theta,parameters,experiment)[S]
	x_new
}
	
mcmc_acceptance_prob = function(x,x_new,theta,mcmc_params,parameters,experiment){
	w_new = -gradient(x_new,theta,parameters,experiment)/parameters$sigma
	thetah_new = thetahat(x_new,w_new,parameters,experiment)
	check = check_SSOSP(x_new,w_new,theta,parameters,experiment)&(sum(((thetah_new)-theta)^2)<=parameters$tol)
	if(check){
		acceptance_prob = min(1,
			c(exp(-sum(gradient(x_new,theta,parameters,experiment)^2)/(2*parameters$sigma^2/parameters$d))/
				exp(-sum(gradient(x,theta,parameters,experiment)^2)/(2*parameters$sigma^2/parameters$d)) *
				det(hessian(x_new,theta,parameters,experiment))/det(hessian(x,theta,parameters,experiment))))
	}else{
		acceptance_prob = 0
	}
	acceptance_prob
}

choose_mcmc_params_rule = function(x_list,x_new_list,num_accept,num_trial,parameters,experiment){
	avg_accept = num_accept/num_trial
	avg_newsample = avg_accept * unlist(parameters$try_mcmc_params)
	avg_newsample_clip = avg_newsample; avg_newsample_clip[which(avg_accept<parameters$example$acceptance_prob_min)] = -Inf
	iparam_best = which.max(avg_newsample_clip)
	output = list()
	output$mcmc_params = list()
	output$mcmc_params = parameters$try_mcmc_params[[iparam_best]]
	output$mcmc_numstep = min(parameters$example$numstep_max, round(parameters$example$numstep_mult / (avg_newsample[iparam_best]/parameters$example$n)))
	output
}

score_test <- function(x,parameters,experiment){
	thetahat = thetahat(x,0,parameters,experiment)
	teststat = (mean(x[1:parameters$example$n0])-mean(x[parameters$example$n0+1:parameters$example$n1]))^2/(thetahat[2]/parameters$example$n0+thetahat[3]/parameters$example$n1)
	pval = 1-pchisq(teststat,1)
}


