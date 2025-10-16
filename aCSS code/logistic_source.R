


generate_parameters = function(){
	parameters = list()
	parameters$d = 5
	parameters$theta0 = rep(0.2,5)
	parameters$sigma = sqrt(10)
	parameters$tol = 1e-4
	parameters$signal = seq(from=0,to=1,length.out=11)
	parameters$try_mcmc_params = lapply(seq(from=5,to=25,by=5),list) 
	parameters$try_mcmc_params_ntrial = 50
	
	# example specific parameters
	parameters$example = list()
	parameters$example$n = 100
	parameters$example$coefyZ = rep(0.1,5)
	parameters$example$numstep_mult = 2
	parameters$example$numstep_max = 500
	parameters$example$acceptance_prob_min = 0.05
	parameters$example$nonlinear_fn = function(t){t + 0.5*t^3}
	parameters$example$elementwise_Z_fn = function(t){matrix(pmax(0,t),nrow(as.matrix(t)),ncol(as.matrix(t)))}
	parameters$example$coefyZ_baseline = rep(0.5,5)
	parameters$example$coefyZ_group0add = c(1,0,0,0,0)
	parameters$example$coefyZ_group1add = c(0,0,0,0,1)
	
	parameters
}

generate_experiment = function(isignal,parameters){
	experiment = list()
	
	Z = matrix(rnorm(parameters$example$n*parameters$d),parameters$example$n,parameters$d)
	X = rbinom(parameters$example$n,1,c(1/(1+exp(-Z%*%parameters$theta0))))
	coefyZ_X0 = parameters$signal[isignal] * parameters$example$coefyZ_group0add
	coefyZ_X1 = parameters$signal[isignal] * parameters$example$coefyZ_group1add
	baseline = parameters$example$elementwise_Z_fn(Z)%*%parameters$example$coefyZ_baseline
	Y = parameters$example$nonlinear_fn(baseline + (Z%*%coefyZ_X0)*(X==0) + (Z%*%coefyZ_X1)*(X==1)) + rnorm(parameters$example$n)
	experiment$X = X
	
	# example specific random variables
	experiment$example = list()
	experiment$example$Z = Z
	experiment$example$Y = Y

	experiment
}

Tfun = function(x,parameters,experiment){
	y0_lo = quantile(experiment$example$Y[x==0],1/3)
	y0_hi = quantile(experiment$example$Y[x==0],2/3)
	y1_lo = quantile(experiment$example$Y[x==1],1/3)
	y1_hi = quantile(experiment$example$Y[x==1],2/3)
	Z0_lo = colMeans(experiment$example$Z[(x==0)&(experiment$example$Y<=y0_lo),]) - colMeans(experiment$example$Z[x==0,])
	Z0_med = colMeans(experiment$example$Z[(x==0)&(experiment$example$Y>y0_lo)&(experiment$example$Y<y0_hi),]) - colMeans(experiment$example$Z[x==0,])
	Z0_hi = colMeans(experiment$example$Z[(x==0)&(experiment$example$Y>=y0_hi),]) - colMeans(experiment$example$Z[x==0,])
	Z1_lo = colMeans(experiment$example$Z[(x==1)&(experiment$example$Y<=y1_lo),]) - colMeans(experiment$example$Z[x==1,])
	Z1_med = colMeans(experiment$example$Z[(x==1)&(experiment$example$Y>y1_lo)&(experiment$example$Y<y1_hi),]) - colMeans(experiment$example$Z[x==1,])
	Z1_hi = colMeans(experiment$example$Z[(x==1)&(experiment$example$Y>=y1_hi),]) - colMeans(experiment$example$Z[x==1,])

	betah0 = svd(cbind(Z0_lo,Z0_med,Z0_hi),nu=1)$u[,1]
	betah1 = svd(cbind(Z1_lo,Z1_med,Z1_hi),nu=1)$u[,1]
	1 - abs(sum(betah0*betah1))
}
generate_X_null = function(theta,parameters,experiment){
	rbinom(parameters$example$n,1,c(1/(1+exp(-experiment$example$Z%*%theta))))
}
negloglik = function(x,theta,parameters,experiment){
	probs = c(1/(1+exp(-experiment$example$Z%*%theta)))
	-sum(x*log(probs)+(1-x)*log(1-probs))
}
gradient = function(x,theta,parameters,experiment){
	probs = c(1/(1+exp(-experiment$example$Z%*%theta)))
	c((probs-x)%*%experiment$example$Z)
}
hessian = function(x,theta,parameters,experiment){
	probs = c(1/(1+exp(-experiment$example$Z%*%theta)))
	t(experiment$example$Z)%*%diag(probs*(1-probs))%*%experiment$example$Z
}
thetahat = function(x,w,parameters,experiment){
	Lw = function(theta){negloglik(x,theta,parameters,experiment)+parameters$sigma*sum(w*theta)}
	gradLw = function(theta){gradient(x,theta,parameters,experiment)+parameters$sigma*w}
	theta_init = glm(x~0+experiment$example$Z,family='binomial')$coef
	optim(theta_init,Lw,gradLw,control=list(reltol=1e-12))$par
}

mcmc_proposal_draw = function(x,theta,mcmc_params,parameters,experiment){
	S = sample(parameters$example$n,mcmc_params[[1]])
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

score_test = function(x,parameters,experiment){NA}


