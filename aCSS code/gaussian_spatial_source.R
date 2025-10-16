


generate_parameters = function(){
	parameters = list()
	parameters$d = 1
	parameters$theta0 = 0.2
	parameters$sigma = 1
	parameters$tol = 1e-5
	parameters$signal = seq(from=0,to=1,by=0.2)
	parameters$try_mcmc_params = lapply(seq(from=0.99,to=0.75,by=-0.02),list)
	parameters$try_mcmc_params_ntrial = 50
	
	# example specific parameters
	parameters$example = list()
	parameters$example$numstep_mult = 20
	parameters$example$numstep_max = 500
	parameters$example$acceptance_prob_min = 0.05
	parameters$example$theta_min = 0.001
	grid_size = 10; n = grid_size^2
	parameters$example$n = n
	points_loc = cbind(rep(1:grid_size,grid_size),rep(1:grid_size,each=grid_size))
	Dmat_dim1 = sqrt(points_loc[,1]^2%*%t(rep(1,n)) + rep(1,n)%*%t(points_loc[,1]^2) 
			- 2*outer(points_loc[,1],points_loc[,1]))
	Dmat_dim2 = sqrt(points_loc[,2]^2%*%t(rep(1,n)) + rep(1,n)%*%t(points_loc[,2]^2) 
			- 2*outer(points_loc[,2],points_loc[,2]))
	parameters$example$Dmat0 = sqrt(Dmat_dim1^2 + Dmat_dim2^2)
	parameters$example$dist_thresh = 1
	cutoff = 5
	parameters$example$split1 = which(points_loc[,2]<=cutoff)
	parameters$example$split2 = which(points_loc[,2]>cutoff)

	parameters
}

generate_experiment = function(isignal,parameters){
	experiment = list()
	Sigma = exp(-parameters$example$Dmat0*parameters$theta0)
	Sigma[parameters$example$split1,parameters$example$split2] = (1-parameters$signal[isignal])*Sigma[parameters$example$split1,parameters$example$split2]
	Sigma[parameters$example$split2,parameters$example$split1] = (1-parameters$signal[isignal])*Sigma[parameters$example$split2,parameters$example$split1]
	experiment$X = c(t(chol(Sigma))%*%rnorm(parameters$example$n))

	# example specific random variables
	experiment$example = list()

	experiment
}

Tfun = function(x,parameters,experiment){
        W = exp(-abs(x%*%t(rep(1,parameters$example$n))-t(x%*%t(rep(1,parameters$example$n)))))
	W[parameters$example$Dmat0>parameters$example$dist_thresh] = 0
	D = diag(rowSums(W))
	cluster1 = which(kmeans(eigen(D-W)$vectors[,parameters$example$n-1],2)$cluster==1)
	c(-sum(W[cluster1,-cluster1])*(1/sum(W[cluster1,]) + 1/sum(W[-cluster1,])))
}
generate_X_null = function(theta,parameters,experiment){
	c(t(chol(exp(-parameters$example$Dmat0*theta)))%*%rnorm(parameters$example$n))
}
negloglik = function(x,theta,parameters,experiment){
	0.5 * log(det(exp(-parameters$example$Dmat0*theta))) + 
		0.5 * t(x)%*%solve(exp(-parameters$example$Dmat0*theta),x)
}
gradient = function(x,theta,parameters,experiment){
	Siginv_x = solve(exp(-parameters$example$Dmat0*theta),x)
	Sig = exp(-parameters$example$Dmat0*theta)
	c(0.5 * sum(diag(solve(Sig,-parameters$example$Dmat0*Sig))) - 
		0.5 * t(Siginv_x)%*%(-parameters$example$Dmat0*Sig)%*%Siginv_x)
}
hessian = function(x,theta,parameters,experiment){
	Siginv_x = solve(exp(-parameters$example$Dmat0*theta),x)
	Sig = exp(-parameters$example$Dmat0*theta)
	matrix(c(0.5 * t(Siginv_x)%*%(-parameters$example$Dmat0^2*Sig 
			+ 2*(parameters$example$Dmat0*Sig)%*%solve(Sig,parameters$example$Dmat0*Sig))%*%Siginv_x
			+ 0.5 * sum(diag(solve(Sig,parameters$example$Dmat0^2*Sig))) 
			- 0.5 * sum(diag(solve(Sig,parameters$example$Dmat0*Sig)%*%solve(Sig,parameters$example$Dmat0*Sig)))),1,1)
}
thetahat = function(x,w,parameters,experiment){
	tmp = sum(outer(x,x)*(parameters$example$Dmat0==1))/sum(parameters$example$Dmat0==1)
	if(tmp > 0){
		theta_init = max(-log(tmp),parameters$example$theta_min)
	}else{
		theta_init = parameters$example$theta_min
	}
	optim(theta_init,function(theta){
		negloglik(x,theta,parameters,experiment) + parameters$sigma*w*theta},
		gr=function(theta){gradient(x,theta,parameters,experiment)+parameters$sigma*w},
		method='L-BFGS-B',lower=parameters$example$theta_min)$par
}

mcmc_proposal_draw = function(x,theta,mcmc_params,parameters,experiment){
	rho = mcmc_params[[1]]
	c(rho*x + sqrt(1-rho^2)*t(chol(exp(-parameters$example$Dmat0*theta)))%*%rnorm(parameters$example$n))
}
	
mcmc_acceptance_prob = function(x,x_new,theta,mcmc_params,parameters,experiment){
	w_new = -gradient(x_new,theta,parameters,experiment)/parameters$sigma
	thetah_new = thetahat(x_new,w_new,parameters,experiment)
	check = check_SSOSP(x_new,w_new,theta,parameters,experiment)&
		(((thetah_new)-theta)^2<=parameters$tol)
	if(check){
		acceptance_prob = min(1,
			c(exp(-gradient(x_new,theta,parameters,experiment)^2/2/parameters$sigma^2)/
			exp(-gradient(x,theta,parameters,experiment)^2/2/parameters$sigma^2) *
			(hessian(x_new,theta,parameters,experiment))/(hessian(x,theta,parameters,experiment))))
	}else{
		acceptance_prob = 0
	}
	acceptance_prob
}

choose_mcmc_params_rule = function(x_list,x_new_list,num_accept,num_trial,parameters,experiment){

	avg_accept = num_accept/num_trial
	avg_corr = rep(0,length(parameters$try_mcmc_params))
	for(iparam in 1:length(parameters$try_mcmc_params)){
		for(itrial in 1:num_trial){
			avg_corr[iparam] = avg_corr[iparam] +
				cor(x_list[[itrial]],x_new_list[[iparam]][[itrial]])/num_trial
		}
	}
	avg_corr_clip = avg_corr
	avg_corr_clip[which(avg_accept<parameters$example$acceptance_prob_min)] = Inf
	iparam_best = which.min(avg_corr_clip)
	output = list()
	output$mcmc_params = list()
	output$mcmc_params = parameters$try_mcmc_params[[iparam_best]]
	output$mcmc_numstep = min(parameters$example$numstep_max, 
		round(parameters$example$numstep_mult / (1-max(0,avg_corr[iparam_best]))))
	output
}

score_test = function(x,parameters,experiment){NA}