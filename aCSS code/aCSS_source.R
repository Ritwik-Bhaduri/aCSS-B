# aCSS code for all examples

# Sources example-specific parameters & functions, 
# runs one trial of the simulation for each signal value, 
# and saves the resulting p-values (oracle method & aCSS method)
run_one_trial = function(M,seed,example, parameters= NULL){
	set.seed(seed)
	if(is.null(parameters)){
	  parameters = generate_parameters()
	  source(paste0("../aCSS - original code/", example,'_source.R'))
	}
	nsignal = length(parameters$signal)
	pval_aCSS = pval_oracle = pval_score = rep(0,nsignal)
	for(isignal in 1:nsignal){
		experiment = generate_experiment(isignal,parameters)
		Tfun_ = function(x){Tfun(x,parameters,experiment)}

		# oracle method
		Xcopies = list()
		for(m in 1:M){
			Xcopies[[m]] = generate_X_null(parameters$theta0,parameters,experiment)
		}
		pval_oracle[isignal] = 
			(1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+M)

		# score test (if applicable)
		pval_score[isignal] = score_test(experiment$X,parameters,experiment)

		# aCSS method
		W = rnorm(parameters$d)/sqrt(parameters$d)
		thetah = thetahat(experiment$X,W,parameters,experiment)
		if(check_SSOSP(experiment$X,W,thetah,parameters,experiment)){
			mcmc_params_tmp = choose_mcmc_params(thetah,parameters,experiment)
			mcmc_params = mcmc_params_tmp$mcmc_params
			mcmc_numstep = mcmc_params_tmp$mcmc_numstep
			cat(paste0("mcmc_numstep = ", mcmc_numstep, " for isignal = ", isignal, "\n"))
			Xcopies = generate_copies(M,thetah,mcmc_params,mcmc_numstep,parameters,experiment)
			pval_aCSS[isignal] = 
				(1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+M)
		}else{
			pval_aCSS[isignal] = 1
		}
		
	}
	
	pvals = cbind(pval_oracle,pval_aCSS,pval_score)
	colnames(pvals) = c('oracle','aCSS','score')
	rownames(pvals) = parameters$signal
	pvals
}


# Checks whether a parameter value theta is a strict second-order stationary point 
# of the perturbed likelihood determined by data x, noise parameter sigma, and Gaussian noise w
check_SSOSP = function(x,w,theta,parameters,experiment){
	# Checks the minimum eigenvalue of the Hessian is positive 
	# and that the gradient is within a small tolerance of zero
	(min(eigen(hessian(x,theta,parameters,experiment))$values)>0) & 
		(sum((gradient(x,theta,parameters,experiment)+parameters$sigma*w)^2)<=parameters$tol)
}



# Generates M copies (stored as entries in a list object) of X 
# using the hub-and-spoke sampler with each spoke just applying the mcmc function
generate_copies = function(M,thetah,mcmc_params,mcmc_numstep,parameters,experiment){
	Xhub = mcmc(experiment$X,thetah,mcmc_params,mcmc_numstep,parameters,experiment)$x
	Xcopies = list()
	for(m in 1:M){
		Xcopies[[m]] = mcmc(Xhub,thetah,mcmc_params,mcmc_numstep,parameters,experiment)$x
	}
	Xcopies
}


# Runs a Metropolis-Hastings Markov chain mcmc_numstep steps starting from state x_init,
# and returns the final state and how many times the MH sampler accepted the proposal
mcmc = function(x_init,theta,mcmc_params,mcmc_numstep,parameters,experiment){
	x = x_init
	num_moves = 0
	for(step in 1:mcmc_numstep){
	 	x_new = mcmc_proposal_draw(x,theta,mcmc_params,parameters,experiment)
	 	acceptance_prob = mcmc_acceptance_prob(x,x_new,theta,mcmc_params,parameters,experiment)
	 	accept = (runif(1) <= acceptance_prob)
	 	if(accept){
	 		x = x_new
	 		num_moves = num_moves + 1
	 	}
	}
	output = list()
	output$x = x
	output$num_moves = num_moves
	output
}

# Generates true_mcmc_params_ntrial fake X's using the input theta as the parameter value
# and for each one, estimates theta and uses the estimated theta (if SSOSP) 
# to run a single mcmc step for each of the mcmc tuning parameter values in try_mcmc_params;
# applies the example-specific choose_mcmc_params_rule function to outputs of the SSOSP-successful runs
choose_mcmc_params = function(theta,parameters,experiment){
	num_accept = rep(0,length(parameters$try_mcmc_params))
	num_trial = 0
	x_list = x_new_list = list()
	for(iparam in 1:length(parameters$try_mcmc_params)){
		x_new_list[[iparam]] = list()
	}
	for(itrial in 1:parameters$try_mcmc_params_ntrial){
		x = generate_X_null(theta,parameters,experiment)
		w = rnorm(parameters$d)/sqrt(parameters$d)
		thetah = thetahat(x,w,parameters,experiment)
		if(check_SSOSP(x,w,thetah,parameters,experiment)){
			num_trial = num_trial + 1
			x_list[[num_trial]] = x
			for(iparam in 1:length(parameters$try_mcmc_params)){
				out = mcmc(x,thetah,parameters$try_mcmc_params[[iparam]],1,parameters,experiment)
				num_accept[iparam] = num_accept[iparam] + out$num_moves
				x_new_list[[iparam]][[num_trial]] = out$x
			}
		}
	}
	choose_mcmc_params_rule(x_list,x_new_list,num_accept,num_trial,parameters,experiment)
}

