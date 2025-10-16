


generate_parameters = function(){
	parameters = list()
	parameters$d = 3
	parameters$theta0 = c(1,-.5,2)
	parameters$sigma = 1
	parameters$tol = 1e-6
	parameters$signal = seq(from=2,to=10,by=2)
	parameters$try_mcmc_params = lapply(seq(from=5,to=25,by=5),list) 
	parameters$try_mcmc_params_ntrial = 50
	
	# example specific parameters
	parameters$example = list()
	parameters$example$df0 = 2
	parameters$example$n = 100
	parameters$example$k = 2
	parameters$example$numstep_mult = 2
	parameters$example$numstep_max = 500
	parameters$example$acceptance_prob_min = 0.2

	parameters
}

generate_experiment = function(isignal,parameters){
	experiment = list()
	df = parameters$signal[isignal]
	theta_ = matrix(0,parameters$example$k,parameters$example$k)
	theta_[!lower.tri(theta_)] = parameters$theta0
	theta_ = theta_ + t(theta_) - diag(diag(theta_))

	experiment$X = matrix(rnorm(parameters$example$n*parameters$example$k),parameters$example$n,parameters$example$k)%*%chol(solve(theta_)) / sqrt(rchisq(parameters$example$n,df)/df)

	# example specific random variables
	experiment$example = list()

	experiment
}

Tfun = function(x,parameters,experiment){
	c(qchisq(1-score_test(x,parameters,experiment),1))
}
generate_X_null = function(theta,parameters,experiment){
	theta_ = matrix(0,parameters$example$k,parameters$example$k)
	theta_[!lower.tri(theta_)] = theta
	theta_ = theta_ + t(theta_) - diag(diag(theta_))
	matrix(rnorm(parameters$example$n*parameters$example$k),parameters$example$n,parameters$example$k)%*%chol(solve(theta_)) / sqrt(rchisq(parameters$example$n,parameters$example$df0)/parameters$example$df0)
}
negloglik = function(x,theta,parameters,experiment){
	theta_ = matrix(0,parameters$example$k,parameters$example$k)
	theta_[!lower.tri(theta_)] = theta
	theta_ = theta_ + t(theta_) - diag(diag(theta_))

	-parameters$example$n/2*log(det(theta_))+(parameters$example$df0+parameters$example$k)/2*sum(log(parameters$example$df0+diag(x%*%theta_%*%t(x))))
}
gradient = function(x,theta,parameters,experiment){
	theta_ = matrix(0,parameters$example$k,parameters$example$k)
	theta_[!lower.tri(theta_)] = theta
	theta_ = theta_ + t(theta_) - diag(diag(theta_))
	grad_ = -parameters$example$n/2*solve(theta_) + (parameters$example$df0+parameters$example$k)/2*t(x)%*%diag(1/(parameters$example$df0+diag(x%*%theta_%*%t(x))))%*%x
	(grad_*(upper.tri(grad_)+(!lower.tri(grad_))))[!lower.tri(grad_)]
}
hessian = function(x,theta,parameters,experiment){
	theta_ = matrix(0,parameters$example$k,parameters$example$k)
	theta_[!lower.tri(theta_)] = theta
	theta_ = theta_ + t(theta_) - diag(diag(theta_))
	hess_ = array(0,c(parameters$example$k,parameters$example$k,parameters$example$k,parameters$example$k))
	for(i1 in 1:parameters$example$k){for(j1 in 1:parameters$example$k){for(i2 in 1:parameters$example$k){for(j2 in 1:parameters$example$k){
		hess_[i1,j1,i2,j2] = parameters$example$n/2*solve(theta_)[i1,i2]*solve(theta_)[j1,j2] - 
			(parameters$example$df0+parameters$example$k)/2 * sum(x[,i1]*x[,i2]*x[,j1]*x[,j2] / (parameters$example$df0 + diag(x%*%theta_%*%t(x)))^2)
	}}}}
	hess = matrix(0,parameters$d,parameters$d)
	ind1 = 0
	for(j1 in 1:parameters$example$k){for(i1 in 1:j1){
		ind1 = ind1+1
		ind2 = 0
		for(j2 in 1:parameters$example$k){for(i2 in 1:j2){
			ind2 = ind2+1
			hess[ind1,ind2] = hess_[i1,j1,i2,j2] * (1+(i1!=j1)) * (1+(i2!=j2))
		}}
	}}
	hess
}
thetahat = function(x,w,parameters,experiment){
	sd_est = (apply(x,2,function(t){quantile(t,0.75)}) - 
			   apply(x,2,function(t){quantile(t,0.25)}))/(qt(.75,parameters$example$df0)-qt(.25,parameters$example$df0))
	kendalltau = matrix(0,parameters$example$k,parameters$example$k)
	for(i in 1:(parameters$example$n-1)){
		for(j in (i+1):parameters$example$n){
			kendalltau = kendalltau + sign(outer(x[i,]-x[j,],x[i,]-x[j,]))
		}
	}
	kendalltau = kendalltau / choose(parameters$example$n,2)
	theta_init = solve(diag(sd_est)%*%sin(pi/2*kendalltau)%*%diag(sd_est))
	if(min(eigen(theta_init)$values)<=0){
		theta_init = diag(1/sd_est^2)
	}
	theta_init = theta_init[!lower.tri(theta_init)]
	
	thetamat = function(theta){
		theta_ = matrix(0,parameters$example$k,parameters$example$k)
		theta_[!lower.tri(theta_)] = theta
		theta_ = theta_ + t(theta_) - diag(diag(theta_))
		theta_
	}
	fn = function(theta){negloglik(x,theta,parameters,experiment) + parameters$sigma*sum(w*theta)}
	gr = function(theta){gradient(x,theta,parameters,experiment)+parameters$sigma*w}
	err = 1;iter = 0
	theta = theta_init
	stepsize = 1
	while((err > 1e-8)&(iter<5000)){
		theta_ = theta - stepsize * gr(theta)
		if((min(eigen(thetamat(theta_))$values)<=0)){
			stepsize = stepsize / 2
		}else{
			if((fn(theta_)>=fn(theta))){
				stepsize = stepsize / 2
			}else{
				theta = theta_
				err = sum((gr(theta)^2))
				stepsize = stepsize*1.1				
			}
		}
		iter = iter+1
	}
	theta
}

mcmc_proposal_draw = function(x,theta,mcmc_params,parameters,experiment){
	S = sample(parameters$example$n,mcmc_params[[1]])
	x_new = x
	x_new[S,] = generate_X_null(theta,parameters,experiment)[S,]
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

score_test = function(x,parameters,experiment){
	thetamat = function(theta){
                theta_ = matrix(0,parameters$example$k,parameters$example$k)
                theta_[!lower.tri(theta_)] = theta
                theta_ = theta_ + t(theta_) - diag(diag(theta_))
                theta_
        }
	thetahat = thetamat(thetahat(x,0,parameters,experiment))
	s = diag(x%*%thetahat%*%t(x))
	score = (1/2)*sum(digamma((parameters$example$df0+parameters$example$k)/2) - digamma(parameters$example$df0/2) + log(parameters$example$df0) + 1 - log(parameters$example$df0+s) - (parameters$example$df0+parameters$example$k)/(parameters$example$df0+s))
	ofi11 = -(1/2)*sum(trigamma((parameters$example$df0+parameters$example$k)/2)/2 - trigamma(parameters$example$df0/2)/2 + 1/parameters$example$df0 - 1/(parameters$example$df0+s) + (parameters$example$df0+parameters$example$k)/((parameters$example$df0+s)^2) - 1/(parameters$example$df0+s))
	ofi12 = matrix(0,parameters$example$k,parameters$example$k)
	for(i in 1:parameters$example$n){
  	        ofi12 = ofi12 - (1/2)*((parameters$example$df0+parameters$example$k)/(parameters$example$df0+s[i])-1)*t(x[i,,drop=F])%*%x[i,,drop=F]/(parameters$example$df0+s[i])
	}
	ofi12 = ofi12[!lower.tri(ofi12)]
	ofi22 = matrix(0,parameters$example$k*(parameters$example$k+1)/2,parameters$example$k*(parameters$example$k+1)/2)
	for(i in 1:parameters$example$k){
	        for(j in i:parameters$example$k){
	      	        for(k in 1:parameters$example$k){
				for(l in k:parameters$example$k){
					ofi22[i+j*(j-1)/2,k+l*(l-1)/2] = (parameters$example$n/2)*(solve(thetahat)%*%diag(parameters$example$k)[,k]%*%diag(parameters$example$k)[l,]%*%solve(thetahat))[i,j] - ((parameters$example$df0+parameters$example$k)/2)*sum(x[,i]*x[,j]*x[,k]*x[,l]/((parameters$example$df0+s)^2))
				}
			}
		}
	}
	ofi = ofi11 - t(ofi12)%*%solve(ofi22)%*%ofi12
	teststat = score^2/ofi
	pval = 1-pchisq(teststat,1)
}


