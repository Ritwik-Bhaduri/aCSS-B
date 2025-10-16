library('sigmoid')
library('ggplot2')
library('dplyr')
library('tidyr')
library('matrixStats')
library('mvtnorm')
library('MASS')

prior_beta = function(beta, tau = rep(1,m), log = FALSE){
  if(length(unique(tau) == 1)){
    if(log == TRUE){
      return(sum(dnorm(beta, mean=0, sd=unique(tau), log = TRUE)))
    }else{
      return(prod(dnorm(beta, mean=0, sd=unique(tau), log = FALSE)))
    }
  }else{
    if(log == TRUE){
      return(sum(sapply(1:m, function(j) dnorm(beta[j], mean=0, sd=tau[j], log = TRUE))))
    }else{
      return(prod(sapply(1:m, function(j) dnorm(beta[j], mean=0, sd=tau[j]))))
    }
  }
}

likelihood = function(beta, Y, X, log = FALSE){ ## fast
  beta = as.matrix(beta)
  if(log == TRUE){
    return(Y %*% X %*% beta - colSums(log(1+exp(X %*% beta))))
  }else{
    # only works when beta is a vector
    return(exp(sum(Y * X %*% beta)) / prod(1+exp(X %*% beta)))
  }
}
run_one_trial = function(M,seed,example = 'logistic', print_progress = FALSE, parameters = NULL, bayes_acss_parameters = NULL){
  set.seed(seed)
  source(paste0('../aCSS code/', example,'_source.R'))
  if(is.null(parameters)) parameters = generate_parameters()
  nsignal = length(parameters$signal)
  pval_BaCSS = rep(0,nsignal)
  pval_oracle = rep(0, nsignal)
  for(isignal in 1:nsignal){
    experiment = generate_experiment(isignal,parameters)
    Tfun_ = function(x){Tfun(x,parameters,experiment)}
    
    if(is.null(bayes_acss_parameters)){
      bayes_acss_parameters = list(tau = rep(1, parameters$d), burnin = 500, thinning = 10,
                                   B = 250, M = M, n = parameters$example$n)
    }
    
    n_samples_for_beta_posterior = bayes_acss_parameters$burnin + bayes_acss_parameters$thinning * bayes_acss_parameters$B
    
    # generate posterior samples for beta
    beta_posterior_samples = matrix(NA, nrow = parameters$d, ncol = n_samples_for_beta_posterior)
    beta_posterior_samples[,1] = current_beta = rep(0,parameters$d)
    
    proposal_parameters = function(Y, X, tau, starting_point = rep(0, parameters$d)){
      # precomputations
      Z = t(X) %*% Y
      const = - sum(log(sqrt(2*pi)*tau))
      
      neg_psi_beta = function(beta) - (sum(-0.5*(beta/tau)^2) + const + Y %*% X %*% beta - colSums(log(1+exp(X %*% beta))))
      neg_psi_grad = function(beta) - Z + t(X) %*% as.vector(exp(X %*% beta) / (1+exp(X %*% beta))) + 1/bayes_acss_parameters$tau^2 * beta
      optimization = optim(par = starting_point, fn = neg_psi_beta, gr = neg_psi_grad,
                           method = "L-BFGS-B", hessian = FALSE, lower = rep(-1,5), upper = rep(4, 5))
      # method = "BFGS", hessian = FALSE)
      beta_hat = optimization$par
      hessian = t(X) %*% diag(as.vector(exp(X %*% beta_hat) / (1+exp(X %*% beta_hat))^2)) %*% X + diag(1/tau^2)
      hessian_inv = solve(hessian)
      psi_beta_hat = - neg_psi_beta(beta_hat)
      return(list("mean" = beta_hat, "var" = hessian_inv, "psi_beta_hat" = psi_beta_hat, "hessian" = hessian))
    }
    
    acceptance_prob = c()
    for(i in 2:n_samples_for_beta_posterior){
      prop_params = proposal_parameters(Y=experiment$X, X=experiment$example$Z, tau = bayes_acss_parameters$tau)
      proposed_beta = mvrnorm(n = 1, mu = prop_params[["mean"]], Sigma = prop_params[["var"]])
      log_likelihood_diff = likelihood(proposed_beta, Y=experiment$X, X=experiment$example$Z, log = TRUE) -
        likelihood(current_beta, Y=experiment$X, X=experiment$example$Z, log = TRUE)
      log_prior_diff = prior_beta(proposed_beta, tau = bayes_acss_parameters$tau, log = TRUE) - 
        prior_beta(current_beta, tau = bayes_acss_parameters$tau, log = TRUE)
      proposal_diff = dmvnorm(current_beta, mean = prop_params[["mean"]], sigma = prop_params[["var"]], log = TRUE) - 
        dmvnorm(proposed_beta, mean = prop_params[["mean"]], sigma = prop_params[["var"]], log = TRUE)
      acceptance_prob[i] = min(1, exp(log_likelihood_diff + log_prior_diff + proposal_diff))
      if(runif(1) <= acceptance_prob[i]) current_beta = proposed_beta
      beta_posterior_samples[,i] = current_beta
    }
    
    burnin = bayes_acss_parameters$burnin
    thinning = bayes_acss_parameters$thinning
    retained_samples = burnin + (1:bayes_acss_parameters$B)*thinning
    beta_posterior_samples = beta_posterior_samples[,retained_samples]
    beta_plugin = rowMeans(beta_posterior_samples)
    
    # generate X tilde
    
    gibbs_conditional_ratio = function(Y, X, beta_posterior_samples, beta_plugin, i){ #P(Y_i=1 | everything)/P(Y_i=0 | everything)
      plugin_term = sum(X[i,]*beta_plugin)
      first_term =  sum(X[i,]*rowSums(beta_posterior_samples))
      
      # second term calculations
      Y_0 = Y_1 = Y
      Y_0[i]=0
      Y_1[i]=1
      
      prop_params_0 = proposal_parameters(Y=Y_0, X=experiment$example$Z, tau = bayes_acss_parameters$tau)
      prop_params_1 = proposal_parameters(Y=Y_1, X=experiment$example$Z, tau = bayes_acss_parameters$tau)
      log_f_Y_0 = prop_params_0$psi_beta_hat + parameters$d/2 * log(2*pi) - 1/2 * as.numeric(determinant.matrix(prop_params_0$hessian, logarithm = TRUE)$modulus)
      log_f_Y_1 = prop_params_1$psi_beta_hat + parameters$d/2 * log(2*pi) - 1/2 * as.numeric(determinant.matrix(prop_params_1$hessian, logarithm = TRUE)$modulus)
      second_term = log_f_Y_0-log_f_Y_1
      
      ans = first_term + (bayes_acss_parameters$B-1) * second_term
      return(exp(ans))
    }
    
    m_star = sample(0:bayes_acss_parameters$M, 1)
    m_star = min(m_star, M-m_star)
    
    if(m_star == 0){
      Xcopies = rep(list(NA), bayes_acss_parameters$M)
      X_current = experiment$X
      prob_X_matrix <- matrix(NA, nrow = bayes_acss_parameters$M, ncol = bayes_acss_parameters$n)
      for(iter in 1:(bayes_acss_parameters$L * bayes_acss_parameters$M)) {
        for(i in 1:bayes_acss_parameters$n) {
          p_x_ratio = gibbs_conditional_ratio(Y=X_current, X=experiment$example$Z, beta_posterior_samples, beta_plugin, i)#P(Y_i=1 | everything)/P(Y_i=0 | everything)
          p_x = p_x_ratio/(1+p_x_ratio)
          if(p_x_ratio == Inf) p_x = 1
          if(p_x_ratio == 0) p_x = 0
          
          X_current[i] <- rbinom(1, size = 1, prob = p_x)
          prob_X_matrix[iter, i] <- p_x
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies[[iter / bayes_acss_parameters$L]] <- X_current
      }
    }else{
      Xcopies1 = rep(list(NA), m_star)
      Xcopies2 = rep(list(NA), bayes_acss_parameters$M - m_star)
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * m_star)) {
        for(i in bayes_acss_parameters$n:1) {
          p_x_ratio = gibbs_conditional_ratio(Y=X_current, X=experiment$example$Z, beta_posterior_samples, beta_plugin, i)#P(Y_i=1 | everything)/P(Y_i=0 | everything)
          p_x = p_x_ratio/(1+p_x_ratio)
          if(p_x_ratio == Inf) p_x = 1
          if(p_x_ratio == 0) p_x = 0
          
          X_current[i] <- rbinom(1, size = 1, prob = p_x)
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies1[[iter / bayes_acss_parameters$L]] <- X_current
      }
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * (bayes_acss_parameters$M - m_star))) {
        for(i in 1:bayes_acss_parameters$n) {
          p_x_ratio = gibbs_conditional_ratio(Y=X_current, X=experiment$example$Z, beta_posterior_samples, beta_plugin, i)#P(Y_i=1 | everything)/P(Y_i=0 | everything)
          p_x = p_x_ratio/(1+p_x_ratio)
          if(p_x_ratio == Inf) p_x = 1
          if(p_x_ratio == 0) p_x = 0
          
          X_current[i] <- rbinom(1, size = 1, prob = p_x)
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies2[[iter / bayes_acss_parameters$L]] <- X_current
      }
      Xcopies = c(Xcopies1, Xcopies2)
    }
    
    Xcopies_matrix = matrix(unlist(Xcopies), nrow = bayes_acss_parameters$M, byrow = TRUE)
    
    pval_BaCSS[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)
    if(print_progress == TRUE) cat("\rProgress: Iteration", isignal, "of", nsignal, sep = " ")
  
    # oracle method
    Xcopies = list()
    for(m in 1:bayes_acss_parameters$M){
      Xcopies[[m]] = generate_X_null(parameters$theta0, parameters, experiment)
    }
    pval_oracle[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)
    
    
  }
  
  pvals = cbind(pval_oracle, pval_BaCSS)
  colnames(pvals) = c('oracle', 'Bayesian_aCSS')
  rownames(pvals) = parameters$signal
  pvals
}
