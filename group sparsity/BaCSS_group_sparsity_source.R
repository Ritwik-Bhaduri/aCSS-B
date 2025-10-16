library(matrixStats)
library(MASS)
get_posterior_samples = function(X, Z, parameters, bayes_acss_parameters, precomputations){
  A = precomputations$A
  A_inv = precomputations$A_inv
  
  b = C = vector("list", parameters$G)
  log_w = numeric(parameters$G)
  post_means = vector("list", parameters$G)
  post_sigma = vector("list", parameters$G)
  for(g in 1:parameters$G){
    Z_g = Z[ ,parameters$grp_indices[[g]]]
    b[[g]] = t(Z_g)%*%X/parameters$sigma^2 + bayes_acss_parameters$mu/bayes_acss_parameters$tau^2
    C[[g]] = sum(X^2)/parameters$sigma^2 + (parameters$grp_sizes[g]*bayes_acss_parameters$mu^2)/bayes_acss_parameters$tau^2
    post_sigma[[g]] = A_inv[[g]]
    post_means[[g]] = post_sigma[[g]]%*%b[[g]]
  }
  for (g in 1:parameters$G){
    log_w[g] = -parameters$grp_sizes[g]*log(bayes_acss_parameters$tau) - 
      (1/2)*as.numeric(determinant(A[[g]], logarithm = TRUE)$modulus) + 
      (1/2)*t(b[[g]])%*%post_means[[g]] - (1/2)*C[[g]]
  }
  norm_const = logSumExp(log_w)
  log_w = log_w - norm_const
  post_samples = matrix(0,nrow = bayes_acss_parameters$B,ncol = parameters$p)
  w = exp(log_w)
  for (i in 1:bayes_acss_parameters$B){
    component = sample(parameters$G,size = 1, prob = w)
    post_samples[i,parameters$grp_indices[[component]]] = mvrnorm(1, post_means[[component]], post_sigma[[component]])
  }
  return(post_samples)
}

log_sampling_density = function(X, Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations){
  A = precomputations$A
  A_inv = precomputations$A_inv
  
  b = C = log_D = vector("list", parameters$G)
  for (g in 1:parameters$G){
    Z_g = Z[ ,parameters$grp_indices[[g]]]
    b[[g]] = t(Z_g)%*%X/parameters$sigma^2 + bayes_acss_parameters$mu/bayes_acss_parameters$tau^2
    C[[g]] = sum(X^2)/parameters$sigma^2 + (parameters$grp_sizes[g]*bayes_acss_parameters$mu^2)/bayes_acss_parameters$tau^2
    log_D[[g]] = -parameters$grp_sizes[g]*log(bayes_acss_parameters$tau) -1/2 * as.numeric(determinant(A[[g]], logarithm = TRUE)$modulus) +
      1/2 * t(b[[g]]) %*% A_inv[[g]] %*% b[[g]] - 1/2 * C[[g]]
  }
  log_D = unlist(log_D)
  t1 = -1/(2 * parameters$sigma^2) * sum((X[i] - Z[i, ] %*% t(theta_posterior_samples))^2)
  t2 = -(bayes_acss_parameters$B-1) * (-log(parameters$G) + logSumExp(log_D))
  return(t1+t2)
}

d_log_sampling_density = function(X, Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations){
  A = precomputations$A
  A_inv = precomputations$A_inv
  
  b = C = vector("list", parameters$G)
  log_D = D_prime_by_D = numeric(parameters$G)
  for (g in 1:parameters$G){
    Z_g = Z[ ,parameters$grp_indices[[g]]]
    b[[g]] = t(Z_g)%*%X/parameters$sigma^2 + bayes_acss_parameters$mu/bayes_acss_parameters$tau^2
    C[[g]] = sum(X^2)/parameters$sigma^2 + (parameters$grp_sizes[g]*bayes_acss_parameters$mu^2)/bayes_acss_parameters$tau^2
    
    log_D[[g]] = -parameters$grp_sizes[g]*log(bayes_acss_parameters$tau) -1/2 * as.numeric(determinant(A[[g]], logarithm = TRUE)$modulus) +
      1/2 * t(b[[g]]) %*% A_inv[[g]] %*% b[[g]] - 1/2 * C[[g]]
    D_prime_by_D[[g]] = 1/parameters$sigma^2 * (Z_g[i,] %*% A_inv[[g]] %*% b[[g]] - X[i])
  }
  t1 = -1/parameters$sigma^2 * sum(X[i] - Z[i, ] %*% t(theta_posterior_samples))
  
  D_by_sum_D = exp(log_D - logSumExp(log_D))
  t2 = -(bayes_acss_parameters$B-1) * sum(D_prime_by_D * D_by_sum_D)
  
  return(t1 + t2)
}

d2_log_sampling_density = function(X, Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations) {
  A = precomputations$A
  A_inv = precomputations$A_inv
  
  b = C = vector("list", parameters$G)
  log_D = D_prime_by_D = D_double_prime_by_D = numeric(parameters$G)
  
  for (g in 1:parameters$G) {
    Z_g = Z[, parameters$grp_indices[[g]]]
    b[[g]] = t(Z_g) %*% X / parameters$sigma^2 + bayes_acss_parameters$mu / bayes_acss_parameters$tau^2
    C[[g]] = sum(X^2) / parameters$sigma^2 + (parameters$grp_sizes[g] * bayes_acss_parameters$mu^2) / bayes_acss_parameters$tau^2
    
    log_D[g] = -parameters$grp_sizes[g]*log(bayes_acss_parameters$tau) -1/2 * as.numeric(determinant(A[[g]], logarithm = TRUE)$modulus) +
      1/2 * t(b[[g]]) %*% A_inv[[g]] %*% b[[g]] - 1/2 * C[[g]]
    D_prime_by_D[g] = 1/parameters$sigma^2 * (Z_g[i,] %*% A_inv[[g]] %*% b[[g]] - X[i])
    
    t1 = 1 / parameters$sigma^2 * (1 / parameters$sigma^2 * Z_g[i,] %*% A_inv[[g]] %*% as.matrix(Z_g[i,]) - 1)
    t2 = 1 / parameters$sigma^2 * (Z_g[i,] %*% A_inv[[g]] %*% b[[g]] - X[i]) * D_prime_by_D[g]
    D_double_prime_by_D[g] = t1 + t2
  }
  D_by_sum_D = exp(log_D - logSumExp(log_D))
  
  d2_log_density = - bayes_acss_parameters$B / parameters$sigma^2 - 
    (bayes_acss_parameters$B-1) * (sum(D_double_prime_by_D*D_by_sum_D) - sum(D_prime_by_D * D_by_sum_D)^2)
  
  return(d2_log_density)
}

get_laplace_params = function(X, Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations){
  optim_function <- function(X_i){
    X_vec = X
    X_vec[i] = X_i
    return(-log_sampling_density(X_vec, Z, i, theta_posterior_samples, parameters=parameters, bayes_acss_parameters=bayes_acss_parameters, precomputations))
  }
  
  optim_gradient <- function(X_i){
    X_vec = X
    X_vec[i] = X_i
    return(-d_log_sampling_density(X_vec, Z, i, theta_posterior_samples, parameters=parameters, bayes_acss_parameters=bayes_acss_parameters, precomputations))
  }
  
  optim_output <- optim(X[i], fn = optim_function, gr = optim_gradient, method = "BFGS", hessian = FALSE)$par
  
  X_optim = X
  X_optim[i] = optim_output
  var = - 1/d2_log_sampling_density(X_optim, Z, i, theta_posterior_samples, parameters=parameters, bayes_acss_parameters=bayes_acss_parameters, precomputations)
  output = list("mean" = optim_output, "var" = var)
  return(output)
}


run_one_trial = function(seed,example = 'group_lasso', print_progress = FALSE, parameters = NULL, bayes_acss_parameters = NULL){
  set.seed(seed)
  source(paste0(example,'_source.R'))
  if(is.null(parameters)) parameters = generate_parameters()
  nsignal = length(parameters$signal)
  pval_BaCSS = rep(0, nsignal)
  pval_oracle = rep(0, nsignal)
  for(isignal in 1:nsignal){
    experiment = generate_experiment(isignal,parameters)
    Tfun_ = function(x){Tfun(x,parameters,experiment)}
    
    if(is.null(bayes_acss_parameters)){
      bayes_acss_parameters = list(mu = 0, tau = 1, # prior parameters
                                   M=250, B = 25, n = parameters$n, burnin_X = 25)
    }
    
    # precomputation
    A = A_inv = vector("list", parameters$G)
    for(g in 1:parameters$G){
      Z_g = experiment$example$Z[, parameters$grp_indices[[g]]]
      A[[g]] = t(Z_g) %*% Z_g / parameters$sigma^2 + diag(parameters$grp_sizes[g]) / bayes_acss_parameters$tau^2
      A_inv[[g]] = solve(A[[g]])
    }
    precomputations= list("A" = A, "A_inv" = A_inv)
    
    theta_posterior_samples = get_posterior_samples(X=experiment$X, Z=experiment$example$Z, 
                                                    parameters, bayes_acss_parameters, precomputations)
    
    m_star = sample(0:bayes_acss_parameters$M, 1)
    m_star = min(m_star, bayes_acss_parameters$M - m_star)
    
    if(m_star == 0){
      Xcopies = rep(list(NA), bayes_acss_parameters$M)
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$M * bayes_acss_parameters$L)) {
        for(i in 1:parameters$n) {
          laplace_params = get_laplace_params(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          likelihood_diff = log_sampling_density(X_proposed, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations) -
            log_sampling_density(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
          
          acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
          
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) 
            cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",", "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies[[iter / bayes_acss_parameters$L]] <- X_current
      }
    }else{
      Xcopies1 = rep(list(NA), m_star)
      Xcopies2 = rep(list(NA), bayes_acss_parameters$M - m_star)
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * m_star)) {
        for(i in parameters$n:1) {
          laplace_params = get_laplace_params(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          likelihood_diff = log_sampling_density(X_proposed, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations) -
            log_sampling_density(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
          
          acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
          
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) 
            cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",", "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies1[[iter / bayes_acss_parameters$L]] <- X_current
      }
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * (bayes_acss_parameters$M - m_star))) {
        for(i in 1:parameters$n) {
          laplace_params = get_laplace_params(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          likelihood_diff = log_sampling_density(X_proposed, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations) -
            log_sampling_density(X_current, experiment$example$Z, i, theta_posterior_samples, parameters, bayes_acss_parameters, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
          
          acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
          # print(acceptance_prob)
          
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) 
            cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",", "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies2[[iter / bayes_acss_parameters$L]] <- X_current
      }
      Xcopies = c(Xcopies1, Xcopies2)
    }
    
    pval_BaCSS[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)
    if(print_progress == TRUE) cat("\n\rProgress: Iteration", isignal, "of", nsignal, sep = " ")
    
    # oracle method
    Xcopies = list()
    for(m in 1:bayes_acss_parameters$M){
      Xcopies[[m]] = generate_X_null(parameters, experiment)$X
    }
    pval_oracle[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)
    
  }
  
  pvals = cbind(pval_oracle, pval_BaCSS)
  colnames(pvals) = c('oracle', 'Bayesian_aCSS')
  rownames(pvals) = parameters$signal
  pvals
}
