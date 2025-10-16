library(matrixStats)
library(MASS)
library(segmented)
library(truncnorm)
library(mvnfast)
library(mixtools)
library(fBasics)

sourceCpp("functions.cpp")

get_posterior_beta = function(X,Z,parameters,bayes_acss_parameters,t_vec){
  k = length(t_vec) + 1 
  h_t <- cbind(1,Z,pmax(outer(Z,t_vec,"-"),0))
  beta_post_var <- solve(1/parameters$sigma^2 * t(h_t)%*%h_t + 1/bayes_acss_parameters$tau_1^2 * diag(k+1))
  beta_post_mean <- beta_post_var%*%(1/parameters$sigma^2 * t(h_t) %*% X + 
                                       bayes_acss_parameters$mu_1/bayes_acss_parameters$tau_1^2*rep(1,k+1))
  beta_post_samples <- mvrnorm(1,mu = beta_post_mean,Sigma = beta_post_var)
  return(beta_post_samples)
}

log_diff_exp <- function(x, y) {
  if(length(x) != length(y)) stop("x and y must be of the same length")
  result <- numeric(length(x))
  result[y == -Inf] <- x[y == -Inf]
  mask <- y != -Inf
  result[mask] <- x[mask] + log1p(-exp(y[mask] - x[mask]))
  return(result)
}

get_posterior_t = function(X,Z,j,parameters,bayes_acss_parameters,beta,t_vec){
  t_minus_j <- t_vec[-j]
  h_t_minus_j <- cbind(1,Z,pmax(outer(Z,t_minus_j,"-"),0))
  X_prime = X - h_t_minus_j %*% beta[-(j+2)]
  ordering = order(Z)
  Z_ordered = Z[ordering]
  X_prime_ordered = X_prime[ordering]
  X_prime_ordered_reverse = rev(X_prime_ordered)
  Z_ordered_reverse = rev(Z_ordered)
  
  n = length(X_prime)
  log_constants = rep(NA, n)
  a = rep(NA, n)
  b = rep(NA, n)
  for(i in 1:n){
    log_constants[i] = beta[j+2]/parameters$sigma^2*sum(X_prime_ordered[i:n] * Z_ordered[i:n]) -
      beta[j+2]^2/(2*parameters$sigma^2)*sum(Z_ordered[i:n]^2)
    
    a[i] = beta[j+2]^2/(2*parameters$sigma^2)*(parameters$n-i+1)+1/(2*bayes_acss_parameters$tau_2^2)
    b[i] = -beta[j+2]/parameters$sigma^2 * sum(X_prime_ordered[i:n]) +
      beta[j+2]^2/parameters$sigma^2*sum(Z_ordered[i:n]) +
      bayes_acss_parameters$mu_2/bayes_acss_parameters$tau_2^2
  }
  log_constants_2 = (1/2)*(log(pi) - log(a)) + b^2/(4*a)
  log_constants = log_constants + log_constants_2
  t_j_var = 1/(2*a)
  t_j_mean = b*t_j_var
  
  log_prob_upper <- pnorm((Z_ordered - t_j_mean)/sqrt(t_j_var),log.p = TRUE)
  log_prob_lower <- pnorm((c(-Inf,Z_ordered[-parameters$n]) - t_j_mean)/sqrt(t_j_var),log.p = TRUE)
  log_prob_last_class = pnorm(Z_ordered[parameters$n],mean = bayes_acss_parameters$mu_2,
                              sd = bayes_acss_parameters$tau_2,lower.tail = FALSE,log.p = TRUE)
  log_probs = c(log_constants + log_diff_exp(log_prob_upper,log_prob_lower),log_prob_last_class)
  log_probs_normalized = log_probs - max(log_probs)
  
  i = sample(0:parameters$n,1,prob = exp(log_probs_normalized))     ### active group
  t_j_mean = c(t_j_mean, bayes_acss_parameters$mu_2)
  t_j_var = c(t_j_var, bayes_acss_parameters$tau_2^2)
  
  if (i==0){
    t_j <- rtruncnorm(1,a = -Inf, b = Z_ordered[1], mean = t_j_mean[i+1],sd = sqrt(t_j_var[i+1]))
  }else if(i==parameters$n){
    t_j <- rtruncnorm(1,a = Z_ordered[parameters$n], b = Inf, mean = t_j_mean[i+1],sd = sqrt(t_j_var[i+1]))
  }else{
    t_j <- rtruncnorm(1,a = Z_ordered[i], b = Z_ordered[i+1],  mean = t_j_mean[i+1],sd = sqrt(t_j_var[i+1]))
  }
  return(t_j)
}

get_posterior_samples = function(X, Z, beta_init=NULL, t_init=NULL, parameters, bayes_acss_parameters){
  
  if(is.null(t_init) || is.null(beta_init)){
    initial_model <- lm(X ~ Z)
    segmented_model <- segmented(initial_model, seg.Z = ~Z, npsi = parameters$k-1)
  }
  if(is.null(t_init)){
    t_init <- segmented_model$psi[, "Est."]
  }
  if(is.null(beta_init)){
    h_t_init <- cbind(1,Z,pmax(outer(Z,t_init,"-"),0))
    spline_model <- lm(X ~ 0 + h_t_init)
    beta_init <- coef(spline_model)
  }
  
  beta_post_samples <- matrix(0,nrow = (bayes_acss_parameters$B+bayes_acss_parameters$burnin_posterior),ncol = parameters$k+1)
  t_post_samples <- matrix(0,nrow = (bayes_acss_parameters$B+bayes_acss_parameters$burnin_posterior),ncol = parameters$k-1)
  t_current = t_init
  beta_current = beta_init
  for (i in 1:(bayes_acss_parameters$B+bayes_acss_parameters$burnin_posterior)){
    beta_post_samples[i,] <- get_posterior_beta(X,Z,parameters, bayes_acss_parameters, t_current)
    beta_current = beta_post_samples[i,]
    for (j in 1:(parameters$k-1)){
      t_post_samples[i,j] = get_posterior_t(X,Z,j,parameters,bayes_acss_parameters,beta_current, t_current)
      t_current[j] = t_post_samples[i,j]
    }
    t_current = sort(t_current)
  }
  return(list("t_post_samples" = t_post_samples,"beta_post_samples" = beta_post_samples,
              "beta_init" = beta_init, "t_init" = t_init))
}

log_posterior_t_given_X = function(X, t_vec, bayes_acss_parameters, parameters,experiment){
  # log post calculations
  Z = experiment$example$Z
  h_t <- cbind(1,Z,pmax(outer(Z,t_vec,"-"),0))
  X_t_mean <- bayes_acss_parameters$mu_1*rowSums(h_t)
  
  A = (1/bayes_acss_parameters$tau_1^2)*diag(parameters$k+1)+1/parameters$sigma^2*t(h_t)%*%h_t
  A_inv = chol2inv(chol(A))
  Sigma_inv = 1/parameters$sigma^2*diag(parameters$n) - 1/parameters$sigma^4*h_t%*%A_inv%*%t(h_t)
  
  return( #- parameters$n * log(parameters$sigma) -
    -0.5*as.numeric(determinant(A * bayes_acss_parameters$tau_1^2,logarithm=TRUE)$modulus)-
      0.5*t(X-X_t_mean)%*%Sigma_inv%*%(X-X_t_mean) - 
      sum((t_vec - bayes_acss_parameters$mu_2)^2/(2*bayes_acss_parameters$tau_2^2)))  }


log_numerator_X_tilde <- function(X_i, i, theta_posterior_samples, parameters,experiment){
  Z = experiment$example$Z
  t_mat <- theta_posterior_samples$t_post_samples
  beta_mat <- theta_posterior_samples$beta_post_samples
  answer_vec = - ((X_i - rowSums(cbind(1, Z[i], pmax(Z[i] - t_mat, 0))*beta_mat))^2) / (2 * parameters$sigma^2)
  answer <- sum(answer_vec)
  return(answer)
}

log_sampling_density <- function(X, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations){
  log_numerator <- log_numerator_X_tilde(X[i], i, theta_posterior_samples, parameters, experiment)
  log_denominator <- log_marginal_cpp(X, bayes_acss_parameters, parameters,experiment, precomputations)
  return(log_numerator - (bayes_acss_parameters$B - 1) * log_denominator)
}

get_laplace_params = function(X, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations){
  optim_function <- function(X_i){
    X_vec = X
    X_vec[i] = X_i
    return(-log_sampling_density(X_vec, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations))
  }
  optim_output <- optim(X[i], fn = optim_function, method = "BFGS", hessian = TRUE)
  output = list("mean" = optim_output$par, "var" = 1/optim_output$hessian)
  return(output)
}

run_one_trial = function(seed,example = 'rank_one_matrix', print_progress = FALSE, parameters = NULL, bayes_acss_parameters = NULL){
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
      bayes_acss_parameters = list(mu_1 = 0, mu_2 = 0, tau_1 = 1,tau_2 = 1, # prior parameters
                                   burnin_posterior = 100,
                                   M=5, B = 15, n = parameters$n, L = 1, var_min = 1e-8)
    }

    theta_posterior_samples = get_posterior_samples(experiment$X, experiment$example$Z,
                                                    parameters = parameters,bayes_acss_parameters = bayes_acss_parameters)
    theta_posterior_samples$beta_post_samples = theta_posterior_samples$beta_post_samples[bayes_acss_parameters$burnin_posterior + 1:bayes_acss_parameters$B, ]
    theta_posterior_samples$t_post_samples = theta_posterior_samples$t_post_samples[bayes_acss_parameters$burnin_posterior + 1:bayes_acss_parameters$B, ]

    Z = experiment$example$Z
    Z_sorted = sort(Z)
    precomputations = list('Z_sorted' = Z_sorted)
    
    m_star = sample(0:bayes_acss_parameters$M, 1)
    m_star = min(m_star, bayes_acss_parameters$M - m_star)
    
    if(m_star == 0){
      Xcopies = rep(list(NA), bayes_acss_parameters$M)
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$M * bayes_acss_parameters$L)) {
        for(i in 1:parameters$n) {
          laplace_params = get_laplace_params(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations)
          X_i_var = ifelse(laplace_params$var >= 0, laplace_params$var, bayes_acss_parameters$var_min)
          
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(X_i_var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
  
          likelihood_diff = log_sampling_density(X_proposed, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations) -
            log_sampling_density(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE) -
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE)
  
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
          laplace_params = get_laplace_params(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations)
          X_i_var = ifelse(laplace_params$var >= 0, laplace_params$var, bayes_acss_parameters$var_min)
          
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(X_i_var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          likelihood_diff = log_sampling_density(X_proposed, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations) -
            log_sampling_density(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE) -
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE)
          
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
          laplace_params = get_laplace_params(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment,precomputations)
          X_i_var = ifelse(laplace_params$var >= 0, laplace_params$var, bayes_acss_parameters$var_min)
          
          X_i_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(X_i_var))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          likelihood_diff = log_sampling_density(X_proposed, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations) -
            log_sampling_density(X_current, i, theta_posterior_samples, parameters, bayes_acss_parameters,experiment, precomputations)
          proposal_diff = dnorm(X_proposed[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE) -
            dnorm(X_current[i], mean = laplace_params$mean, sd = sqrt(X_i_var), log=TRUE)
          
          acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
          
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
      Xcopies[[m]] = generate_X_null(parameters,experiment)$X
    }
    pval_oracle[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)

  }

  pvals = cbind(pval_oracle, pval_BaCSS)
  colnames(pvals) = c('oracle', 'Bayesian_aCSS')
  rownames(pvals) = parameters$signal
  pvals
}
