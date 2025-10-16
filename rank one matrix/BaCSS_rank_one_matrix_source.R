get_posterior_samples <- function(n_samples, experiment, prior_mean, prior_sd, parameters){
  X = experiment$X
  svd_X = svd(X)
  U_current <- sqrt(svd_X$d[1])*svd_X$u[,1]
  V_current <- sqrt(svd_X$d[1])*svd_X$v[,1]
  post_samples <- list("U" = matrix(0,n_samples,nrow(X)), "V" = matrix(0,n_samples,ncol(X)))
  for (i in 1:n_samples){
    U_var <- 1/(sum(V_current^2)/parameters$sigma^2 + 1/prior_sd[1]^2)
    U_mean <- U_var*(X%*%V_current/parameters$sigma^2 + prior_mean[1]/prior_sd[1]^2)
    U_current = rnorm(nrow(X), 0, sd = sqrt(U_var)) + U_mean
    post_samples$U[i,] = U_current
    
    V_var <- 1/(sum(U_current^2)/parameters$sigma^2 + 1/prior_sd[2]^2)
    V_mean <- V_var*(t(X)%*%U_current/parameters$sigma^2 + prior_mean[2]/prior_sd[2]^2)
    V_current = rnorm(ncol(X), 0, sd = sqrt(V_var)) + V_mean
    post_samples$V[i,] = V_current
  }
  return(post_samples)
}

log_posterior_t_given_X <- function(X, t, bayes_acss_parameters, parameters) {
  n <- parameters$example$n
  d2 <- svd(X)$d^2
  p  <- length(d2)
  W   <- exp(t)
  S   <- sum(W)
  U   <- sum(W * d2)
  den <- 1 + 4 * S
  
  # Likelihood
  log_like  <- -(n/2) * log(den) - 2 * sum(X * X) + 8 * U / den
  
  # Prior on t = log(W) induced by prior on W = V^2
  log_prior <- sum(-0.5 * log(2*pi) + 0.5 * t - 0.5 * W)
  
  log_like + log_prior
}

# Gradient wrt t:  d/dt log p(t | X)
d_log_posterior_t_given_X <- function(X, t, bayes_acss_parameters, parameters) {
  n      <- parameters$example$n
  d2 <- svd(X)$d^2
  W    <- exp(t)
  S    <- sum(W)
  U    <- sum(W * d2)
  den  <- 1 + 4 * S
  den2 <- den * den
  
  # gradient of log-likelihood wrt W
  grad_ll_W <- rep(-(n/2) * 4 / den, n) + 8 * ((den * d2 - 4 * U) / den2)
  
  # Jacobian of W -> log(W) = t
  grad_t <- (grad_ll_W * W) + (0.5 - 0.5 * W)
  grad_t
}

# log determinant of the negative Hessian wrt t
log_det_hessian_log_posterior <- function(X, t, bayes_acss_parameters, parameters, jitter = 1e-10) {
  n      <- parameters$example$n
  d2 <- (svd(X, nu = 0, nv = 0)$d)^2
  W    <- exp(t)
  S    <- sum(W)
  U    <- sum(W * d2)
  den  <- 1 + 4 * S
  den2 <- den * den
  den3 <- den2 * den
  one <- rep(1.0, n)
  
  # Hessian of log-likelihood wrt W
  A <- ((n/2) * 16 / den2) * (one %o% one)
  B <- 8 * (-4 * (outer(d2, one) + outer(one, d2)) / den2 + (8 * U / den3) * (one %o% one))
  
  H_like_W <- A + B
  
  # Transform to t-space: H_like_T = J^T H_like_W J + diag( (dL/dW) **Hadamard product** W )
  grad_ll_W <- rep(-(n/2) * 4 / den, n) + 8 * ((den * d2 - 4 * U) / den2)
  
  JW <- diag(W, n, n)
  H_like_T <- JW %*% H_like_W %*% JW + diag(grad_ll_W * W, n, n)
  
  # Add prior(t) Hessian: d^2/dt_i^2 [ 0.5 t_i - 0.5 e^{t_i} ] = -0.5 e^{t_i}
  H_T <- H_like_T + diag(-0.5 * W, n, n)
  
  # Negative Hessian (should be p.d. at the mode) + tiny jitter
  Hneg <- -(H_T + t(H_T)) / 2
  diag(Hneg) <- diag(Hneg) + jitter
  
  d <- determinant(Hneg, logarithm = TRUE)
  if (d$sign <= 0 || !is.finite(d$modulus)) return(NA_real_)
  as.numeric(d$modulus)  # log det(-H)
}

log_marginal <- function(X, bayes_acss_parameters, parameters) {
  svd_X = svd(X)
  t_init <- log(svd_X$d[1]*svd_X$v[,1]^2)
  
  # objective and gradient in t-space (maximize log posterior)
  obj_t <- function(t) -log_posterior_t_given_X(X, t, bayes_acss_parameters, parameters)
  gr_t  <- function(t) -d_log_posterior_t_given_X(X, t, bayes_acss_parameters, parameters)
  
  opt <- optim(par = t_init, fn = obj_t, gr = gr_t, method = "BFGS")
  
  t_hat <- opt$par
  lp_hat <- log_posterior_t_given_X(X, t_hat, bayes_acss_parameters, parameters)
  logdet_negH <- log_det_hessian_log_posterior(X, t_hat, bayes_acss_parameters, parameters)
  
  # Laplace approximation: dimension = length(t_hat)
  lp_hat + (length(t_hat) / 2) * log(2 * pi) - 0.5 * logdet_negH
}

log_numerator_X_tilde <- function(X_ij, i, j, theta_posterior_samples, parameters){
  U_i <- theta_posterior_samples$U[,i]
  V_j <- theta_posterior_samples$V[,j]
  return(-1/(2*parameters$sigma^2)*sum((X_ij - U_i*V_j)^2))
}


log_sampling_density <- function(X, i, j, theta_posterior_samples, parameters, bayes_acss_parameters){
  log_numerator <- log_numerator_X_tilde(X[i,j], i,j , theta_posterior_samples, parameters)
  log_denominator <- log_marginal(X, bayes_acss_parameters, parameters)
  return(log_numerator - (bayes_acss_parameters$B - 1) * log_denominator)
}

get_laplace_params = function(X, i, j, theta_posterior_samples, parameters, bayes_acss_parameters){
  optim_function <- function(X_ij){
    X_mat = X
    X_mat[i,j] = X_ij
    return(-log_sampling_density(X_mat, i, j, theta_posterior_samples, parameters, bayes_acss_parameters))
  }
  optim_output <- optim(X[i,j], fn = optim_function, method = "BFGS", hessian = TRUE)
  output = list("mean" = optim_output$par, "var" = 1/optim_output$hessian)
  output$var = max(1e-2,output$var)
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
      bayes_acss_parameters = list(prior_mean = c(0, 0), prior_sd = c(1, 1), 
                                   burnin_posterior = 500,thinning_posterior = 10,
                                   M=300, B = 25, n = parameters$example$n, L = 1)
    }
    
    theta_posterior_samples = get_posterior_samples(n_samples = (bayes_acss_parameters$burnin_posterior+1) + (bayes_acss_parameters$B-1)*bayes_acss_parameters$thinning_posterior, 
                                                    experiment = experiment, prior_mean = bayes_acss_parameters$prior_mean, 
                                                    prior_sd = bayes_acss_parameters$prior_sd,
                                                    parameters = parameters)
    theta_posterior_samples$U = theta_posterior_samples$U[(bayes_acss_parameters$burnin_posterior+1) + (0:(bayes_acss_parameters$B-1))*bayes_acss_parameters$thinning_posterior, ]
    theta_posterior_samples$V = theta_posterior_samples$V[(bayes_acss_parameters$burnin_posterior+1) + (0:(bayes_acss_parameters$B-1))*bayes_acss_parameters$thinning_posterior, ]
    
    m_star = sample(0:bayes_acss_parameters$M, 1)
    m_star = min(m_star, bayes_acss_parameters$M - m_star)
    
    if(m_star == 0){
      Xcopies = rep(list(NA), bayes_acss_parameters$M)
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L*bayes_acss_parameters$M)) {
        for(i in 1:parameters$example$m) {
          for(j in 1:parameters$example$n) {
            laplace_params = get_laplace_params(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            X_ij_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
            X_proposed = X_current; X_proposed[i,j] = X_ij_proposed
            
            likelihood_diff = log_sampling_density(X_proposed, i, j, theta_posterior_samples, parameters, bayes_acss_parameters) -
              log_sampling_density(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            proposal_diff = dnorm(X_proposed[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
              dnorm(X_current[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
            
            acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
            # print(acceptance_prob)
            
            X_ij_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i,j], X_current[i,j])
            X_current[i,j] = X_ij_retained
            if(print_progress == TRUE) 
              cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",",j, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
          }
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies[[iter / bayes_acss_parameters$L]] <- X_current
      }
    }else{
      Xcopies1 = rep(list(NA), m_star)
      Xcopies2 = rep(list(NA), bayes_acss_parameters$M - m_star)
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * m_star)) {
        for(i in 1:parameters$example$m) {
          for(j in 1:parameters$example$n) {
            laplace_params = get_laplace_params(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            X_ij_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
            X_proposed = X_current; X_proposed[i,j] = X_ij_proposed
            
            likelihood_diff = log_sampling_density(X_proposed, i, j, theta_posterior_samples, parameters, bayes_acss_parameters) -
              log_sampling_density(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            proposal_diff = dnorm(X_proposed[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
              dnorm(X_current[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
            
            acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
            # print(acceptance_prob)
            
            X_ij_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i,j], X_current[i,j])
            X_current[i,j] = X_ij_retained
            if(print_progress == TRUE) 
              cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",",j, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
          }
        }
        if(iter %% bayes_acss_parameters$L == 0) Xcopies1[[iter / bayes_acss_parameters$L]] <- X_current
      }
      
      X_current = experiment$X
      for(iter in 1:(bayes_acss_parameters$L * (bayes_acss_parameters$M - m_star))) {
        for(i in 1:parameters$example$m) {
          for(j in 1:parameters$example$n) {
            laplace_params = get_laplace_params(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            X_ij_proposed = rnorm(1, mean = laplace_params$mean, sd = sqrt(laplace_params$var))
            X_proposed = X_current; X_proposed[i,j] = X_ij_proposed
            
            likelihood_diff = log_sampling_density(X_proposed, i, j, theta_posterior_samples, parameters, bayes_acss_parameters) -
              log_sampling_density(X_current, i, j, theta_posterior_samples, parameters, bayes_acss_parameters)
            proposal_diff = dnorm(X_proposed[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE) - 
              dnorm(X_current[i,j], mean = laplace_params$mean, sd = sqrt(laplace_params$var), log=TRUE)
            
            acceptance_prob = min(1, exp(likelihood_diff - proposal_diff))
            # print(acceptance_prob)
            
            X_ij_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i,j], X_current[i,j])
            X_current[i,j] = X_ij_retained
            if(print_progress == TRUE) 
              cat("\rProgress: Iteration", iter, "of", bayes_acss_parameters$M, ", Index: ", i, ",",j, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
          }
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
      Xcopies[[m]] = generate_X_null(parameters)$X
    }
    pval_oracle[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+bayes_acss_parameters$M)
    
  }
  
  pvals = cbind(pval_oracle, pval_BaCSS)
  colnames(pvals) = c('oracle', 'Bayesian_aCSS')
  rownames(pvals) = parameters$signal
  pvals
}