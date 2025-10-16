library('dplyr')
library(DPQ)
library("LaplacesDemon")
library(Rcpp)
sourceCpp("mixture.cpp")

prior_theta = function(theta, prior_params, log = FALSE){
  if(log == TRUE){
    prior_w = LaplacesDemon::ddirichlet(theta$w, prior_params$w$pi, log = TRUE) 
    prior_mu = sum(dnorm(theta$mu, mean = prior_params$mu$mean, sd = sqrt(theta$sigma2 / prior_params$mu$kappa), log = TRUE))
    prior_sigma2 = sum(dinvgamma(theta$sigma2, shape = prior_params$sigma2$alpha, scale = prior_params$sigma2$beta, log = TRUE))
    return(prior_w + prior_mu + prior_sigma2)
  }else{
    prior_w = LaplacesDemon::ddirichlet(theta$w, prior_params$w$pi, log = FALSE) 
    prior_mu = prod(dnorm(theta$mu, mean = prior_params$mu$mean, sd = sqrt(theta$sigma2 / prior_params$mu$kappa), log = FALSE))
    prior_sigma2 = prod(dinvgamma(theta$sigma2, shape = prior_params$sigma2$alpha, scale = prior_params$sigma2$beta, log = FALSE))
    return(prior_w * prior_mu * prior_sigma2)
  }
}

library(matrixStats)


log_marginal <- function(X, bayes_acss_parameters, parameters){
  J <- parameters$J
  if(J != 2) stop("This optimized version only works for J = 2")
  
  # softmax for J=2
  softmax2 <- function(z){
    ez <- exp(z)
    denom <- 1 + ez
    c(ez / denom, 1 / denom)
  }
  
  # objective in transformed space
  optim_function <- function(theta_t){
    mu <- theta_t[1:2]
    eta <- theta_t[3:4]
    z   <- theta_t[5]
    sigma2 <- exp(eta)
    w <- softmax2(z)
    theta_vec <- c(mu, sigma2, w)
    -log_posterior_theta_given_X(X, theta_vec, bayes_acss_parameters)
  }
  
  # gradient in transformed space
  optim_gradient <- function(theta_t){
    mu <- theta_t[1:2]
    eta <- theta_t[3:4]
    z   <- theta_t[5]
    sigma2 <- exp(eta)
    w <- softmax2(z)
    theta_vec <- c(mu, sigma2, w)
    
    g <- d_log_posterior_theta_given_X(X, theta_vec, bayes_acss_parameters)
    g_mu <- g[1:2]
    g_s2 <- g[3:4]
    g_w  <- g[5:6]
    
    s <- sum(g_w * w)  # <g_w, w>
    g_z <- w[1] * (g_w[1] - s)
    
    c(g_mu, sigma2 * g_s2, g_z) * -1
  }
  
  # initialization
  km <- kmeans(X, centers = 2)
  clusters <- km$cluster
  mu0 <- sapply(1:2, function(i){
    xi <- X[clusters == i]
    if(length(xi) == 0) mean(X) else mean(xi)
  })
  s20 <- sapply(1:2, function(i){
    xi <- X[clusters == i]
    v <- if(length(xi) <= 1) var(X) else var(xi)
    if(!is.finite(v)) var(X) else v
  }) + 1e-5
  w0 <- tabulate(clusters, nbins=2) / length(X)
  w0 <- pmax(w0, 1e-12); w0 <- w0 / sum(w0)
  z0 <- log(w0[1] / w0[2])
  
  base_init <- c(mu0, log(s20), z0)
  
  # two permutations only
  perms <- list(c(1,2), c(2,1))
  B <- matrix(0, nrow=6, ncol=5)
  B[1:2, 1:2] <- diag(2)        # mu
  B[3:4, 3:4] <- diag(2)        # sigma2
  B[5:6, 5]   <- c(1, -1)       # weights
  
  k <- 5
  L_vec <- rep(-Inf, 2)
  
  for(m in 1:2){
    perm <- perms[[m]]
    mu_p <- mu0[perm]
    s2_p <- s20[perm]
    w_p  <- w0[perm]; w_p <- pmax(w_p, 1e-12); w_p <- w_p / sum(w_p)
    z_p  <- log(w_p[1] / w_p[2])
    theta0_t <- c(mu_p, log(s2_p), z_p)
    
    opt <- tryCatch(
      optim(theta0_t, fn=optim_function, gr=optim_gradient,
            method="BFGS", hessian=FALSE,
            control=list(abstol=1e-8, reltol=1e-8)),
      error = function(e) NULL
    )
    if(is.null(opt)) next
    theta_t_hat <- opt$par
    
    mu_hat <- theta_t_hat[1:2]
    eta_hat <- theta_t_hat[3:4]
    z_hat <- theta_t_hat[5]
    sigma2_hat <- exp(eta_hat)
    w_hat <- softmax2(z_hat)
    theta_vec_hat <- c(mu_hat, sigma2_hat, w_hat)
    
    H_full <- -d2_log_posterior_theta_given_X(X, theta_vec_hat, bayes_acss_parameters)
    H_eff  <- t(B) %*% H_full %*% B
    
    log_det_hessian <- tryCatch({
      R <- chol(H_eff)
      2*sum(log(diag(R)))
    }, error = function(e){
      determinant(H_eff, logarithm=TRUE)$modulus
    })
    
    L_vec[m] <- log_posterior_theta_given_X(X, theta_vec_hat, bayes_acss_parameters) -
      0.5*log_det_hessian + 0.5*k*log(2*pi)
  }
  
  log_marginal_all <- logSumExp(L_vec[is.finite(L_vec)])
  return(log_marginal_all)
}

log_sampling_density <- function(X, i, theta_posterior_samples, parameters, bayes_acss_parameters, log=TRUE){
  log_numerator <- log_numerator_X_tilde(X[i], theta_posterior_samples)
  log_denominator <- log_marginal(X, bayes_acss_parameters, parameters)
  return(log_numerator - (bayes_acss_parameters$B - 1) * log_denominator)
}

get_mixture_of_normal_proposal_params <- function(
    X, i, theta_posterior_samples, parameters, bayes_acss_parameters,
    normalize_weights = TRUE, guard_negative = TRUE
){
  # Only supports J = 2
  J <- parameters$J
  if(J != 2) stop("This implementation currently supports J = 2 only.")
  
  # target log-density along coordinate i (others fixed to X)
  target_logdens_1d <- function(xi) {
    xvec <- X
    xvec[i] <- xi
    log_sampling_density(xvec, i, theta_posterior_samples, parameters, bayes_acss_parameters)
  }
  
  # piecewise quadratic fitter
  fit_two_quadratics_continuous <- function(x, y, n_grid = 401, min_side = 3, refine = TRUE){
    stopifnot(length(x) == length(y))
    keep <- is.finite(x) & is.finite(y)
    x <- as.numeric(x[keep]); y <- as.numeric(y[keep])
    n <- length(x); if(n < 5) stop("Need at least 5 finite points.")
    
    ord <- order(x); x <- x[ord]; y <- y[ord]
    rng <- range(x)
    eps <- .Machine$double.eps^(1/3) * diff(rng)
    lower <- rng[1] + eps; upper <- rng[2] - eps
    if(lower >= upper) stop("x has no spread.")
    
    make_X <- function(c){
      left  <- as.numeric(x <= c)
      right <- 1 - left
      xc <- x - c
      cbind(1,
            left  * xc,
            left  * xc^2,
            right * xc,
            right * xc^2)
    }
    sse_at <- function(c){
      nL <- sum(x <= c); nR <- n - nL
      if(nL < min_side || nR < min_side) return(Inf)
      Xd <- make_X(c)
      fit <- lm.fit(Xd, y)
      if(fit$rank < ncol(Xd)) return(Inf)
      sum(fit$residuals^2)
    }
    cs <- seq(lower, upper, length.out = n_grid)
    sse_grid <- vapply(cs, sse_at, numeric(1))
    c_start <- cs[which.min(sse_grid)]
    
    if(refine){
      step <- if(n_grid > 1) diff(cs[1:2]) else (upper - lower)/100
      a <- max(lower, c_start - 3*step)
      b <- min(upper, c_start + 3*step)
      opt <- tryCatch(optimize(sse_at, lower = a, upper = b), error = function(e) NULL)
      c_hat <- if(is.null(opt)) c_start else opt$minimum
    } else c_hat <- c_start
    
    Xhat <- make_X(c_hat)
    fit  <- lm.fit(Xhat, y)
    coef <- as.numeric(fit$coefficients)
    names(coef) <- c("alpha","bL","gL","bR","gR")
    yhat <- drop(Xhat %*% coef)
    predict_fun <- function(xnew){
      xnew <- as.numeric(xnew)
      left  <- as.numeric(xnew <= c_hat)
      right <- 1 - left
      xc <- xnew - c_hat
      alpha <- coef[1]; bL <- coef[2]; gL <- coef[3]; bR <- coef[4]; gR <- coef[5]
      alpha + left*(bL*xc + gL*xc^2) + right*(bR*xc + gR*xc^2)
    }
    list(c = c_hat, coefficients = coef, predict = predict_fun,
         x_sorted = x, y_sorted = y, yhat_sorted = yhat)
  }
  
  # build a grid, evaluate -log f, fit piecewise quadratic
  x_min <- min(X); x_max <- max(X)
  if(!is.finite(x_min) || !is.finite(x_max) || x_min == x_max){
    x_min <- X[i] - 5; x_max <- X[i] + 5
  }
  xseq <- seq(x_min, x_max, length.out = 20)
  g     <- function(xi) -target_logdens_1d(xi)       # negative log-density
  yseq <- sapply(xseq, g)
  
  fit <- fit_two_quadratics_continuous(xseq, yseq)
  co  <- fit$coefficients; c0 <- fit$c
  alpha <- co["alpha"]; bL <- co["bL"]; gL <- co["gL"]; bR <- co["bR"]; gR <- co["gR"]
  
  # extract mixture means/vars from the two quadratics
  tiny <- 1e-12
  # means (guard if curvature is non-positive)
  if(is.finite(gL) && gL > tiny) {
    muL <- c0 - bL/(2*gL)
  } else {
    # fallback: minimizer on left grid
    idxL <- which(xseq <= c0)
    muL <- if(length(idxL)) xseq[idxL][which.min(yseq[idxL])] else min(xseq)
    gL  <- max(abs(gL), tiny) # kept positive for variance
  }
  if(is.finite(gR) && gR > tiny) {
    muR <- c0 - bR/(2*gR)
  } else {
    idxR <- which(xseq > c0)
    muR <- if(length(idxR)) xseq[idxR][which.min(yseq[idxR])] else max(xseq)
    gR  <- max(abs(gR), tiny)
  }
  varL <- 1 / (2 * max(gL, tiny))
  varR <- 1 / (2 * max(gR, tiny))
  
  # heights at the minima
  QL <- function(x) alpha + bL*(x - c0) + gL*(x - c0)^2
  QR <- function(x) alpha + bR*(x - c0) + gR*(x - c0)^2
  yLmin <- QL(muL); yRmin <- QR(muR)
  
  FL <- exp(-yLmin)  # f(muL) approx
  FR <- exp(-yRmin)  # f(muR) approx
  
  # solve for weights using the 2x2 system at the two modes
  dnorm_s <- function(x, m, v) dnorm(x, mean = m, sd = sqrt(v))
  A <- matrix(c(
    dnorm_s(muL, muL, varL), dnorm_s(muL, muR, varR),
    dnorm_s(muR, muL, varL), dnorm_s(muR, muR, varR)),
    nrow = 2, byrow = TRUE
  )
  F <- c(FL, FR)
  ok <- TRUE
  w <- tryCatch(solve(A, F), error = function(e){ ok <<- FALSE; c(NA_real_, NA_real_) })
  
  # fallback or cleanup
  if(!ok || any(!is.finite(w)) || (guard_negative && any(w <= 0))){
    w <- c(FL * sqrt(2*pi*varL), FR * sqrt(2*pi*varR)) # ignore cross-terms
  }
  if(normalize_weights){
    s <- sum(w)
    if(!is.finite(s) || s <= 0) w[] <- 1/2 else w <- w / s
  }
  if(guard_negative){
    w <- pmax(w, .Machine$double.eps)
    w <- w / sum(w)
  }
  
  means   <- c(muL, muR)
  vars    <- c(varL, varR)
  weights <- c(w[1], w[2])

  list(means = means, vars = vars, weights = weights)
}

dnormmix <- function(x, lambda, mu, sigma, log = TRUE) {
  x      <- as.numeric(x)
  lambda <- as.numeric(lambda); mu <- as.numeric(mu); sigma <- as.numeric(sigma)
  w  <- lambda; w <- w / sum(w)
  log_comp <- sapply(1:length(mu), function(j) dnorm(x, mean = mu[j], sd = sigma[j], log = TRUE) + log(w[j]))
  if(log) logSumExp(log_comp) else exp(logSumExp(log_comp))
}

run_one_trial = function(seed,example = 'gaussian_mixture', print_progress = FALSE, parameters = NULL, bayes_acss_parameters = NULL){
  set.seed(seed)
  source(paste0(example,'_source.R'))
  if(is.null(parameters)) parameters = generate_parameters()
  nsignal = length(parameters$signal)
  pval_BaCSS = rep(0,nsignal)
  pval_oracle = rep(0,nsignal)
  M = bayes_acss_parameters$M
  
  for(isignal in 1:nsignal){
    experiment = generate_experiment(isignal,parameters)
    Tfun_ = function(x){Tfun(x,parameters,experiment)}
    
    if(is.null(bayes_acss_parameters)){
      bayes_acss_parameters = list(prior_params = list("mu" = list("mean" = 0, "kappa" = 1),
                                                       "sigma2" = list("alpha" = 1, "beta" = 0.5),
                                                       "w" = list("pi" = rep(2, parameters$J))),
                                   burnin = 500, thinning = 10, B = 20, 
                                   M=250, L = 1, n = parameters$example$n)
    }
    
    # generate posterior samples for beta
    prior_params = bayes_acss_parameters$prior_params
    J = parameters$J
    n = length(experiment$X)
    
    source(paste0("posterior_draw_", example, ".R"))
    theta_posterior_samples = gibbs_theta(x = experiment$X, J=J, niter = bayes_acss_parameters$burnin + bayes_acss_parameters$thinning * bayes_acss_parameters$B, prior_params = prior_params)
    
    retained_samples = seq(bayes_acss_parameters$burnin,  bayes_acss_parameters$burnin + bayes_acss_parameters$thinning * bayes_acss_parameters$B, length.out = bayes_acss_parameters$B+1)[-1]
    theta_posterior_samples = lapply(theta_posterior_samples, function(mat) mat[retained_samples, ])
    
    # generate X and Z tilde
    m_star = sample(0:bayes_acss_parameters$M, 1)
    m_star = min(m_star, bayes_acss_parameters$M - m_star)
    if(m_star == 0){
      Xcopies = rep(list(NA), M)
      X_current = experiment$X 
      
      for(iter in 1:(bayes_acss_parameters$L*bayes_acss_parameters$M)) {
        for(i in 1:n) {
          # update X_i
          mixnorm_params = get_mixture_of_normal_proposal_params(X=X_current, i=i, theta_posterior_samples = theta_posterior_samples, 
                                                                 parameters = parameters, bayes_acss_parameters = bayes_acss_parameters)
          X_i_proposed = rnormmix(1, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          acceptance_prob = exp(
            log_sampling_density(X=X_proposed, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) -
              log_sampling_density(X=X_current, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) +
              dnormmix(X_current[i], lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE) -
              dnormmix(X_i_proposed, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE))
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) cat("\rProgress: Iteration", iter, "of", M, ", Index", i, "of", n, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0)
          Xcopies[[iter / bayes_acss_parameters$L]] <- X_current
      }
    }else{
      Xcopies1 = rep(list(NA), m_star)
      Xcopies2 = rep(list(NA), bayes_acss_parameters$M - m_star)
      X_current = experiment$X 
      
      for(iter in 1:(bayes_acss_parameters$L * m_star)) {
        for(i in n:1) {
          # update X_i
          mixnorm_params = get_mixture_of_normal_proposal_params(X=X_current, i=i, theta_posterior_samples = theta_posterior_samples, 
                                                                 parameters = parameters, bayes_acss_parameters = bayes_acss_parameters)
          X_i_proposed = rnormmix(1, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          acceptance_prob = exp(
            log_sampling_density(X=X_proposed, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) -
              log_sampling_density(X=X_current, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) +
              dnormmix(X_current[i], lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE) -
              dnormmix(X_i_proposed, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE))
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) cat("\rProgress: Iteration", iter, "of", M, ", Index", i, "of", n, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0)
          Xcopies1[[iter / bayes_acss_parameters$L]] <- X_current
      }
      
      X_current = experiment$X 
      
      for(iter in 1:(bayes_acss_parameters$L * (bayes_acss_parameters$M - m_star))) {
        for(i in 1:n) {
          # update X_i
          mixnorm_params = get_mixture_of_normal_proposal_params(X=X_current, i=i, theta_posterior_samples = theta_posterior_samples, 
                                                                 parameters = parameters, bayes_acss_parameters = bayes_acss_parameters)
          X_i_proposed = rnormmix(1, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars))
          X_proposed = X_current; X_proposed[i] = X_i_proposed
          
          acceptance_prob = exp(
            log_sampling_density(X=X_proposed, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) -
              log_sampling_density(X=X_current, i=i, theta_posterior_samples, parameters, bayes_acss_parameters, log = TRUE) +
              dnormmix(X_current[i], lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE) -
              dnormmix(X_i_proposed, lambda=mixnorm_params$weights, mu = mixnorm_params$means, sigma = sqrt(mixnorm_params$vars), log = TRUE))
          X_i_retained = ifelse(runif(1) < acceptance_prob, X_proposed[i], X_current[i])
          X_current[i] = X_i_retained
          if(print_progress == TRUE) cat("\rProgress: Iteration", iter, "of", M, ", Index", i, "of", n, "acceptance prob: ", acceptance_prob, "  ", sep = " ")
        }
        if(iter %% bayes_acss_parameters$L == 0)
          Xcopies2[[iter / bayes_acss_parameters$L]] <- X_current
      }
      Xcopies = c(Xcopies1, Xcopies2)
    }
    
    pval_BaCSS[isignal] = (1+sum(unlist(lapply(Xcopies,Tfun_))>=Tfun_(experiment$X)))/(1+M)
    if(print_progress == TRUE) cat("\n\rCompleted isignal:", isignal, "of", nsignal, "\n", sep = " ")
    
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

