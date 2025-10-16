# Function to normalize probabilities
normalize <- function(x) {
  return(x / sum(x))
}

# Sample latent variables Z
sample_z <- function(x, w, mu, sigma2) {
  J <- length(w)
  n <- length(x)
  log_prob <- matrix(0, nrow = J, ncol = n)
  for (i in 1:J) {
    log_prob[i, ] <- log(w[i]) + dnorm(x, mean = mu[i], sd = sqrt(sigma2[i]), log = TRUE)
  }
  prob <- exp(log_prob - apply(log_prob, 2, max)) # prevent underflow
  prob <- apply(prob, 2, normalize)
  z <- apply(prob, 2, function(p) sample(1:J, size = 1, prob = p, replace = TRUE))
  return(z)
}

# Sample mixture weights w using Dirichlet distribution
sample_w <- function(z, J, w_prior) {
  counts <- table(factor(z, levels = 1:J))
  w <- gtools::rdirichlet(1, counts + w_prior)
  return(w)
}

# Sample means mu
sample_mu <- function(x, z, J, sigma2, mu_prior) {
  mu <- rep(0, J)
  for (j in 1:J) {
    indices <- which(z == j)
    n_j <- length(indices)
    if (n_j > 0) {
      sample_mean <- mean(x[indices])
      # Update the posterior variance and mean calculations to reflect the new prior
      post_var <- 1 / (n_j / sigma2[j] + mu_prior$kappa / sigma2[j])
      post_mean <- (sample_mean * n_j + mu_prior$mean * mu_prior$kappa) / (n_j + mu_prior$kappa)
      mu[j] <- rnorm(1, post_mean, sqrt(post_var))
    } else {
      # When no data points are assigned to this cluster, sample from the prior
      mu[j] <- rnorm(1, mu_prior$mean, sqrt(sigma2[j] / mu_prior$kappa))
    }
  }
  return(mu)
}

# Sample variances sigma2
sample_sigma2 <- function(x, z, mu, J, sigma2_prior, mu_prior) {
  sigma2 <- rep(0, J)
  for (i in 1:J) {
    indices <- which(z == i)
    n_i <- length(indices)
    if (n_i > 0) {
      alpha_post <- sigma2_prior$alpha + n_i / 2 + 1/2
      beta_post <- sigma2_prior$beta + 0.5 * sum((x[indices] - mu[i])^2) + 0.5 * mu_prior$kappa * (mu[i] - mu_prior$mean)^2
      sigma2[i] <- 1 / rgamma(1, alpha_post, rate = beta_post)
    } else {
      sigma2[i] <- 1 / rgamma(1, sigma2_prior$alpha + 1/2, rate = sigma2_prior$beta + 0.5 * mu_prior$kappa * (mu[i] - mu_prior$mean)^2) # No data points in cluster, sample from prior
    }
  }
  return(sigma2)
}

# Gibbs sampler function
gibbs_theta <- function(x, J, niter = 1000, prior_params) {
  mu_prior = prior_params$mu
  sigma2_prior = prior_params$sigma2
  w_prior = prior_params$w$pi
  
  n <- length(x)
  w <- rep(1/J, J)

  q = quantile(x, seq(0,1,length.out = J+1)); q[1]=-Inf; q[length(q)]=Inf
  z = sapply(x, function(y) sum(y-q>0))
  mu = sapply(1:J, function(j) mean(x[z==j]))
  sigma2 = sapply(1:J, function(j) var(x[z==j]))
  
  # Containers for storing samples
  res <- list(mu = matrix(0, nrow = niter, ncol = J), 
              sigma2 = matrix(0, nrow = niter, ncol = J),
              w = matrix(0, nrow = niter, ncol = J), 
              z = matrix(0, nrow = niter, ncol = n))
  
  for (iter in 1:niter) {
    z <- sample_z(x, w, mu, sigma2)
    w <- sample_w(z, J, w_prior)
    mu <- sample_mu(x, z, J, sigma2, mu_prior)
    sigma2 <- sample_sigma2(x, z, mu, J, sigma2_prior, mu_prior)
    
    res$mu[iter, ] <- mu
    res$sigma2[iter, ] <- sigma2
    res$w[iter, ] <- w
    res$z[iter, ] <- z
  }
  
  return(res)
}
