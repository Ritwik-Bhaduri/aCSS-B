// cpp/mixture.cpp
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

// ---------------- helpers ----------------

inline double logsumexp(const NumericVector& v) {
  int n = v.size();
  if (n == 0) return R_NegInf;
  double m = R_NegInf;
  for (int i = 0; i < n; ++i) if (R_finite(v[i]) && v[i] > m) m = v[i];
  if (!R_finite(m)) return m; // all -Inf
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += std::exp(v[i] - m);
  return m + std::log(s);
}

inline NumericVector row_logsumexp(const NumericMatrix& M) {
  int n = M.nrow(), p = M.ncol();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double m = R_NegInf;
    for (int j = 0; j < p; ++j) {
      double v = M(i, j);
      if (R_finite(v) && v > m) m = v;
    }
    if (!R_finite(m)) { out[i] = m; continue; }
    double s = 0.0;
    for (int j = 0; j < p; ++j) s += std::exp(M(i, j) - m);
    out[i] = m + std::log(s);
  }
  return out;
}

inline NumericMatrix as_matrix(SEXP obj) {
  if (Rf_isMatrix(obj)) {
    return as<NumericMatrix>(obj);
  } else {
    NumericVector v = as<NumericVector>(obj);
    NumericMatrix M(1, v.size());
    for (int j = 0; j < v.size(); ++j) M(0, j) = v[j];
    return M;
  }
}

inline bool any_nonpos(const NumericVector& v) {
  for (int i = 0; i < v.size(); ++i) if (!(v[i] > 0.0)) return true;
  return false;
}

// ---------------- likelihood (single theta) ----------------

// [[Rcpp::export]]
double likelihood(List theta, NumericVector X, bool log_ = true) {
  NumericVector mu = as<NumericVector>(theta["mu"]);
  NumericVector s2 = as<NumericVector>(theta["sigma2"]);
  NumericVector w  = as<NumericVector>(theta["w"]);
  int J = mu.size();
  int n = X.size();
  
  if (s2.size() != J || w.size() != J)
    stop("theta$mu, theta$sigma2, and theta$w must have identical lengths J.");
  
  if (any_nonpos(s2) || any_nonpos(w))
    return R_NegInf;
  
  double ll = 0.0;
  NumericVector comp(J);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < J; ++j) {
      double dx = X[i] - mu[j];
      // -0.5*((x-mu)^2/s2) - 0.5*log(s2) + log(w)
      comp[j] = -0.5 * (dx*dx) / s2[j] - 0.5 * std::log(s2[j]) + std::log(w[j]);
    }
    ll += logsumexp(comp) - 0.5 * LOG2PI;
  }
  
  return (log_ ? ll : std::exp(ll));
}

// ---------------- likelihood (multiple thetas, B x J) ----------------

// [[Rcpp::export]]
NumericVector likelihood_multiple(List theta, NumericVector X, bool log_ = true) {
  NumericMatrix mu   = as_matrix(theta["mu"]);      // B x J (or 1 x J)
  NumericMatrix s2   = as_matrix(theta["sigma2"]);  // B x J
  NumericMatrix wmat = as_matrix(theta["w"]);       // B x J
  
  int B = mu.nrow();
  int J = mu.ncol();
  int n = X.size();
  
  if (s2.nrow() != B || s2.ncol() != J || wmat.nrow() != B || wmat.ncol() != J)
    stop("theta$mu, theta$sigma2, and theta$w must have identical dimensions B x J.");
  
  NumericVector out(B);
  
  for (int b = 0; b < B; ++b) {
    // quick positivity check per row
    for (int j = 0; j < J; ++j) {
      if (!(s2(b, j) > 0.0) || !(wmat(b, j) > 0.0)) { out[b] = R_NegInf; goto next_b; }
    }
    
    {
      double ll = 0.0;
      NumericVector comp(J);
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < J; ++j) {
          double dx = X[i] - mu(b, j);
          comp[j] = -0.5 * (dx*dx) / s2(b, j) - 0.5 * std::log(s2(b, j)) + std::log(wmat(b, j));
        }
        ll += logsumexp(comp) - 0.5 * LOG2PI;
      }
      out[b] = (log_ ? ll : std::exp(ll));
    }
    
    next_b: ;
  }
  
  return out;
}

// -------- log posterior = log likelihood + log prior (conjugate priors) --------
//
// Priors (matching your gradient/Hessian):
//   μ_j | σ_j^2  ~ Normal(μ0, σ_j^2 / κ)
//   σ_j^2        ~ InvGamma(α, β)  with density proportional to s2^-(α+1) exp(-β/s2)
//   w            ~ Dirichlet(π)
// This includes the proper normalizing constants.

// [[Rcpp::export]]
double log_posterior_theta_given_X(const NumericVector& X,
                                            const NumericVector& theta_vec,
                                            const List& bayes_acss_parameters)
{
  const int len = theta_vec.size();
  if (len % 3 != 0) stop("theta_vec length must be divisible by 3.");
  const int J = len / 3;
  
  // zero-copy "views" into theta_vec (no new allocations)
  const double* mu = &theta_vec[0];
  const double* s2 = &theta_vec[J];
  const double* w  = &theta_vec[2*J];
  
  // quick validity check
  for (int j = 0; j < J; ++j) {
    if (!(s2[j] > 0.0) || !(w[j] > 0.0)) return R_NegInf;
  }
  
  // --------- precompute per-component constants ----------
  // inv_s2[j], log_s2[j], log_w[j], c[j] = -0.5*log(s2_j) + log(w_j)
  std::vector<double> inv_s2(J), log_s2(J), log_w(J), cst(J);
  for (int j = 0; j < J; ++j) {
    inv_s2[j] = 1.0 / s2[j];
    log_s2[j] = std::log(s2[j]);
    log_w[j]  = std::log(w[j]);
    cst[j]    = -0.5 * log_s2[j] + log_w[j];
  }
  
  // --------- log-likelihood (tight two-pass log-sum-exp per x_i) ----------
  double ll = 0.0;
  const int n = X.size();
  
  // Parallelize over observations if OpenMP is available.
  // Each iteration contributes one scalar to "ll", so a sum reduction is safe.
#ifdef _OPENMP
#pragma omp parallel for reduction(+:ll) schedule(static)
#endif
  for (int i = 0; i < n; ++i) {
    const double xi = X[i];
    
    // 1) find max_j v_j
    double m = R_NegInf;
    for (int j = 0; j < J; ++j) {
      const double dx = xi - mu[j];
      // v_j = -0.5 * dx^2 * inv_s2[j] + cst[j]
      const double vj = -0.5 * (dx * dx) * inv_s2[j] + cst[j];
      if (vj > m) m = vj;
    }
    
    // 2) sum exp(v_j - m)
    double sum_exp = 0.0;
    for (int j = 0; j < J; ++j) {
      const double dx = xi - mu[j];
      const double vj = -0.5 * (dx * dx) * inv_s2[j] + cst[j];
      sum_exp += std::exp(vj - m);
    }
    
    ll += (m + std::log(sum_exp)) - 0.5 * LOG2PI;
  }
  
  // --------- priors ----------
  // Structure:
  //   mu | s2  ~ Normal(mu0, s2 / kappa)
  //   s2      ~ InvGamma(alpha, beta)
  //   w       ~ Dirichlet(pi)
  const List  pp   = bayes_acss_parameters["prior_params"];
  const List  p_mu = pp["mu"];
  const List  p_s2 = pp["sigma2"];
  const List  p_w  = pp["w"];
  
  const double mu0   = as<double>(p_mu["mean"]);
  const double kappa = as<double>(p_mu["kappa"]);
  const double alpha = as<double>(p_s2["alpha"]);
  const double beta  = as<double>(p_s2["beta"]);
  NumericVector pi_vec = p_w["pi"];
  if (pi_vec.size() != J) stop("Length of prior_params$w$pi must equal J.");
  
  // log p(mu | s2): sum_j [ -0.5*log(2π*s2/kappa) - 0.5*kappa*(mu-mu0)^2/s2 ]
  double lprior = 0.0;
  for (int j = 0; j < J; ++j) {
    const double dmu = mu[j] - mu0;
    lprior += -0.5 * (LOG2PI + log_s2[j] - std::log(kappa))
      -0.5 * kappa * dmu * dmu * inv_s2[j];
  }
  
  // log p(s2): J*(alpha*log(beta) - lgamma(alpha)) + sum_j [ -(alpha+1)log s2 - beta/s2 ]
  lprior += J * (alpha * std::log(beta) - lgamma(alpha));
  for (int j = 0; j < J; ++j) {
    lprior += -(alpha + 1.0) * log_s2[j] - beta * inv_s2[j];
  }
  
  // log p(w): lgamma(sum pi) - sum lgamma(pi_j) + sum (pi_j - 1) * log w_j
  double pi_sum = 0.0;
  for (int j = 0; j < J; ++j) pi_sum += pi_vec[j];
  lprior += lgamma(pi_sum);
  for (int j = 0; j < J; ++j) lprior -= lgamma(pi_vec[j]);
  for (int j = 0; j < J; ++j) lprior += (pi_vec[j] - 1.0) * log_w[j];
  
  return ll + lprior;
}

// ---------------- gradient of log posterior ----------------

// [[Rcpp::export]]
NumericVector d_log_posterior_theta_given_X(NumericVector X,
                                                NumericVector theta_vec,
                                                List bayes_acss_parameters) {
  int len = theta_vec.size();
  if (len % 3 != 0) stop("theta_vec length must be divisible by 3.");
  int J = len / 3;
  int n = X.size();
  
  NumericVector mu  = theta_vec[ seq(0,          J - 1) ];
  NumericVector s2  = theta_vec[ seq(J,      2 * J - 1) ];
  NumericVector w   = theta_vec[ seq(2 * J,  3 * J - 1) ];
  
  // priors
  List pp    = bayes_acss_parameters["prior_params"];
  List p_mu  = pp["mu"];
  List p_s2  = pp["sigma2"];
  List p_w   = pp["w"];
  
  double mu0   = as<double>(p_mu["mean"]);
  double kappa = as<double>(p_mu["kappa"]);
  double alpha = as<double>(p_s2["alpha"]);
  double beta  = as<double>(p_s2["beta"]);
  NumericVector pi_vec = p_w["pi"];
  if (pi_vec.size() != J) stop("Length of prior_params$w$pi must equal J.");
  
  // responsibilities r_{i,j}
  NumericMatrix r(n, J);
  // we also need dx, a' (for each (i,j))
  NumericMatrix dx(n, J);
  NumericMatrix a_prime(n, J);
  
  NumericVector log_num(J);
  
  for (int i = 0; i < n; ++i) {
    // build log_num for all j, then normalize
    for (int j = 0; j < J; ++j) {
      double d = X[i] - mu[j];
      dx(i, j) = d;
      // log φ = -0.5*log(2π s2) - (dx^2)/(2 s2)
      double log_phi = -0.5 * (LOG2PI + std::log(s2[j])) - 0.5 * (d * d) / s2[j];
      log_num[j] = std::log(w[j]) + log_phi;
      
      // a' = ∂ log φ / ∂ s2 = -0.5/s2 + 0.5*(dx^2)/s2^2
      a_prime(i, j) = -0.5 / s2[j] + 0.5 * (d * d) / (s2[j] * s2[j]);
    }
    double lse = logsumexp(log_num);
    for (int j = 0; j < J; ++j) {
      r(i, j) = std::exp(log_num[j] - lse);
    }
  }
  
  // gradient parts
  NumericVector dmu(J), ds2(J), dw(J);
  for (int j = 0; j < J; ++j) {
    double s2j = s2[j];
    double sum_mu = 0.0, sum_s2 = 0.0, sum_w = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_mu += r(i, j) * (dx(i, j) / s2j);
      sum_s2 += r(i, j) * a_prime(i, j);
      sum_w  += r(i, j);
    }
    dmu[j] = sum_mu - kappa * (mu[j] - mu0) / s2j;
    
    ds2[j] = sum_s2
    - 1.0 / (2.0 * s2j)
      + kappa * (mu[j] - mu0) * (mu[j] - mu0) / (2.0 * s2j * s2j)
      - (alpha + 1.0) / s2j + beta / (s2j * s2j);
      
      dw[j]  = sum_w / w[j] + (pi_vec[j] - 1.0) / w[j];
  }
  
  NumericVector out(3 * J);
  for (int j = 0; j < J; ++j) {
    out[j] = dmu[j];
    out[J + j] = ds2[j];
    out[2 * J + j] = dw[j];
  }
  return out;
}

// ---------------- Hessian (second derivatives) of log posterior ----------------

// [[Rcpp::export]]
NumericMatrix d2_log_posterior_theta_given_X(NumericVector X,
                                                 NumericVector theta_vec,
                                                 List bayes_acss_parameters) {
  int len = theta_vec.size();
  if (len % 3 != 0) stop("theta_vec length must be divisible by 3.");
  int J = len / 3;
  int n = X.size();
  
  NumericVector mu  = theta_vec[ seq(0,          J - 1) ];
  NumericVector s2  = theta_vec[ seq(J,      2 * J - 1) ];
  NumericVector w   = theta_vec[ seq(2 * J,  3 * J - 1) ];
  
  // priors
  List pp    = bayes_acss_parameters["prior_params"];
  List p_mu  = pp["mu"];
  List p_s2  = pp["sigma2"];
  List p_w   = pp["w"];
  
  double mu0   = as<double>(p_mu["mean"]);
  double kappa = as<double>(p_mu["kappa"]);
  double alpha = as<double>(p_s2["alpha"]);
  double beta  = as<double>(p_s2["beta"]);
  NumericVector pi_vec = p_w["pi"];
  if (pi_vec.size() != J) stop("Length of prior_params$w$pi must equal J.");
  
  // compute r, 1-r, dx, a', a''
  NumericMatrix r(n, J), omr(n, J);
  NumericMatrix dx(n, J), a_prime(n, J), a_dblprime(n, J);
  NumericVector log_num(J);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < J; ++j) {
      double d = X[i] - mu[j];
      dx(i, j) = d;
      double log_phi = -0.5 * (LOG2PI + std::log(s2[j])) - 0.5 * (d * d) / s2[j];
      log_num[j] = std::log(w[j]) + log_phi;
      
      a_prime(i, j)    = -0.5 / s2[j] + 0.5 * (d * d) / (s2[j] * s2[j]);
      a_dblprime(i, j) =  0.5 / (s2[j] * s2[j]) - (d * d) / (s2[j] * s2[j] * s2[j]);
    }
    double lse = logsumexp(log_num);
    for (int j = 0; j < J; ++j) {
      r(i, j) = std::exp(log_num[j] - lse);
      omr(i, j) = 1.0 - r(i, j);
    }
  }
  
  int D = 3 * J;
  NumericMatrix H(D, D);
  IntegerVector idx_mu  = seq(0,     J - 1);
  IntegerVector idx_s2  = seq(J,     2 * J - 1);
  IntegerVector idx_w   = seq(2 * J, 3 * J - 1);
  
  // μ–μ block
  for (int j = 0; j < J; ++j) {
    double s2j = s2[j];
    double diagv = 0.0;
    for (int i = 0; i < n; ++i) {
      diagv += r(i, j) * omr(i, j) * (dx(i, j) * dx(i, j)) / (s2j * s2j) - r(i, j) / s2j;
    }
    diagv -= kappa / s2j;
    H(idx_mu[j], idx_mu[j]) = diagv;
    
    for (int k = 0; k < J; ++k) if (k != j) {
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
        val += - r(i, j) * r(i, k) * (dx(i, j) / s2j) * (dx(i, k) / s2[k]);
      }
      H(idx_mu[j], idx_mu[k]) = val;
      H(idx_mu[k], idx_mu[j]) = val;
    }
  }
  
  // σ2–σ2 block
  for (int j = 0; j < J; ++j) {
    double s2j = s2[j];
    double diagv = 0.0;
    for (int i = 0; i < n; ++i) {
      diagv += r(i, j) * omr(i, j) * (a_prime(i, j) * a_prime(i, j)) + r(i, j) * a_dblprime(i, j);
    }
    diagv +=  1.0 / (2.0 * s2j * s2j)
      - kappa * (mu[j] - mu0) * (mu[j] - mu0) / (s2j * s2j * s2j)
      + (alpha + 1.0) / (s2j * s2j)
      - 2.0 * beta / (s2j * s2j * s2j);
      H(idx_s2[j], idx_s2[j]) = diagv;
      
      for (int k = 0; k < J; ++k) if (k != j) {
        double val = 0.0;
        for (int i = 0; i < n; ++i) {
          val += - r(i, j) * r(i, k) * a_prime(i, j) * a_prime(i, k);
        }
        H(idx_s2[j], idx_s2[k]) = val;
        H(idx_s2[k], idx_s2[j]) = val;
      }
  }
  
  // w–w block
  for (int j = 0; j < J; ++j) {
    double diagv = 0.0;
    for (int i = 0; i < n; ++i) {
      double term = r(i, j) / w[j];
      diagv += - term * term;
    }
    diagv += - (pi_vec[j] - 1.0) / (w[j] * w[j]);
    H(idx_w[j], idx_w[j]) = diagv;
    
    for (int k = 0; k < J; ++k) if (k != j) {
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
        val += - (r(i, j) / w[j]) * (r(i, k) / w[k]);
      }
      H(idx_w[j], idx_w[k]) = val;
      H(idx_w[k], idx_w[j]) = val;
    }
  }
  
  // μ–σ2 cross
  for (int j = 0; j < J; ++j) {
    double s2j = s2[j];
    double val_jj = 0.0;
    for (int i = 0; i < n; ++i) {
      val_jj += r(i, j) * omr(i, j) * a_prime(i, j) * (dx(i, j) / s2j)
      - r(i, j) * (dx(i, j) / (s2j * s2j));
    }
    val_jj += kappa * (mu[j] - mu0) / (s2j * s2j);
    H(idx_mu[j], idx_s2[j]) = val_jj;
    H(idx_s2[j], idx_mu[j]) = val_jj;
    
    for (int k = 0; k < J; ++k) if (k != j) {
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
        val += - r(i, j) * r(i, k) * a_prime(i, k) * (dx(i, j) / s2j);
      }
      H(idx_mu[j], idx_s2[k]) = val;
      H(idx_s2[k], idx_mu[j]) = val;
    }
  }
  
  // μ–w cross
  for (int j = 0; j < J; ++j) {
    double s2j = s2[j];
    for (int k = 0; k < J; ++k) {
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
        double indicator = (j == k) ? 1.0 : 0.0;
        val += (1.0 / w[k]) * ((indicator * r(i, k) - r(i, j) * r(i, k)) * (dx(i, j) / s2j));
      }
      H(idx_mu[j], idx_w[k]) = val;
      H(idx_w[k], idx_mu[j]) = val;
    }
  }
  
  // σ2–w cross
  for (int j = 0; j < J; ++j) {
    for (int k = 0; k < J; ++k) {
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
        double indicator = (j == k) ? 1.0 : 0.0;
        val += (1.0 / w[k]) * ( r(i, k) * (indicator - r(i, j)) * a_prime(i, j) );
      }
      H(idx_s2[j], idx_w[k]) = val;
      H(idx_w[k], idx_s2[j]) = val;
    }
  }
  
  return H;
}

// ---------------- log_numerator_X_tilde ----------------
//
// Mirrors your R function. Accepts X_i (scalar) and theta_posterior_samples
// where $mu, $sigma2, $w can be either vectors (J) or matrices (B x J).
// Returns a single scalar: sum over rows of rowLogSumExp(...) - 0.5*log(2π).

// [[Rcpp::export]]
double log_numerator_X_tilde(double X_i, List theta_posterior_samples) {
  NumericMatrix mu = as_matrix(theta_posterior_samples["mu"]);
  NumericMatrix s2 = as_matrix(theta_posterior_samples["sigma2"]);
  NumericMatrix w  = as_matrix(theta_posterior_samples["w"]);
  
  int B = mu.nrow();
  int J = mu.ncol();
  
  if (s2.nrow() != B || s2.ncol() != J || w.nrow() != B || w.ncol() != J)
    stop("theta$mu, theta$sigma2, and theta$w must have identical dimensions.");
  
  double out = 0.0;
  NumericVector comps(J);
  for (int b = 0; b < B; ++b) {
    for (int j = 0; j < J; ++j) {
      if (!(s2(b, j) > 0.0) || !(w(b, j) > 0.0)) return R_NegInf;
      double dx = X_i - mu(b, j);
      comps[j] = -0.5 * (dx*dx) / s2(b, j) - 0.5 * std::log(s2(b, j)) + std::log(w(b, j));
    }
    out += logsumexp(comps) - 0.5 * LOG2PI;
  }
  return out;
}
