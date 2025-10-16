// log_marginal_cpp.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double log_posterior_t_given_X_cpp(const arma::vec& X, const arma::vec& t_vec, 
                                   const List& bayes_acss_parameters, const List& parameters, const List& experiment) {
  // Extract parameters
  double mu_1 = bayes_acss_parameters["mu_1"];
  double tau_1 = bayes_acss_parameters["tau_1"];
  double mu_2 = bayes_acss_parameters["mu_2"];
  double tau_2 = bayes_acss_parameters["tau_2"];
  double sigma = parameters["sigma"];
  int k = parameters["k"];
  int n = parameters["n"];
  
  // Retrieve the nested list and access Z
  List example = experiment["example"];  // Get the "example" list within "experiment"
  arma::vec Z = as<arma::vec>(example["Z"]);  // Now access "Z" within "example"
  
  int p = t_vec.n_elem;
  
  // Compute h_t
  arma::mat h_t(Z.n_elem, k + 1);
  h_t.col(0).ones();
  h_t.col(1) = Z;
  for (int i = 0; i < Z.n_elem; i++) {
    for (int j = 0; j < p; j++) {
      h_t(i, j + 2) = std::max(Z[i] - t_vec[j], 0.0);
    }
  }
  
  // Calculate X_t_mean
  arma::vec X_t_mean = mu_1 * sum(h_t, 1);
  
  // Calculate A
  arma::mat A = (1.0 / (tau_1 * tau_1)) * arma::eye(k + 1, k + 1) + (1.0 / (sigma * sigma)) * h_t.t() * h_t;
  
  // Compute A_inv using Cholesky decomposition
  arma::mat A_inv = inv_sympd(A);
  
  // Compute Sigma_inv
  arma::mat Sigma_inv = (1.0 / (sigma * sigma)) * arma::eye(n, n) - (1.0 / (sigma * sigma * sigma * sigma)) * h_t * A_inv * h_t.t();
  
  // Calculate log determinant of A scaled by tau_1^2
  double log_det_A_tau = log(det(A * tau_1 * tau_1)) * -0.5;
  
  // Calculate the final expression
  arma::vec diff = X - X_t_mean;
  double result = log_det_A_tau - 0.5 * as_scalar(diff.t() * Sigma_inv * diff) - sum(square(t_vec - mu_2) / (2.0 * tau_2 * tau_2));
  
  return result;
}


// [[Rcpp::export]]
double log_marginal_cpp(const arma::vec& X, const List& bayes_acss_parameters, 
                        const List& parameters, const List& experiment, const List& precomputations,
                        double min_val = -10.0, double max_val = 10.) {
  
  // Extract Z_sorted from precomputations
  arma::vec Z_sorted = as<arma::vec>(precomputations["Z_sorted"]);
  int n = as<int>(parameters["n"]);
  double integration = 0.0;
  
  // Build extended grid including tails
  arma::vec Z_ext(n + 2);
  Z_ext[0] = min_val;
  Z_ext.subvec(1, n) = Z_sorted;
  Z_ext[n + 1] = max_val;
  
  // Loop over intervals to calculate the integral
  for (int i = 0; i < Z_ext.n_elem - 1; i++) {
    // Generate points for integration within each interval
    arma::vec points = linspace<vec>(Z_ext[i], Z_ext[i + 1], 21);
    arma::vec fun_vals(points.n_elem);
    
    // Compute the values of log_posterior_t_given_X_cpp at each point
    for (int j = 0; j < points.n_elem; j++) {
      arma::vec t_vec = { points[j] };  // Convert point to a one-element vector
      fun_vals[j] = exp(log_posterior_t_given_X_cpp(X, t_vec, bayes_acss_parameters, parameters, experiment));
    }
    
    // Perform the trapezoidal integration for this segment
    double trapz_integral = 0.0;
    for (int j = 0; j < fun_vals.n_elem - 1; j++) {
      double dx = points[j + 1] - points[j];
      trapz_integral += (fun_vals[j] + fun_vals[j + 1]) / 2.0 * dx;
    }
    
    // Accumulate the result
    integration += trapz_integral;
  }
  
  // Return the logarithm of the integration result
  return log(integration);
}
