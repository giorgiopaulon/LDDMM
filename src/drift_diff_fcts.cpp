
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <R_ext/Utils.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <rgen.h> 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(rgen)]]
// [[Rcpp::plugins(cpp11)]]


double dinvgaussian_cpp(arma::vec tau, arma::vec mu, arma::vec lambda, 
                        arma::vec delta, bool log){
  
  double out; 
  
  out = arma::accu(0.5 * (arma::log(lambda) - std::log(2.0 * arma::datum::pi)) - 
    3.0/2.0 * arma::log(tau - delta) - 
    lambda % arma::square(tau - delta - mu) / (2.0 * arma::square(mu) % (tau - delta)));
  
  return(out); 
}


double pinvgaussian_cpp(arma::vec tau, arma::vec mu, arma::vec lambda, 
                        arma::vec delta, bool log){
  
  double out; 
  out = arma::accu(arma::log(arma::abs(
    1.0 - 
      arma::normcdf(  arma::sqrt(lambda / (tau - delta)) % ((tau - delta) / mu - 1.0) ) - 
      arma::exp( 2.0 * lambda / mu +   
      arma::log(arma::normcdf( - arma::sqrt(lambda / (tau - delta)) % ((tau - delta) / mu + 1.0) )))
  )));
  
  return(out); 
}


//' Log-likelihood computation
//' 
//' Compute the log-likelihood for the drift-diffusion model, including the
//' censored data contribution.
//'
//' @param tau vector of size n containing the response times
//' @param mu matrix of size (n x d1) containing the drift parameters 
//' corresponding to the n response times for each possible d1 decision
//' @param b matrix of size (n x d1) containing the boundary parameters 
//' corresponding to the n response times for each possible d1 decision
//' @param delta vector of size n containing the offset parameters corresponding
//' to the n response times
//' @param cens vector of size n containing censoring indicators (1 censored, 0
//' not censored) corresponding to the n response times
//' @param D (n x 2) matrix whose first column has the n input stimuli, and whose second column has the n decision categories
//' @param log should the results be returned on the log scale?
// [[Rcpp::export]]
double log_likelihood(arma::vec tau, arma::mat mu, arma::mat b, 
                      arma::vec delta, arma::uvec cens, arma::umat D, 
                      bool log){
  
  unsigned int d_1 = mu.n_cols;
  arma::uvec one_to_d = arma::regspace<arma::uvec>(0, d_1 - 1);
  arma::uvec idx_cens = arma::find(cens == 1);
  arma::uvec idx_noncens = arma::find(cens == 0);
  arma::uvec idx_nsupp = arma::find(tau <= delta);
  arma::uvec yprime;
  arma::uvec idx_succ;
  
  arma::vec out(d_1*(d_1 + 1), arma::fill::zeros);
  out -= arma::datum::inf;
  unsigned int pos = 0;
  
  // Computation is done on the log scale to prevent underflow/overflow
  if (idx_nsupp.n_elem == 0){
    
    for (unsigned int y_temp = 0; y_temp < d_1; y_temp++){
      idx_succ = arma::find( (D.col(1) == (y_temp + 1)) && (cens == 0) );
      
      if (idx_succ.n_elem > 0){
        out(pos) = dinvgaussian_cpp(tau.elem(idx_succ), b.submat(idx_succ, arma::regspace<arma::uvec>(y_temp, y_temp))/mu.submat(idx_succ, arma::regspace<arma::uvec>(y_temp, y_temp)), 
            arma::square(b.submat(idx_succ, arma::regspace<arma::uvec>(y_temp, y_temp))), delta.elem(idx_succ), true);
      }
      else{
        out(pos) = 0.0;
      }
      pos += 1;
      yprime = one_to_d.elem(arma::find(one_to_d != y_temp));
      
      for (unsigned int j = 0; j < d_1 - 1; j++){
        if (idx_succ.n_elem > 0){
          out(pos) = pinvgaussian_cpp(tau.elem(idx_succ), b.submat(idx_succ, arma::regspace<arma::uvec>(yprime(j), yprime(j)))/mu.submat(idx_succ, arma::regspace<arma::uvec>(yprime(j), yprime(j))), 
              arma::square(b.submat(idx_succ, arma::regspace<arma::uvec>(yprime(j), yprime(j)))), delta.elem(idx_succ), true);
        }
        else{
          out(pos) = 0.0;
        }
        pos += 1;
      }
    }
    if (idx_cens.n_elem > 0){
      for (unsigned int y_temp = 0; y_temp < d_1; y_temp++){
        out(pos) = pinvgaussian_cpp(tau.elem(idx_cens), b.submat(idx_cens, arma::regspace<arma::uvec>(y_temp, y_temp))/mu.submat(idx_cens, arma::regspace<arma::uvec>(y_temp, y_temp)),
            arma::square(b.submat(idx_cens, arma::regspace<arma::uvec>(y_temp, y_temp))), delta.elem(idx_cens), true);;
        
        pos += 1;
      }
    }
    else{
      for (unsigned int y_temp = 0; y_temp < d_1; y_temp++){
        out(pos) = 0.0;
        pos += 1;
      }
    }
  }
  
  if (out.has_nan()){
    out.zeros();
    out -= arma::datum::inf;
  }
  
  if (log){
    return (arma::accu(out));
  } 
  else{
    return (exp(arma::accu(out)));
  }
} 


// Count Frequencies
// 
// Analogous of the table function in R: computes the counts of the elements 
// of an integer valued vector.
//
// @param x input vector of length n
// @param K maximum value that x can have (the support of x is /{1, ..., K/})
// @return vector of counts for every category
// [[Rcpp::export]]
arma::uvec table_int(arma::uvec x, unsigned int K){
  
  unsigned int n = x.n_elem;
  arma::uvec count(K, arma::fill::zeros); // here we consider that x has values in {1, ..., max(x)}
  
  for (unsigned int i = 0; i < n; i++){
    count[x[i] - 1] += 1;
  }
  return (count);
}


// Stable Sum of the Rows of a Matrix
// 
// Sums the rows of a matrix defined on the log scale, returns the vector on 
// the log scale.
//
// @param Q input matrix (with unnormalized log probabilities on the rows)
// @return vector with the sums (on the log scale)
// [[Rcpp::export]]
arma::vec sum_rows_log(arma::mat Q){
  
  unsigned int S = Q.n_cols;
  arma::vec out(S, arma::fill::zeros);
  double m_star; 
  
  for (unsigned hh = 0; hh < S; hh++){
    m_star = max(Q.row(hh));
    if (m_star != - arma::datum::inf){
      out(hh) = m_star + log(arma::accu(exp(Q.row(hh) - m_star)));
    }
    else{
      out(hh) = - arma::datum::inf;
    }
  }
  
  return (out);
}


// Stable Normalization
// 
// Normalize a vector on the log scale with the log-sum-exp trick.
// 
// @param x vector (on the log-scale) to be normalized
// @return normalized vector on the exponentiated scale
// [[Rcpp::export]]
arma::vec normalise_log(arma::vec x){
  unsigned int n = x.n_elem;
  
  double c = x.max();
  x -= c;
  
  arma::vec out(n);
  out = arma::exp(x)/arma::accu(arma::exp(x));
  
  return (out);
}


// Dirichlet Distribution
// 
// Samples from the Dirichlet distribution.
//
// @param deltas vector of scale parameters
// @return one sample from the Dirichlet distribution
// [[Rcpp::export]]
arma::vec rdirichlet(arma::vec deltas){
  arma::vec out = rgen::rdirichlet(deltas);
  
  return (out);
}



// Count the Clustering Assignments
// 
// Given the cluster memberships over time, count the clustering assignments 
// to update the transition matrix.
//
// @param z matrix with cluster membership labels over time 
// @param M max number of clusters 
// @return count the transitions between states and return it into a (M x M) 
// matrix
// [[Rcpp::export]]
arma::umat count_assign(arma::umat z, unsigned int M){
  
  unsigned int d_1 = z.n_rows;
  unsigned int K = z.n_cols;
  
  arma::umat count(M, M, arma::fill::zeros); 
  for (unsigned int i = 0; i < d_1; i++){
    for (unsigned int j = 1; j < K; j++){
      count(z(i,j-1) - 1,z(i,j) - 1) += 1;
    }
  }
  return (count);
}


// Random effects sampling
//
// Sample the random effect parameters for the drift parameters using 
// Metropolis Hastings
//
// @param tau vector of size n containing the response times
// @param D (n x 2) matrix whose first column has the n input stimuli, and whose second column has the n decision categories
// @param cens vector of size n containing censoring indicators (1 censored, 0
// not censored) corresponding to the n response times
// @param beta_u_old old random effects spline coefficients
// @param delta_dat vector of size n containing the offset parameters corresponding
// to the n response times
// @param b_dat matrix of size (n x d1) containing the boundary parameters 
// corresponding to the n response times for each possible d1 decision
// @param B_beta_dat old random effects spline coefficients expanded on the B Spline bases
// @param mu_dat_old matrix of size (n x d1) containing the old drift parameters 
// corresponding to the n response times for each possible d1 decision
// @param B expansion of the blocks on the spline basis
// @param P smoothness inducing matrix (first order differences)
// @param ind vector denoting the individuals
// @param time vector denoting the time index
// @param sigma2_us smoothness parameter for the random effects
// @param sigma2_ua variance parameter for the random effects
// @param sd_beta_u standard deviation for the random walk proposal
// @param acc_beta_u acceptance rate for the random effects spline coefficients
// @return new random effects spline coefficients
// [[Rcpp::export]]
Rcpp::List sample_reff_mu(arma::vec tau, arma::umat D, arma::uvec cens, 
                          arma::mat beta_u_old, 
                          arma::vec delta_dat, arma::mat b_dat, 
                          arma::mat B_beta_dat, arma::mat mu_dat_old, 
                          arma::mat B, arma::mat P, arma::uvec ind, 
                          arma::uvec time, double sigma2_us, 
                          double sigma2_ua, arma::mat sd_beta_u, 
                          arma::mat acc_beta_u){
  
  arma::uvec un_ind = arma::unique(ind);
  unsigned int n_ind = un_ind.n_elem;
  unsigned int J = B.n_cols;
  unsigned int n = tau.n_elem;
  unsigned int d_1 = b_dat.n_cols;
  
  arma::mat beta_u_prop(n_ind, J, arma::fill::zeros);
  arma::mat B_beta_u_dat(n, d_1, arma::fill::zeros);
  arma::mat B_beta_u_prop_dat, mu_dat_prop, temp;
  arma::uvec idx_it, idx_i;
  double l_u, alpha_acc, logpost_prop, logpost_old, logpr_prop, logpr_old;
  for (unsigned int i = 1; i <= n_ind; i++){
    for (unsigned int k = 0; k < J; k++){
      // Propose new value
      beta_u_prop.row(i-1) = beta_u_old.row(i-1);
      beta_u_prop(i - 1,k) = beta_u_old(i - 1,k) + sd_beta_u(i - 1,k) * arma::randn();
      
      if (k == 0){ // only data at t = 1 influence the first coefficient
        idx_it = arma::find( (ind == i) && (time == 1) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k + 1) - beta_u_prop(i - 1,k), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k + 1) - beta_u_old(i - 1,k), 2.0) / sigma2_us);
      }
      else if (k == (J - 1)){ // only data at t = T influence the last coefficient
        idx_it = arma::find( (ind == i) && (time == (J - 1)) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k) - beta_u_prop(i - 1,k - 1), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k) - beta_u_old(i - 1,k - 1), 2.0) / sigma2_us);
      }
      else { // data at t = {k-1,k} influence the kth coefficient
        idx_it = arma::find( (ind == i) && ((time == k) || (time == k+1)) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k + 1) - beta_u_prop(i - 1,k), 2.0) / sigma2_us + 
          std::pow(beta_u_prop(i - 1,k) - beta_u_prop(i - 1,k - 1), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k + 1) - beta_u_old(i - 1,k), 2.0) / sigma2_us + 
          std::pow(beta_u_old(i - 1,k) - beta_u_old(i - 1,k - 1), 2.0) / sigma2_us);
      }
      
      // Compute the implied \mu
      B_beta_u_prop_dat = B.rows(idx_it) * beta_u_prop.row(i-1).t();
      temp = B_beta_dat.rows(idx_it);
      temp.each_col() += B_beta_u_prop_dat;
      mu_dat_prop = arma::exp(temp);
      
      // Evaluate log posterior
      logpost_prop = log_likelihood(tau.elem(idx_it), 
                                    mu_dat_prop, 
                                    b_dat.rows(idx_it), 
                                    delta_dat.elem(idx_it),
                                    cens.elem(idx_it),
                                    D.rows(idx_it), true) + 
                                      logpr_prop;
      logpost_old = log_likelihood(tau.elem(idx_it), 
                                   mu_dat_old.rows(idx_it), 
                                   b_dat.rows(idx_it), 
                                   delta_dat.elem(idx_it), 
                                   cens.elem(idx_it),
                                   D.rows(idx_it), true) + 
                                     logpr_old;
      
      // Acceptance rate
      alpha_acc = logpost_prop - logpost_old;
      l_u = std::log(arma::randu());
      if (l_u < alpha_acc){
        beta_u_old(i - 1,k) = beta_u_prop(i - 1,k);
        mu_dat_old.rows(idx_it) = mu_dat_prop;
        acc_beta_u(i - 1,k) += 1;
      }
    }
    idx_i = arma::find(ind == i);
    for (unsigned j = 0; j < d_1; j++){
      B_beta_u_dat.submat(idx_i, arma::regspace<arma::uvec>(j, j)) = B.rows(idx_i) * beta_u_old.row(i-1).t();
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("B_beta_u_dat")=B_beta_u_dat,
    Rcpp::Named("beta_u_old")=beta_u_old, 
    Rcpp::Named("acc_beta_u")=acc_beta_u
  );
}


// Random effects sampling
// 
// Sample the random effect parameters for the boundary parameters using 
// Metropolis Hastings
//
// @param tau vector of size n containing the response times
// @param D (n x 2) matrix whose first column has the n input stimuli, and whose second column has the n decision categories
// @param cens vector of size n containing censoring indicators (1 censored, 0
// not censored) corresponding to the n response times
// @param beta_u_old old random effects spline coefficients
// @param delta_dat vector of size n containing the offset parameters corresponding
// to the n response times
// @param B_beta_dat old random effects spline coefficients expanded on the B Spline bases
// @param b_dat_old matrix of size (n x d1) containing the old boundary parameters 
// corresponding to the n response times for each possible d1 decision
// @param mu_dat matrix of size (n x d1) containing the drift parameters 
// corresponding to the n response times for each possible d1 decision
// @param B expansion of the blocks on the spline basis
// @param P smoothness inducing matrix (first order differences)
// @param ind vector denoting the individuals
// @param time vector denoting the time index
// @param sigma2_us smoothness parameter for the random effects
// @param sigma2_ua variance parameter for the random effects
// @param sd_beta_u standard deviation for the random walk proposal
// @param acc_beta_u acceptance rate for the random effects spline coefficients
// @return new random effects spline coefficients
// [[Rcpp::export]]
Rcpp::List sample_reff_b(arma::vec tau, arma::umat D, arma::uvec cens, 
                         arma::mat beta_u_old, 
                         arma::vec delta_dat, arma::mat B_beta_dat, 
                         arma::mat b_dat_old, arma::mat mu_dat, arma::mat B, 
                         arma::mat P, arma::uvec ind, arma::uvec time, 
                         double sigma2_us, double sigma2_ua, arma::mat sd_beta_u, 
                         arma::mat acc_beta_u){
  
  arma::uvec un_ind = arma::unique(ind);
  unsigned int n_ind = un_ind.n_elem;
  unsigned int J = B.n_cols;
  unsigned int n = tau.n_elem;
  unsigned int d_1 = mu_dat.n_cols;
  
  arma::mat beta_u_prop(n_ind, J, arma::fill::zeros);
  arma::mat B_beta_u_dat(n, d_1, arma::fill::zeros);
  arma::mat B_beta_u_prop_dat, b_dat_prop, temp;
  arma::uvec idx_it, idx_i;

  double l_u, alpha_acc, logpost_prop, logpost_old, logpr_prop, logpr_old;
  for (unsigned int i = 1; i <= n_ind; i++){
    for (unsigned int k = 0; k < J; k++){
      // Propose new value
      beta_u_prop.row(i-1) = beta_u_old.row(i-1);
      beta_u_prop(i - 1,k) = beta_u_old(i - 1,k) + sd_beta_u(i - 1,k) * arma::randn();
      
      if (k == 0){ // only data at t = 1 influence the first coefficient
        idx_it = arma::find( (ind == i) && (time == 1) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k + 1) - beta_u_prop(i - 1,k), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k + 1) - beta_u_old(i - 1,k), 2.0) / sigma2_us);
      }
      else if (k == (J - 1)){ // only data at t = T influence the last coefficient
        idx_it = arma::find( (ind == i) && (time == (J - 1)) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k) - beta_u_prop(i - 1,k - 1), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k) - beta_u_old(i - 1,k - 1), 2.0) / sigma2_us);
      }
      else { // data at t = {k-1,k} influence the kth coefficient
        idx_it = arma::find( (ind == i) && ((time == k) || (time == k+1)) );
        logpr_prop = - 0.5 * ( std::pow(beta_u_prop(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_prop(i - 1,k + 1) - beta_u_prop(i - 1,k), 2.0) / sigma2_us + 
          std::pow(beta_u_prop(i - 1,k) - beta_u_prop(i - 1,k - 1), 2.0) / sigma2_us);
        logpr_old = - 0.5 * ( std::pow(beta_u_old(i - 1,k), 2.0) / sigma2_ua + 
          std::pow(beta_u_old(i - 1,k + 1) - beta_u_old(i - 1,k), 2.0) / sigma2_us + 
          std::pow(beta_u_old(i - 1,k) - beta_u_old(i - 1,k - 1), 2.0) / sigma2_us);
      }
      
      // Compute the implied \mu
      B_beta_u_prop_dat = B.rows(idx_it) * beta_u_prop.row(i-1).t();
      temp = B_beta_dat.rows(idx_it);
      temp.each_col() += B_beta_u_prop_dat;
      b_dat_prop = arma::exp(temp);
      
      
      // Evaluate log posterior
      logpost_prop = log_likelihood(tau.elem(idx_it), 
                                    mu_dat.rows(idx_it), 
                                    b_dat_prop, 
                                    delta_dat.elem(idx_it), 
                                    cens.elem(idx_it),
                                    D.rows(idx_it), true) + 
                                      logpr_prop;
      logpost_old = log_likelihood(tau.elem(idx_it), 
                                   mu_dat.rows(idx_it), 
                                   b_dat_old.rows(idx_it), 
                                   delta_dat.elem(idx_it), 
                                   cens.elem(idx_it),
                                   D.rows(idx_it), true) + 
                                     logpr_old;
      
      alpha_acc = logpost_prop - logpost_old;
      l_u = std::log(arma::randu());
      if (l_u < alpha_acc){
        beta_u_old(i - 1,k) = beta_u_prop(i - 1,k);
        b_dat_old.rows(idx_it) = b_dat_prop;
        acc_beta_u(i - 1,k) += 1;
      }
    }
    idx_i = arma::find(ind == i);
    for (unsigned j = 0; j < d_1; j++){
      B_beta_u_dat.submat(idx_i, arma::regspace<arma::uvec>(j, j)) = B.rows(idx_i) * beta_u_old.row(i-1).t();
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("B_beta_u_dat")=B_beta_u_dat,
    Rcpp::Named("beta_u_old")=beta_u_old, 
    Rcpp::Named("acc_beta_u")=acc_beta_u
  );
}

