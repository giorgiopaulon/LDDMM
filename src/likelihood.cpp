// Exported log-likelihood functions and the per-(trial, accumulator) kernels
// used by the samplers to build and update their likelihood cache.

#include "likelihood.h"

using lddmm::contrib;

// Decision codes used internally: 0-based decision, or -1 for censored trials.
static arma::ivec decision_codes(const arma::vec& cens, const arma::mat& D) {
  const arma::uword n = cens.n_elem;
  if (D.n_rows != n || D.n_cols < 2)
    Rcpp::stop("D must be an n x 2 matrix (stimulus, decision)");
  arma::ivec dec(n);
  for (arma::uword t = 0; t < n; t++) {
    dec(t) = (cens(t) == 1) ? -1 : static_cast<int>(D(t, 1)) - 1;
  }
  return dec;
}

//' Log-likelihood of each observation
//'
//' Log-likelihood of the inverse Gaussian race model for every observation,
//' including the contribution of censored observations.
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
//' @param D (n x 2) matrix whose first column has the n input stimuli, and whose
//' second column has the n decision categories
//' @return Vector of n log-likelihood values. An uncensored response time that
//' is not larger than its offset has log-likelihood \code{-Inf}.
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector log_likelihood_ind(const arma::vec& tau, const arma::mat& mu,
                                       const arma::mat& b, const arma::vec& delta,
                                       const arma::vec& cens, const arma::mat& D) {
  const arma::ivec dec = decision_codes(cens, D);
  const arma::uword n = tau.n_elem, d1 = mu.n_cols;
  Rcpp::NumericVector out(n);
  for (arma::uword t = 0; t < n; t++) {
    for (arma::uword a = 0; a < d1; a++) {
      out(t) += contrib(tau(t), delta(t), dec(t), a,
                        std::log(mu(t, a)), std::log(b(t, a)));
    }
  }
  return out;
}

//' Log-likelihood
//'
//' Log-likelihood of the inverse Gaussian race model (article eq. 2), including
//' the contribution of censored observations.
//'
//' @inheritParams log_likelihood_ind
//' @param log should the result be returned on the log scale?
//' @return The (log-)likelihood of all n observations.
//' @export
// [[Rcpp::export(rng = false)]]
double log_likelihood(const arma::vec& tau, const arma::mat& mu,
                      const arma::mat& b, const arma::vec& delta,
                      const arma::vec& cens, const arma::mat& D, bool log) {
  const double out = Rcpp::sum(log_likelihood_ind(tau, mu, b, delta, cens, D));
  return log ? out : std::exp(out);
}

// Matrix of per-(trial, accumulator) contributions; used to (re)build the
// sampler's likelihood cache.
// [[Rcpp::export(rng = false)]]
arma::mat ll_matrix(const arma::vec& tau, const arma::vec& delta,
                    const arma::ivec& dec, const arma::mat& log_mu,
                    const arma::mat& log_b) {
  const arma::uword n = tau.n_elem, d0 = log_mu.n_cols;
  arma::mat out(n, d0);
  for (arma::uword a = 0; a < d0; a++) {
    for (arma::uword t = 0; t < n; t++) {
      out(t, a) = contrib(tau(t), delta(t), dec(t), a, log_mu(t, a), log_b(t, a));
    }
  }
  return out;
}

// Contributions of selected (trial, accumulator) entries at proposed values.
// `rows` and `cols` are 1-based, as produced in R.
// [[Rcpp::export(rng = false)]]
arma::vec ll_entries(const arma::vec& tau, const arma::vec& delta,
                     const arma::ivec& dec, const arma::uvec& rows,
                     const arma::uvec& cols, const arma::vec& log_mu,
                     const arma::vec& log_b) {
  const arma::uword m = rows.n_elem;
  arma::vec out(m);
  for (arma::uword j = 0; j < m; j++) {
    const arma::uword t = rows(j) - 1;
    out(j) = contrib(tau(t), delta(t), dec(t), cols(j) - 1, log_mu(j), log_b(j));
  }
  return out;
}
