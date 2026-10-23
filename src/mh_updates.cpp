// Metropolis-Hastings kernels for the offsets delta and for the random effects.
//
// The sampler keeps a cache LL (n x d0) of per-(trial, accumulator)
// log-likelihood terms. A proposal recomputes only the entries it touches; the
// current value of those entries is read from the cache. Updated caches are
// returned to R rather than modified in place.

#include "likelihood.h"

using lddmm::contrib;

// Offsets delta_s^{(i)} (supplement §S.5.2 step 1). Prior Unif(0, delta_max)
// with delta_max the smallest uncensored response time of subject i to
// stimulus s; log-normal random-walk proposal. Proposals outside the support
// are rejected (not redrawn), which keeps the proposal symmetric on the log
// scale.
// [[Rcpp::export]]
Rcpp::List update_delta_cpp(const arma::vec& tau, const arma::ivec& dec,
                            const arma::mat& log_mu, const arma::mat& log_b,
                            arma::mat LL, arma::mat delta,
                            const Rcpp::List& rows_si, const arma::mat& delta_max,
                            const arma::mat& sd) {
  const int n_stim = delta.n_rows, n_ind = delta.n_cols, d0 = LL.n_cols;
  arma::mat accepted(n_stim, n_ind, arma::fill::zeros);
  for (int s = 0; s < n_stim; s++) {
    const Rcpp::List rows_s = rows_si[s];
    for (int i = 0; i < n_ind; i++) {
      const Rcpp::IntegerVector rows = rows_s[i];
      const double cur = delta(s, i);
      const double prop = cur * std::exp(sd(s, i) * R::norm_rand());
      if (prop >= delta_max(s, i)) continue;
      arma::mat buf(rows.size(), d0);
      double diff = std::log(prop) - std::log(cur);  // log-normal proposal Jacobian
      for (int j = 0; j < rows.size(); j++) {
        const int t = rows[j] - 1;
        for (int c = 0; c < d0; c++) {
          buf(j, c) = contrib(tau(t), prop, dec(t), c, log_mu(t, c), log_b(t, c));
          diff += buf(j, c) - LL(t, c);
        }
      }
      if (std::log(R::unif_rand()) < diff) {
        delta(s, i) = prop;
        for (int j = 0; j < rows.size(); j++) LL.row(rows[j] - 1) = buf.row(j);
        accepted(s, i) = 1.0;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("delta") = delta, Rcpp::Named("LL") = LL,
                            Rcpp::Named("accepted") = accepted);
}

// Random effects (article eq. 10; supplement §S.5.2 step 7), one coefficient at
// a time with a normal random walk.
//
// U is n_ind x K x R: K basis coefficients (K = J for spline curves, K = 1 for
// the scalar random effects of constant boundaries) for R curves per subject.
// With mode == 1 (correct/incorrect) curve 0 acts on the accumulator of the
// correct decision (column == stimulus) and curve 1 on the others; with
// mode == 0 the single curve acts on every accumulator.
//
// `re` is the n x d0 cache of the random-effect part of this parameter on the
// log scale, `fix` its fixed-effect part, and `log_other` the full log value of
// the other parameter (boundary when updating drifts and vice versa).
// rows_ik[i][k] lists the (1-based) trials of subject i where basis function k
// is non-zero. Prec is the K x K prior precision of each curve.
// [[Rcpp::export]]
Rcpp::List update_re_cpp(const arma::vec& tau, const arma::vec& delta_t,
                         const arma::ivec& dec, const arma::ivec& stim,
                         const arma::mat& fix, arma::mat re,
                         const arma::mat& log_other, bool is_drift, arma::mat LL,
                         arma::cube U, const arma::mat& Bu,
                         const Rcpp::List& rows_ik, const arma::mat& Prec,
                         const arma::cube& sd, int mode) {
  const int n_ind = U.n_rows, K = U.n_cols, R = U.n_slices, d0 = re.n_cols;
  arma::cube accepted(n_ind, K, R, arma::fill::zeros);
  std::vector<double> buf;
  for (int i = 0; i < n_ind; i++) {
    const Rcpp::List rows_i = rows_ik[i];
    for (int r = 0; r < R; r++) {
      for (int k = 0; k < K; k++) {
        const Rcpp::IntegerVector rows = rows_i[k];
        const double cur = U(i, k, r);
        const double prop = cur + sd(i, k, r) * R::norm_rand();
        const double du = prop - cur;
        // Prior: -1/2 u' Prec u, restricted to the terms involving u_k.
        double cross = 0.0;
        for (int j = 0; j < K; j++) if (j != k) cross += Prec(k, j) * U(i, j, r);
        double diff = -0.5 * Prec(k, k) * (prop * prop - cur * cur) - du * cross;
        buf.clear();
        for (int j = 0; j < rows.size(); j++) {
          const int t = rows[j] - 1;
          const double shift = Bu(t, k) * du;
          for (int c = 0; c < d0; c++) {
            if (mode == 1 && ((c == stim(t)) != (r == 0))) continue;
            const double lthis = fix(t, c) + re(t, c) + shift;
            const double val = is_drift ?
              contrib(tau(t), delta_t(t), dec(t), c, lthis, log_other(t, c)) :
              contrib(tau(t), delta_t(t), dec(t), c, log_other(t, c), lthis);
            buf.push_back(val);
            diff += val - LL(t, c);
          }
        }
        if (std::log(R::unif_rand()) < diff) {
          U(i, k, r) = prop;
          size_t pos = 0;
          for (int j = 0; j < rows.size(); j++) {
            const int t = rows[j] - 1;
            const double shift = Bu(t, k) * du;
            for (int c = 0; c < d0; c++) {
              if (mode == 1 && ((c == stim(t)) != (r == 0))) continue;
              re(t, c) += shift;
              LL(t, c) = buf[pos++];
            }
          }
          accepted(i, k, r) = 1.0;
        }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("re") = re,
                            Rcpp::Named("LL") = LL,
                            Rcpp::Named("accepted") = accepted);
}

// Inverse Gaussian draws (Michael, Schucany and Haas, 1976) with the given
// means and shapes, used for posterior predictive draws.
// [[Rcpp::export]]
arma::vec rinvgauss_cpp(const arma::vec& mean, const arma::vec& shape) {
  const arma::uword n = mean.n_elem;
  arma::vec out(n);
  for (arma::uword j = 0; j < n; j++) {
    const double m = mean(j), l = shape(j);
    const double y = std::pow(R::norm_rand(), 2.0);
    const double x = m + m * m * y / (2.0 * l) -
      m / (2.0 * l) * std::sqrt(4.0 * m * l * y + m * m * y * y);
    out(j) = (R::unif_rand() <= m / (m + x)) ? x : m * m / x;
  }
  return out;
}
