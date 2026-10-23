// Cluster-label update for the locally clustered fixed effects (article §3.2.1,
// supplement §S.5.2 step 3): locally informed Hamming ball sampling of the
// factorial HMM z^{(s)} = (z^{(1,s)}, ..., z^{(d0,s)}) of one input stimulus s.
//
// Notation (0-based in the code):
//   chains c = 0..d0-1      one per decision d (chain c == s is the "correct" one)
//   knots k = 0..J-1        spline coefficients; data block t sits between
//                           knots t and t+1 with weights 1/2, 1/2
//   labels 0..Z-1           atoms beta*_{k,h}
//
// Exactness. The sampler alternates
//   v_k ~ q(v_k | z_k, z_{k+1}) = g{l_k(v_k, z_{k+1})} 1{v_k in H(z_k)} / C_k(z_k, z_{k+1})
//   z   ~ p(z | v, rest) via forward filtering / backward sampling,
// where l_k is the log-likelihood of block k. For the z-marginal to remain the
// posterior, the z | v conditional must carry 1 / C_k(z_k, z_{k+1}); and since
// the prior on beta* depends on z (article eq. 9), p(beta* | z) must enter as
// well. Both terms are included below (`log_norm` and `psi`).
//
// Unused labels. Coefficients of labels that no curve uses at a knot enter
// neither the likelihood nor any other prior, so they integrate out and their
// distribution does not affect the posterior (a pseudo-prior; Carlin and Chib,
// 1995). They are drawn close to the coefficients of the labels in use at the
// same knot, a mixture of N(beta*_{k,j}, tau2) over used labels j, so that
// moving a curve into an unused label is a plausible split.
//
// Speed. The block log-likelihood is a sum over accumulators of terms that
// depend only on that accumulator's two labels:
//   l_k(a, b) = sum_c tab_k[c, a_c, b_c].
// One table per block therefore answers every likelihood query by lookup.

#include "likelihood.h"
#include <cstdint>
#include <vector>

using lddmm::contrib;
using lddmm::logsumexp2;

namespace {

inline double log_dnorm(double x, double mean, double var) {
  const double r = x - mean;
  return -0.5 * std::log(2.0 * M_PI * var) - 0.5 * r * r / var;
}

// Hamming ball of radius 1 around `center`. Member 0 is the center; member
// 1 + j (Z - 1) + q changes chain j to the q-th label that differs from
// center[j].
struct Ball {
  std::vector<int> center;
  int Z, size;
  Ball() : Z(0), size(0) {}
  Ball(const std::vector<int>& c, int Z_)
    : center(c), Z(Z_), size(1 + static_cast<int>(c.size()) * (Z_ - 1)) {}
  std::vector<int> member(int m) const {
    std::vector<int> out(center);
    if (m > 0) {
      const int j = (m - 1) / (Z - 1), q = (m - 1) % (Z - 1);
      out[j] = (q < center[j]) ? q : q + 1;
    }
    return out;
  }
};

int sample_log_weights(const std::vector<double>& lw) {
  double mx = R_NegInf;
  for (double x : lw) if (x > mx) mx = x;
  if (!R_FINITE(mx)) Rcpp::stop("label update: no configuration has positive probability");
  std::vector<double> w(lw.size());
  double total = 0.0;
  for (size_t i = 0; i < lw.size(); i++) { w[i] = std::exp(lw[i] - mx); total += w[i]; }
  double u = R::unif_rand() * total;
  for (size_t i = 0; i < w.size(); i++) { u -= w[i]; if (u <= 0.0) return static_cast<int>(i); }
  return static_cast<int>(w.size()) - 1;
}

// Everything the update needs about stimulus s, computed once per call.
struct StimulusModel {
  int s, d0, Z, J;
  bool flexible, informed;
  std::vector<std::vector<double>> tab;  // per block: [(c * Z + a) * Z + b]
  std::vector<std::vector<uint32_t>> other_parents;  // per knot k >= 1: parent mask of each label from other stimuli
  const arma::mat &beta_mu, &beta_b;
  std::vector<uint32_t> other_used;                   // per knot: labels used by other stimuli
  double sigma2_mu, sigma2_b, mu0_mu, s20_mu, mu0_b, s20_b, tau2;
  const arma::mat &logQ_C, &logQ_I;
  const arma::vec &logpi0_C, &logpi0_I;

  StimulusModel(int s_, const arma::imat& z_all, int Z_,
                const arma::mat& beta_mu_, const arma::mat& beta_b_, bool flexible_,
                const Rcpp::List& rows_by_block, const arma::vec& tau,
                const arma::vec& delta_t, const arma::ivec& dec,
                const arma::mat& lmu_re, const arma::mat& lb_fix,
                const arma::mat& lb_re, const arma::mat& logQ_C_,
                const arma::mat& logQ_I_, const arma::vec& logpi0_C_,
                const arma::vec& logpi0_I_, double sigma2_mu_, double sigma2_b_,
                double mu0_mu_, double s20_mu_, double mu0_b_, double s20_b_,
                double tau2_, bool informed_)
    : s(s_ - 1), d0(lmu_re.n_cols), Z(Z_), J(z_all.n_cols),
      flexible(flexible_), informed(informed_), beta_mu(beta_mu_), beta_b(beta_b_),
      sigma2_mu(sigma2_mu_), sigma2_b(sigma2_b_), mu0_mu(mu0_mu_), s20_mu(s20_mu_),
      mu0_b(mu0_b_), s20_b(s20_b_), tau2(tau2_), logQ_C(logQ_C_), logQ_I(logQ_I_),
      logpi0_C(logpi0_C_), logpi0_I(logpi0_I_) {
    if (Z > 32) Rcpp::stop("z_max must be at most 32");
    // Likelihood tables, one per block t (knots t and t+1).
    tab.assign(J - 1, std::vector<double>(d0 * Z * Z, 0.0));
    for (int t = 0; t < J - 1; t++) {
      const Rcpp::IntegerVector rows = rows_by_block[t];
      std::vector<double>& tb = tab[t];
      for (int r : rows) {
        const int i = r - 1;
        for (int c = 0; c < d0; c++) {
          for (int a = 0; a < Z; a++) {
            for (int b = 0; b < Z; b++) {
              const double lmu = 0.5 * (beta_mu(t, a) + beta_mu(t + 1, b)) + lmu_re(i, c);
              const double lb = flexible ?
                0.5 * (beta_b(t, a) + beta_b(t + 1, b)) + lb_re(i, c) :
                lb_fix(i, c) + lb_re(i, c);
              tb[(c * Z + a) * Z + b] += contrib(tau(i), delta_t(i), dec(i), c, lmu, lb);
            }
          }
        }
      }
    }
    // Parent sets and used labels contributed by the chains of the other
    // stimuli (held fixed).
    other_parents.assign(J, std::vector<uint32_t>(Z, 0u));
    other_used.assign(J, 0u);
    for (arma::uword r = 0; r < z_all.n_rows; r++) {
      if (static_cast<int>(r) / d0 == s) continue;
      for (int k = 0; k < J; k++) {
        other_used[k] |= (1u << (z_all(r, k) - 1));
        if (k > 0) other_parents[k][z_all(r, k) - 1] |= (1u << (z_all(r, k - 1) - 1));
      }
    }
  }

  double ell(int t, const std::vector<int>& a, const std::vector<int>& b) const {
    double out = 0.0;
    for (int c = 0; c < d0; c++) out += tab[t][(c * Z + a[c]) * Z + b[c]];
    return out;
  }
  double g(double x) const { return informed ? 0.5 * x : 0.0; }

  // log C_t(a, b) = log sum_{w in H(a)} exp g{l_t(w, b)}; constant when g == 0.
  double log_norm(int t, const std::vector<int>& a, const std::vector<int>& b,
                  double ell_ab) const {
    if (!informed) return 0.0;
    double out = g(ell_ab);
    for (int c = 0; c < d0; c++) {
      const double rest = ell_ab - tab[t][(c * Z + a[c]) * Z + b[c]];
      for (int v = 0; v < Z; v++) {
        if (v == a[c]) continue;
        out = logsumexp2(out, g(rest + tab[t][(c * Z + v) * Z + b[c]]));
      }
    }
    return out;
  }

  // log p(beta*_{k, .} | beta*_{k-1, .}, z) for k >= 1, where stimulus s has
  // labels a at knot k-1 and b at knot k (article eq. 9, normalized reading):
  // a label in use with parent set P is N(mean(beta*_{k-1, P}), sigma2 / |P|),
  // an unused label has the pseudo-prior.
  // log pseudo-prior of the unused label h at knot k, given the set of used
  // labels (bit mask).
  double log_pseudo(int k, int h, uint32_t used) const {
    double acc = R_NegInf, n = 0.0;
    for (int j = 0; j < Z; j++) {
      if (!(used & (1u << j))) continue;
      double v = log_dnorm(beta_mu(k, h), beta_mu(k, j), tau2);
      if (flexible) v += log_dnorm(beta_b(k, h), beta_b(k, j), tau2);
      acc = logsumexp2(acc, v);
      n += 1.0;
    }
    return acc - std::log(n);
  }

  // Knot 1: labels in use have the N(mu0, s20) prior, unused labels the
  // pseudo-prior; stimulus s has labels a.
  double psi_start(const std::vector<int>& a) const {
    uint32_t used = other_used[0];
    for (int c = 0; c < d0; c++) used |= (1u << a[c]);
    double out = 0.0;
    for (int h = 0; h < Z; h++) {
      if (used & (1u << h)) {
        out += log_dnorm(beta_mu(0, h), mu0_mu, s20_mu);
        if (flexible) out += log_dnorm(beta_b(0, h), mu0_b, s20_b);
      } else {
        out += log_pseudo(0, h, used);
      }
    }
    return out;
  }

  double psi(int k, const std::vector<int>& a, const std::vector<int>& b) const {
    std::vector<uint32_t> mask(other_parents[k]);
    for (int c = 0; c < d0; c++) mask[b[c]] |= (1u << a[c]);
    uint32_t used = 0u;
    for (int h = 0; h < Z; h++) if (mask[h] != 0u) used |= (1u << h);
    double out = 0.0;
    for (int h = 0; h < Z; h++) {
      if (mask[h] == 0u) {
        out += log_pseudo(k, h, used);
        continue;
      }
      double m = 0.0, mean_mu = 0.0, mean_b = 0.0;
      for (int j = 0; j < Z; j++) {
        if (mask[h] & (1u << j)) {
          m += 1.0;
          mean_mu += beta_mu(k - 1, j);
          if (flexible) mean_b += beta_b(k - 1, j);
        }
      }
      out += log_dnorm(beta_mu(k, h), mean_mu / m, sigma2_mu / m);
      if (flexible) out += log_dnorm(beta_b(k, h), mean_b / m, sigma2_b / m);
    }
    return out;
  }

  // A simpler smoothness term in place of psi (used only by the tests): squared
  // differences between the two knots, summed over this stimulus's chains,
  // without normalizing constants or the prior of unused labels.
  double psi_simple(int k, const std::vector<int>& a, const std::vector<int>& b) const {
    double out = 0.0;
    for (int c = 0; c < d0; c++) {
      out -= 0.5 / sigma2_mu * std::pow(beta_mu(k - 1, a[c]) - beta_mu(k, b[c]), 2.0);
      if (flexible) out -= 0.5 / sigma2_b * std::pow(beta_b(k - 1, a[c]) - beta_b(k, b[c]), 2.0);
    }
    return out;
  }

  double log_trans(const std::vector<int>& a, const std::vector<int>& b) const {
    double out = 0.0;
    for (int c = 0; c < d0; c++)
      out += (c == s) ? logQ_C(a[c], b[c]) : logQ_I(a[c], b[c]);
    return out;
  }
  double log_init(const std::vector<int>& a) const {
    double out = 0.0;
    for (int c = 0; c < d0; c++) out += (c == s) ? logpi0_C(a[c]) : logpi0_I(a[c]);
    return out;
  }
};

}  // namespace

// Exact log full conditional (up to a constant) of the labels of stimulus s,
// `z_s` (d0 x J, 1-based); used to validate the sampler by enumeration.
// [[Rcpp::export(rng = false)]]
double hb_log_target(int s, const arma::imat& z_all, int Z, const arma::imat& z_s,
                     const arma::mat& beta_mu, const arma::mat& beta_b, bool flexible,
                     const Rcpp::List& rows_by_block, const arma::vec& tau,
                     const arma::vec& delta_t, const arma::ivec& dec,
                     const arma::mat& lmu_re, const arma::mat& lb_fix,
                     const arma::mat& lb_re, const arma::mat& logQ_C,
                     const arma::mat& logQ_I, const arma::vec& logpi0_C,
                     const arma::vec& logpi0_I, double sigma2_mu, double sigma2_b,
                     double mu0_mu, double s20_mu, double mu0_b, double s20_b,
                     double tau2) {
  StimulusModel M(s, z_all, Z, beta_mu, beta_b, flexible, rows_by_block, tau,
                  delta_t, dec, lmu_re, lb_fix, lb_re, logQ_C, logQ_I, logpi0_C,
                  logpi0_I, sigma2_mu, sigma2_b, mu0_mu, s20_mu, mu0_b, s20_b, tau2,
                  false);
  std::vector<std::vector<int>> zz(M.J, std::vector<int>(M.d0));
  for (int k = 0; k < M.J; k++)
    for (int c = 0; c < M.d0; c++) zz[k][c] = z_s(c, k) - 1;
  double out = M.log_init(zz[0]) + M.psi_start(zz[0]);
  for (int k = 0; k < M.J - 1; k++)
    out += M.ell(k, zz[k], zz[k + 1]) + M.psi(k + 1, zz[k], zz[k + 1]) +
      M.log_trans(zz[k], zz[k + 1]);
  return out;
}

// One locally informed Hamming ball update of the labels of stimulus s.
// Returns the new labels as a d0 x J matrix (1-based).
// `drop_norm` (leave out the normalizing constant C_k) and `simple_prior` (use
// psi_simple instead of p(beta* | z)) give approximate kernels that do not
// leave the posterior invariant; the tests use them to check that the
// invariance test can detect this. LDDMM() never sets them.
// [[Rcpp::export]]
arma::imat hb_update_stimulus(int s, const arma::imat& z_all, int Z,
                              const arma::mat& beta_mu, const arma::mat& beta_b,
                              bool flexible, const Rcpp::List& rows_by_block,
                              const arma::vec& tau, const arma::vec& delta_t,
                              const arma::ivec& dec, const arma::mat& lmu_re,
                              const arma::mat& lb_fix, const arma::mat& lb_re,
                              const arma::mat& logQ_C, const arma::mat& logQ_I,
                              const arma::vec& logpi0_C, const arma::vec& logpi0_I,
                              double sigma2_mu, double sigma2_b, double mu0_mu,
                              double s20_mu, double mu0_b, double s20_b,
                              double tau2, bool informed, bool drop_norm = false,
                              bool simple_prior = false) {
  StimulusModel M(s, z_all, Z, beta_mu, beta_b, flexible, rows_by_block, tau,
                  delta_t, dec, lmu_re, lb_fix, lb_re, logQ_C, logQ_I, logpi0_C,
                  logpi0_I, sigma2_mu, sigma2_b, mu0_mu, s20_mu, mu0_b, s20_b, tau2,
                  informed);
  const int J = M.J, d0 = M.d0;
  std::vector<std::vector<int>> cur(J, std::vector<int>(d0));
  for (int k = 0; k < J; k++)
    for (int c = 0; c < d0; c++) cur[k][c] = z_all(M.s * d0 + c, k) - 1;

  // (a) Auxiliary variables: v_k in H(z_k) with weights g{l_k(v_k, z_{k+1})};
  //     v_J uniform on H(z_J).
  std::vector<Ball> V(J);
  for (int k = 0; k < J; k++) {
    Ball H(cur[k], Z);
    std::vector<double> lw(H.size, 0.0);
    if (k < J - 1)
      for (int m = 0; m < H.size; m++) lw[m] = M.g(M.ell(k, H.member(m), cur[k + 1]));
    V[k] = Ball(H.member(sample_log_weights(lw)), Z);
  }

  // (b) Pairwise log potentials phi_k(a, b) on H(v_k) x H(v_{k+1}) and
  //     backward messages.
  const int HB = V[0].size;
  std::vector<arma::mat> phi(J - 1, arma::mat(HB, HB));
  std::vector<arma::vec> msg(J, arma::vec(HB, arma::fill::zeros));
  for (int k = J - 2; k >= 0; k--) {
    std::vector<std::vector<int>> A(HB), B(HB);
    for (int m = 0; m < HB; m++) { A[m] = V[k].member(m); B[m] = V[k + 1].member(m); }
    std::vector<double> g_v(HB);
    for (int ib = 0; ib < HB; ib++) g_v[ib] = M.g(M.ell(k, V[k].center, B[ib]));
    for (int ia = 0; ia < HB; ia++) {
      double acc = R_NegInf;
      for (int ib = 0; ib < HB; ib++) {
        const double l_ab = M.ell(k, A[ia], B[ib]);
        const double val = l_ab + g_v[ib] +
          (drop_norm ? 0.0 : -M.log_norm(k, A[ia], B[ib], l_ab)) +
          (simple_prior ? M.psi_simple(k + 1, A[ia], B[ib]) : M.psi(k + 1, A[ia], B[ib])) +
          M.log_trans(A[ia], B[ib]);
        phi[k](ia, ib) = val;
        acc = logsumexp2(acc, val + msg[k + 1](ib));
      }
      msg[k](ia) = acc;
    }
  }

  // (c) Forward sampling.
  arma::imat out(d0, J);
  std::vector<double> lw(HB);
  for (int m = 0; m < HB; m++) {
    const std::vector<int> a = V[0].member(m);
    lw[m] = msg[0](m) + M.log_init(a) + (simple_prior ? 0.0 : M.psi_start(a));
  }
  int idx = sample_log_weights(lw);
  for (int k = 0; k < J; k++) {
    if (k > 0) {
      for (int m = 0; m < HB; m++) lw[m] = phi[k - 1](idx, m) + msg[k](m);
      idx = sample_log_weights(lw);
    }
    const std::vector<int> z = V[k].member(idx);
    for (int c = 0; c < d0; c++) out(c, k) = z[c] + 1;
  }
  return out;
}
