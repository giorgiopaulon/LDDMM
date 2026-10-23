// Inverse Gaussian race likelihood shared by every sampler kernel.
//
// Each accumulator is a unit-variance Wiener process with drift mu > 0 that
// starts after an offset delta and stops at boundary b > 0. Its first-passage
// time t = tau - delta is inverse Gaussian with mean b / mu and shape b^2
// (article eq. 1). A trial's log-likelihood is a sum over accumulators
// (article eq. 2):
//   * the chosen accumulator of an uncensored trial contributes log g(t),
//   * every other accumulator, and every accumulator of a censored trial,
//     contributes the log survival function log{1 - G(t)}.
// Keeping the per-accumulator terms separate lets the samplers update only the
// terms that a parameter touches.

#ifndef LDDMM_LIKELIHOOD_H
#define LDDMM_LIKELIHOOD_H

#include <RcppArmadillo.h>

namespace lddmm {

// log(1 - exp(x)) for x < 0, accurate for x near 0 and for x very negative.
inline double log1mexp(double x) {
  return (x > -M_LN2) ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

inline double logsumexp2(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  return (a > b) ? a + std::log1p(std::exp(b - a)) : b + std::log1p(std::exp(a - b));
}

// Log density of the first-passage time at t (t <= 0 has zero density).
inline double ig_logpdf(double t, double mu, double b) {
  if (t <= 0.0) return R_NegInf;
  const double r = b - mu * t;
  return std::log(b) - 0.5 * std::log(2.0 * M_PI) - 1.5 * std::log(t) -
    r * r / (2.0 * t);
}

// Log survival function at t:
//   1 - G(t) = Phi((b - mu t) / sqrt(t)) - exp(2 b mu) Phi(-(b + mu t) / sqrt(t)).
// Both terms are kept on the log scale; computing 1 - Phi(.) directly loses all
// precision once the survival probability drops below ~1e-16.
inline double ig_logsurv(double t, double mu, double b) {
  if (t <= 0.0) return 0.0;
  const double st = std::sqrt(t);
  const double l1 = R::pnorm((b - mu * t) / st, 0.0, 1.0, 1, 1);
  const double l2 = 2.0 * b * mu + R::pnorm(-(b + mu * t) / st, 0.0, 1.0, 1, 1);
  if (l2 >= l1) return R_NegInf;  // survival below double precision
  return l1 + log1mexp(l2 - l1);
}

// Contribution of accumulator `acc` (0-based) to one trial, given the log drift
// and log boundary of that accumulator. `dec` is the 0-based decision, or -1
// for a censored trial.
inline double contrib(double tau, double delta, int dec, int acc,
                      double log_mu, double log_b) {
  const double t = tau - delta;
  const double mu = std::exp(log_mu), b = std::exp(log_b);
  return (dec == acc) ? ig_logpdf(t, mu, b) : ig_logsurv(t, mu, b);
}

}  // namespace lddmm

#endif
