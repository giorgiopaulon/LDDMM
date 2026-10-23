#' Spline Basis Functions
#'
#' Quadratic B-spline basis on equally spaced knots (supplement §S.1), evaluated
#' on a grid. Points outside `[min(knots), max(knots)]` get a row of zeros.
#'
#' @param xgrid grid where we want to evaluate the spline functions (vector of length n)
#' @param knots vector of equally spaced knots for the splines (vector of length K)
#' @return n x (K+1) matrix with the value of each basis function on `xgrid`.
#'   At a knot, exactly two adjacent basis functions are non-zero, both equal to 1/2.
#' @export
B_basis = function(xgrid, knots){
  n = length(xgrid)
  K = length(knots)

  B = array(0, dim = c(n, K + 1))
  for (j in 1:(K-1)){
    act_idx = (1:n)[(xgrid >= knots[j]) & (xgrid <= knots[j+1])]
    resc_x = (xgrid[act_idx] - knots[j]) / (knots[j+1] - knots[j])

    B[act_idx,j] = (1/2) * (1 - resc_x)^2
    B[act_idx,j+1] = -(resc_x^2) + resc_x + 1/2
    B[act_idx,j+2] = (resc_x^2)/2
  }

  return(B)
}


#' Spline Penalty Matrix
#'
#' Penalty matrix of the smoothness inducing prior for spline coefficients:
#' `t(beta) %*% P %*% beta` is the sum of squared first differences of `beta`
#' (Eilers and Marx, 1996; article eq. 10).
#'
#' @param K number of spline coefficients
#' @return K x K matrix `P = t(D) %*% D`, with `D` the (K-1) x K first-difference matrix.
#' @export
P_smooth1 = function(K){
  D = diff(diag(K))
  return(crossprod(D))
}


#' Hamming Ball
#'
#' All label vectors within Hamming distance `r` of `x` (supplement §S.4).
#'
#' @param x center of the Hamming Ball (vector of labels in `1:S`)
#' @param S number of states
#' @param r radius of the Hamming Ball
#' @return Matrix with `length(x)` rows; each column is one member of the ball,
#'   starting with `x` itself and ordered by distance from `x`.
#' @export
H_ball = function(x, S, r){
  K = length(x)
  out = matrix(x, K, 1)
  for (m in seq_len(min(r, K))){
    positions = utils::combn(seq_len(K), m, simplify = FALSE)
    for (idx in positions){
      alternatives = expand.grid(lapply(idx, function(p) setdiff(seq_len(S), x[p])))
      block = matrix(x, K, nrow(alternatives))
      block[idx, ] = t(as.matrix(alternatives))
      out = cbind(out, block)
    }
  }
  dimnames(out) = NULL
  return(out)
}


# Small numerical helpers used by the sampler -------------------------------

log_dnorm = function(x, mean, var){
  -0.5 * log(2 * pi * var) - (x - mean)^2 / (2 * var)
}

logsumexp = function(x){
  m = max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Log prior density of a variance v whose square root has a half-Cauchy(0, scale)
# prior; the last term is the Jacobian of sigma -> sigma^2.
log_prior_sd_halfcauchy = function(v, scale){
  log(2 * scale / pi) - log(scale^2 + v) - log(2) - 0.5 * log(v)
}

# Log prior density of a variance v ~ half-Cauchy(0, 1) (article §3.2.2).
log_prior_var_halfcauchy = function(v){
  log(2 / pi) - log1p(v^2)
}

# Log of a Dirichlet(alpha) draw. Gamma(a) = Gamma(a + 1) * U^(1/a) keeps tiny
# concentrations from underflowing to exact zeros.
log_rdirichlet = function(alpha){
  lg = log(stats::rgamma(length(alpha), alpha + 1)) + log(stats::runif(length(alpha))) / alpha
  lg - logsumexp(lg)
}

# One Metropolis-Hastings step for a positive scalar with a log-normal random
# walk. `log_target` is the log full conditional; proposals above `upper` are
# rejected (the prior is truncated there).
mh_positive = function(cur, log_target, sd, upper = Inf, lower = 0){
  prop = as.numeric(cur * exp(sd * stats::rnorm(1)))
  if (prop > upper || prop < lower) return(list(value = cur, accepted = 0))
  log_ratio = log_target(prop) - log_target(cur) + log(prop) - log(cur)
  if (log(stats::runif(1)) < log_ratio) list(value = prop, accepted = 1)
  else list(value = cur, accepted = 0)
}

# Draw from N(Lambda^{-1} eta, Lambda^{-1}).
rnorm_canonical = function(Lambda, eta){
  R = chol(Lambda)
  mean = backsolve(R, forwardsolve(t(R), eta))
  as.numeric(mean + backsolve(R, stats::rnorm(length(eta))))
}
