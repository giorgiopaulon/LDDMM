# Joint moves that improve mixing. Each moves several blocks of parameters
# together along a direction in which the posterior changes little, and is an
# exact Metropolis-Hastings or Gibbs step:
#   * scale_random_effects(): rescale random effects and their variances;
#   * shift_offsets(): shift offsets against drift and boundary levels;
#   * translate(): shift fixed effects against random effects.

# Scale move. Multiplying every random-effect coefficient of a parameter by g
# and its variance parameters by g^2 leaves the quadratic form of the prior
# unchanged, so the random effects and their variances can grow or shrink
# together. One-at-a-time updates do this only very slowly, starting from
# small variances. With log g drawn from a symmetric proposal, the map is a
# valid Metropolis-Hastings move; including its Jacobian, the log acceptance
# ratio is
#   (change in log-likelihood) + (change in log hyperprior) + 2 m log g,
# where m is the number of variance parameters (2 for curves, 1 for scalars).
scale_random_effects = function(st, which, post_burnin){
  D = st$design
  hy = st$hypers
  spline = which == "mu" || D$spline_boundary
  name_a = paste0("sigma2_", which, "_ua")
  name_s = paste0("sigma2_", which, "_us")
  name = paste0("scale_", which)
  log_g = st$mh[[name]]$sd * stats::rnorm(1)
  g = exp(log_g)
  ua = g^2 * st[[name_a]]
  us = g^2 * st[[name_s]]
  if (ua > hy$sigma2_ua_max || (spline && us > hy$sigma2_us_max)){
    return(record_mh(st, name, 0, post_burnin))
  }
  # The random part is recomputed from the scaled coefficients rather than
  # scaled in the cache: scaling the cache would multiply its rounding error by
  # g at every accepted move, a random walk that grows without bound in long runs.
  name_U = paste0("U_", which)
  prop = st
  prop[[name_U]] = g * st[[name_U]]
  re = random_part(prop, which)
  LL = if (which == "mu"){
    ll_matrix(D$tau, st$delta_t, D$dec, st$lmu_fix + re, st$lb_fix + st$lb_re)
  } else {
    ll_matrix(D$tau, st$delta_t, D$dec, st$lmu_fix + st$lmu_re, st$lb_fix + re)
  }
  log_ratio = sum(LL) - sum(st$LL) +
    log_prior_var_halfcauchy(ua) - log_prior_var_halfcauchy(st[[name_a]])
  if (spline){
    log_ratio = log_ratio + log_prior_var_halfcauchy(us) - log_prior_var_halfcauchy(st[[name_s]]) +
      4 * log_g
  } else {
    log_ratio = log_ratio + 2 * log_g
  }
  accepted = is.finite(log_ratio) && log(stats::runif(1)) < log_ratio
  if (accepted){
    st[[name_U]] = prop[[name_U]]
    if (which == "mu") st$lmu_re = re else st$lb_re = re
    st[[name_a]] = ua
    if (spline) st[[name_s]] = us
    st$LL = LL
  }
  record_mh(st, name, as.numeric(accepted), post_burnin)
}


# Offset moves. Lowering the offsets delta while raising drifts and boundaries
# leaves the response-time distributions nearly unchanged. With mean decision
# time m = b / mu, a shift of delta by -eps is matched to first order by
# log mu + eps / (2 m) and log b + 3 eps / (2 m), which keep the mean and the
# variance of tau. One-at-a-time updates move along this direction slowly.
# These moves translate the state along the fixed direction
#   v = (-1 for the offsets, ridge["mu"] for log drifts, ridge["b"] for log boundaries),
# with ridge computed once from the data (st$ridge). They are applied once to
# all offsets and fixed effects together, and once per subject to that
# subject's offsets and random effects. A translation along a fixed direction
# with a symmetric proposal is a valid Metropolis-Hastings move with
# acceptance ratio pi(x + eps v) / pi(x).

# Change in the log prior of the fixed effects when every drift coefficient
# moves by d_mu and every boundary coefficient by d_b. Differences between
# knots, and the pseudo-prior of unused labels, are unchanged; only the
# starting coefficients (and the constant boundaries) carry a prior on their level.
fixed_level_prior_change = function(st, d_mu, d_b){
  D = st$design
  hy = st$hypers
  change = function(x, d, mu0, s20) sum(log_dnorm(x + d, mu0, s20) - log_dnorm(x, mu0, s20))
  used1 = unique(as.vector(st$z[, , 1]))
  out = change(st$beta_mu_star[1, used1], d_mu, hy$mu_beta0_mu, hy$sigma2_beta0_mu)
  b_start = switch(D$boundaries, "flexible" = st$beta_b_star[1, used1], "fixed" = st$beta_b[1, ],
                   as.vector(st$b_const))
  out + change(b_start, d_b, hy$mu_beta0_b, hy$sigma2_beta0_b)
}

shift_offsets = function(st, post_burnin){
  D = st$design
  a = st$ridge[["mu"]]
  c = st$ridge[["b"]]

  # (i) all offsets and the fixed effects
  eps = st$mh$ridge_pop$sd * stats::rnorm(1)
  delta = st$delta - eps
  accepted = 0
  if (all(delta > 0 & delta < D$delta_max)){
    LL = ll_matrix(D$tau, st$delta_t - eps, D$dec, st$lmu_fix + a * eps + st$lmu_re,
                   st$lb_fix + c * eps + st$lb_re)
    log_ratio = sum(LL) - sum(st$LL) + fixed_level_prior_change(st, a * eps, c * eps)
    if (is.finite(log_ratio) && log(stats::runif(1)) < log_ratio){
      st$delta = delta
      st$delta_t = st$delta_t - eps
      st$beta_mu_star = st$beta_mu_star + a * eps
      if (D$flexible) st$beta_b_star = st$beta_b_star + c * eps
      else if (D$boundaries == "fixed") st$beta_b = st$beta_b + c * eps
      else st$b_const = st$b_const + c * eps
      st$lmu_fix = st$lmu_fix + a * eps
      st$lb_fix = st$lb_fix + c * eps
      st$LL = LL
      accepted = 1
    }
  }
  st = record_mh(st, "ridge_pop", accepted, post_burnin)

  # (ii) each subject's offsets and random effects
  prec_mu = re_precision(st, "mu")
  prec_b = re_precision(st, "b")
  # change in -1/2 u' Prec u when every coefficient of u moves by d
  re_prior_change = function(U_i, Prec, d){
    sum(apply(matrix(U_i, ncol(Prec)), 2, function(u) -d * sum(Prec %*% u) - 0.5 * d^2 * sum(Prec)))
  }
  accepted = numeric(D$n_ind)
  for (i in seq_len(D$n_ind)){
    eps = st$mh$ridge_subj$sd[i] * stats::rnorm(1)
    delta_i = st$delta[, i] - eps
    if (any(delta_i <= 0 | delta_i >= D$delta_max[, i])) next
    rows = D$rows_ik_scalar[[i]][[1]]
    cols = rep(seq_len(D$d0), each = length(rows))
    rows_e = rep(rows, D$d0)
    lin = rows_e + (cols - 1) * D$n
    delta_t = st$delta_t
    delta_t[rows] = delta_t[rows] - eps
    ll_new = ll_entries(D$tau, delta_t, D$dec, rows_e, cols,
                        st$lmu_fix[lin] + st$lmu_re[lin] + a * eps,
                        st$lb_fix[lin] + st$lb_re[lin] + c * eps)
    log_ratio = sum(ll_new) - sum(st$LL[lin]) +
      re_prior_change(st$U_mu[i, , ], prec_mu, a * eps) +
      re_prior_change(st$U_b[i, , ], prec_b, c * eps)
    if (is.finite(log_ratio) && log(stats::runif(1)) < log_ratio){
      st$delta[, i] = delta_i
      st$delta_t = delta_t
      st$U_mu[i, , ] = st$U_mu[i, , ] + a * eps
      st$U_b[i, , ] = st$U_b[i, , ] + c * eps
      st$lmu_re[lin] = st$lmu_re[lin] + a * eps
      st$lb_re[lin] = st$lb_re[lin] + c * eps
      st$LL[lin] = ll_new
      accepted[i] = 1
    }
  }
  record_mh(st, "ridge_subj", accepted, post_burnin)
}


# Translation move. Adding c_k to every fixed-effect coefficient at knot k and
# subtracting it from every random-effect coefficient at knot k leaves the
# likelihood unchanged, but fixed and random effects can drift against each
# other along this direction and mix slowly. Drawing c from its exact Gaussian
# full conditional p(c | rest) is a valid MCMC move (a group move; Liu and
# Sabatti, 2000).
#
# Lambda and eta collect the precision and linear term of p(c | rest); the
# helpers below add the prior terms of the fixed effects and of the random
# effects. Translations act on all labels, used or not.
add_start_terms = function(Lambda, eta, k, values, mu0, s20){
  Lambda[k, k] = Lambda[k, k] + length(values) / s20
  eta[k] = eta[k] + sum(mu0 - values) / s20
  list(Lambda = Lambda, eta = eta)
}
# Term (e + c_k - c_{k-1})^2 * w / 2, from a coefficient that is normal around
# its parents' mean; e is the current deviation.
add_difference_terms = function(Lambda, eta, k, e, w){
  idx = c(k - 1, k)
  Lambda[idx, idx] = Lambda[idx, idx] + sum(w) * matrix(c(1, -1, -1, 1), 2)
  eta[k] = eta[k] - sum(w * e)
  eta[k - 1] = eta[k - 1] + sum(w * e)
  list(Lambda = Lambda, eta = eta)
}

# The pseudo-prior of unused labels is a mixture centred on the used labels
# at the same knot, so a common shift leaves it unchanged and it contributes
# nothing here.
translation_terms_clustered = function(beta, z, sigma2, mu0, s20){
  J = nrow(beta)
  terms = add_start_terms(matrix(0, J, J), numeric(J), 1, beta[1, unique(as.vector(z[, , 1]))],
                          mu0, s20)
  for (k in 2:J){
    zk = z[, , k]
    zp = z[, , k - 1]
    used = unique(as.vector(zk))
    for (h in used){
      parents = unique(zp[zk == h])
      terms = add_difference_terms(terms$Lambda, terms$eta, k,
                                   beta[k, h] - mean(beta[k - 1, parents]),
                                   length(parents) / sigma2)
    }
  }
  terms
}

translation_terms_spline = function(beta, sigma2, mu0, s20){
  J = nrow(beta)
  terms = add_start_terms(matrix(0, J, J), numeric(J), 1, beta[1, ], mu0, s20)
  for (k in 2:J){
    terms = add_difference_terms(terms$Lambda, terms$eta, k, beta[k, ] - beta[k - 1, ],
                                 rep(1 / sigma2, ncol(beta)))
  }
  terms
}

translate = function(st, which){
  D = st$design
  hy = st$hypers
  U = st[[paste0("U_", which)]]
  Prec = re_precision(st, which)
  if (which == "mu"){
    terms = translation_terms_clustered(st$beta_mu_star, st$z, st$sigma2_1_mu,
                                        hy$mu_beta0_mu, hy$sigma2_beta0_mu)
  } else if (D$flexible){
    terms = translation_terms_clustered(st$beta_b_star, st$z, st$sigma2_1_b,
                                        hy$mu_beta0_b, hy$sigma2_beta0_b)
  } else if (D$boundaries == "fixed"){
    terms = translation_terms_spline(st$beta_b, st$sigma2_1_b, hy$mu_beta0_b, hy$sigma2_beta0_b)
  } else {
    terms = add_start_terms(matrix(0, 1, 1), 0, 1, as.vector(st$b_const),
                            hy$mu_beta0_b, hy$sigma2_beta0_b)
  }
  # Random effects: each curve u contributes (u - c)' Prec (u - c) / 2.
  n_curves = D$n_ind * D$R
  sum_u = apply(U, 2, sum)
  shift = rnorm_canonical(terms$Lambda + n_curves * Prec, terms$eta + as.numeric(Prec %*% sum_u))

  st[[paste0("U_", which)]] = sweep(U, 2, shift)
  if (which == "mu"){
    st$beta_mu_star = st$beta_mu_star + shift
    by_trial = as.numeric(D$B %*% shift)
    st$lmu_fix = st$lmu_fix + by_trial
    st$lmu_re = st$lmu_re - by_trial
    return(st)
  }
  if (D$flexible) st$beta_b_star = st$beta_b_star + shift
  else if (D$boundaries == "fixed") st$beta_b = st$beta_b + shift
  else st$b_const = st$b_const + shift
  by_trial = if (D$spline_boundary) as.numeric(D$B %*% shift) else shift
  st$lb_fix = st$lb_fix + by_trial
  st$lb_re = st$lb_re - by_trial
  st
}
