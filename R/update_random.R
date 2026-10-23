# Random effects and their variance parameters (supplement §S.5.2 steps 7-8).

# Prior precision of one random-effect curve (article eq. 10), or of a scalar
# random effect for the constant boundary models.
re_precision = function(st, which){
  D = st$design
  ua = st[[paste0("sigma2_", which, "_ua")]]
  us = st[[paste0("sigma2_", which, "_us")]]
  if (which == "mu" || D$spline_boundary) diag(D$J) / ua + D$P / us else matrix(1 / ua)
}

# (7) Random-effect coefficients, one at a time (src/mh_updates.cpp).
update_random_effects = function(st, which, post_burnin){
  D = st$design
  spline = which == "mu" || D$spline_boundary
  if (which == "mu"){
    fix = st$lmu_fix; re = st$lmu_re; other = st$lb_fix + st$lb_re; U = st$U_mu
  } else {
    fix = st$lb_fix; re = st$lb_re; other = st$lmu_fix + st$lmu_re; U = st$U_b
  }
  name = paste0("U_", which)
  res = update_re_cpp(D$tau, st$delta_t, D$dec, D$stim - 1L, fix, re, other,
                      which == "mu", st$LL, U,
                      if (spline) D$B else matrix(1, D$n, 1),
                      if (spline) D$rows_ik_spline else D$rows_ik_scalar,
                      re_precision(st, which), st$mh[[name]]$sd, D$re_mode)
  st[[name]] = res$U
  if (which == "mu") st$lmu_re = res$re else st$lb_re = res$re
  st$LL = res$LL
  record_mh(st, name, res$accepted, post_burnin)
}

# (8) Variance parameters of the random effects, both ~ half-Cauchy(0, 1)
# (article §3.2.2), truncated at sigma2_ua_max (default: none) and
# sigma2_us_max.
# The R curves of every subject share these parameters.
update_re_variances = function(st, which, post_burnin){
  D = st$design
  U = st[[paste0("U_", which)]]
  curves = do.call(rbind, lapply(seq_len(D$R), function(r) matrix(U[, , r], D$n_ind)))
  n_curves = nrow(curves)
  ss_a = sum(curves^2)
  name_a = paste0("sigma2_", which, "_ua")
  name_s = paste0("sigma2_", which, "_us")

  if (which == "b" && !D$spline_boundary){
    target = function(a) -0.5 * n_curves * log(a) - 0.5 * ss_a / a + log_prior_var_halfcauchy(a)
    res = mh_positive(st[[name_a]], target, st$mh[[name_a]]$sd, upper = st$hypers$sigma2_ua_max)
    st[[name_a]] = res$value
    return(record_mh(st, name_a, res$accepted, post_burnin))
  }

  ss_s = sum((curves %*% D$P) * curves)
  # log det(I / a + P / s) via the eigenvalues of P
  log_det = function(a, s) sum(log(1 / a + D$P_eigen / s))
  us = st[[name_s]]
  target_a = function(a) 0.5 * n_curves * log_det(a, us) - 0.5 * ss_a / a + log_prior_var_halfcauchy(a)
  res = mh_positive(st[[name_a]], target_a, st$mh[[name_a]]$sd, upper = st$hypers$sigma2_ua_max)
  st[[name_a]] = res$value
  st = record_mh(st, name_a, res$accepted, post_burnin)

  ua = st[[name_a]]
  target_s = function(s) 0.5 * n_curves * log_det(ua, s) - 0.5 * ss_s / s + log_prior_var_halfcauchy(s)
  res = mh_positive(us, target_s, st$mh[[name_s]]$sd, upper = st$hypers$sigma2_us_max)
  st[[name_s]] = res$value
  record_mh(st, name_s, res$accepted, post_burnin)
}
