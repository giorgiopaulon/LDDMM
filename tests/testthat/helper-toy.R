# A label-update problem small enough to enumerate: 2 decisions, 3 labels,
# 3 knots (2 blocks), so the labels of one stimulus take 3^6 = 729 values.
toy_label_problem = function(flexible, seed = 10){
  set.seed(seed)
  d0 = 2; Z = 3; J = 3
  grid = expand.grid(rep = 1:3, s = 1:d0, t = 1:(J - 1), i = 1:2)
  n = nrow(grid)
  mu = ifelse(outer(grid$s, 1:d0, "=="), 1.5, 0.8)
  times = 0.1 + matrix(statmod::rinvgauss(n * d0, mean = 1.2 / mu, shape = 1.2^2), n)
  prob = list(
    s = 1L, d0 = d0, Z = Z, J = J,
    z_all = matrix(sample.int(Z, d0 * d0 * J, replace = TRUE), d0 * d0, J),
    beta_mu = matrix(stats::rnorm(J * Z, 0, 0.5), J, Z),
    beta_b = matrix(stats::rnorm(J * Z, log(1.2), 0.3), J, Z),
    flexible = flexible,
    tau = apply(times, 1, min), delta_t = rep(0.1, n),
    dec = as.integer(apply(times, 1, which.min) - 1L),
    lmu_re = matrix(stats::rnorm(n * d0, 0, 0.1), n), lb_fix = matrix(log(1.2), n, d0),
    lb_re = matrix(stats::rnorm(n * d0, 0, 0.1), n),
    logQ_C = t(replicate(Z, lddmm:::log_rdirichlet(rep(1, Z)))),
    logQ_I = t(replicate(Z, lddmm:::log_rdirichlet(rep(1, Z)))),
    logpi0_C = lddmm:::log_rdirichlet(rep(1, Z)), logpi0_I = lddmm:::log_rdirichlet(rep(1, Z)),
    sigma2_mu = 0.3, sigma2_b = 0.3, mu0_mu = 0, s20_mu = 0.5, mu0_b = log(1.2), s20_b = 0.3,
    tau2 = 0.1^2)
  storage.mode(prob$z_all) = "integer"
  prob$rows = lapply(1:(J - 1), function(t) which(grid$s == prob$s & grid$t == t))
  prob$configs = as.matrix(expand.grid(rep(list(1:Z), d0 * J)))  # column-major d0 x J
  prob
}

toy_log_target = function(p, z_s){
  lddmm:::hb_log_target(p$s, p$z_all, p$Z, z_s, p$beta_mu, p$beta_b, p$flexible,
                        p$rows, p$tau, p$delta_t, p$dec, p$lmu_re, p$lb_fix, p$lb_re,
                        p$logQ_C, p$logQ_I, p$logpi0_C, p$logpi0_I, p$sigma2_mu,
                        p$sigma2_b, p$mu0_mu, p$s20_mu, p$mu0_b, p$s20_b, p$tau2)
}

toy_exact = function(p){
  lp = apply(p$configs, 1, function(z) toy_log_target(p, matrix(z, p$d0, p$J)))
  exp(lp - lddmm:::logsumexp(lp))
}

# Starting from N exact posterior draws, apply one label update to each and
# count where the chain lands.
toy_one_step = function(p, probs, N, informed, drop_norm = FALSE, simple_prior = FALSE){
  start = sample.int(nrow(p$configs), N, replace = TRUE, prob = probs)
  key = function(z) sum((as.vector(z) - 1) * p$Z^(seq_along(z) - 1)) + 1
  rows_s = (p$s - 1) * p$d0 + seq_len(p$d0)
  counts = numeric(nrow(p$configs))
  for (j in start){
    z_all = p$z_all
    z_all[rows_s, ] = matrix(p$configs[j, ], p$d0, p$J)
    z_new = lddmm:::hb_update_stimulus(
      p$s, z_all, p$Z, p$beta_mu, p$beta_b, p$flexible, p$rows, p$tau, p$delta_t,
      p$dec, p$lmu_re, p$lb_fix, p$lb_re, p$logQ_C, p$logQ_I, p$logpi0_C, p$logpi0_I,
      p$sigma2_mu, p$sigma2_b, p$mu0_mu, p$s20_mu, p$mu0_b, p$s20_b, p$tau2, informed,
      drop_norm, simple_prior)
    k = key(z_new)
    counts[k] = counts[k] + 1
  }
  counts
}

# Chi-square goodness of fit, pooling configurations with small expected counts.
chisq_pvalue = function(counts, probs){
  expected = sum(counts) * probs
  small = expected < 5
  obs = c(counts[!small], sum(counts[small]))
  exp_ = c(expected[!small], sum(expected[small]))
  keep = exp_ > 0
  stat = sum((obs[keep] - exp_[keep])^2 / exp_[keep])
  stats::pchisq(stat, df = sum(keep) - 1, lower.tail = FALSE)
}
