# One-dimensional Metropolis steps: starting from exact posterior draws, one
# step must leave the distribution unchanged (Kolmogorov-Smirnov test).

test_that("the offset update targets its full conditional", {
  set.seed(3)
  n = 15
  tau = 0.3 + stats::rexp(n, 2)
  dec = sample(0:1, n, TRUE)
  log_mu = matrix(log(c(1.5, 0.7)), n, 2, byrow = TRUE)
  log_b = matrix(log(1.2), n, 2)
  dmax = min(tau)
  # exact full conditional on a fine grid (uniform prior on (0, dmax))
  grid = seq(dmax / 4000, dmax * (1 - 1 / 4000), length.out = 4000)
  lp = sapply(grid, function(x) sum(lddmm:::ll_matrix(tau, rep(x, n), dec, log_mu, log_b)))
  cdf = cumsum(exp(lp - max(lp)))
  cdf = cdf / cdf[length(cdf)]
  start = stats::approx(cdf, grid, stats::runif(20000), rule = 2, ties = "ordered")$y
  after = vapply(start, function(x){
    LL = lddmm:::ll_matrix(tau, rep(x, n), dec, log_mu, log_b)
    lddmm:::update_delta_cpp(tau, dec, log_mu, log_b, LL, matrix(x), list(list(seq_len(n))),
                             matrix(dmax), matrix(0.3))$delta[1, 1]
  }, numeric(1))
  ks = suppressWarnings(stats::ks.test(after, stats::approxfun(grid, cdf, rule = 2)))
  expect_gt(ks$p.value, 1e-3)
})

test_that("truncated positive updates target the truncated distribution", {
  set.seed(4)
  upper = 2
  target = function(v) stats::dgamma(v, 3, 2, log = TRUE)
  start = stats::qgamma(stats::runif(20000) * stats::pgamma(upper, 3, 2), 3, 2)
  after = vapply(start, function(v) lddmm:::mh_positive(v, target, 0.8, upper = upper)$value,
                 numeric(1))
  ks = suppressWarnings(stats::ks.test(after, function(q) stats::pgamma(pmin(q, upper), 3, 2) /
                                          stats::pgamma(upper, 3, 2)))
  expect_gt(ks$p.value, 1e-3)
})
