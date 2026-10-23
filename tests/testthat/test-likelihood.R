test_that("inverse Gaussian terms match statmod far into the tails", {
  skip_if_not_installed("statmod")
  g = expand.grid(t = c(0.01, 0.1, 0.5, 1, 3, 10, 30), mu = c(0.2, 1, 3), b = c(0.5, 1.5, 4))
  n = nrow(g)
  log_mu = matrix(log(g$mu))
  log_b = matrix(log(g$b))
  # accumulator 0 wins (density) / accumulator 1 wins (survival of accumulator 0)
  dens = lddmm:::ll_matrix(g$t, rep(0, n), rep(0L, n), log_mu, log_b)[, 1]
  surv = lddmm:::ll_matrix(g$t, rep(0, n), rep(1L, n), log_mu, log_b)[, 1]
  expect_equal(dens, statmod::dinvgauss(g$t, g$b / g$mu, shape = g$b^2, log = TRUE),
               tolerance = 1e-10)
  expect_equal(surv, statmod::pinvgauss(g$t, g$b / g$mu, shape = g$b^2, lower.tail = FALSE,
                                        log.p = TRUE), tolerance = 1e-8)
  expect_true(all(is.finite(surv)))
})

test_that("trial log-likelihood is the race likelihood, including censoring", {
  skip_if_not_installed("statmod")
  set.seed(1)
  n = 60
  d1 = 3
  tau = stats::rexp(n) + 0.2
  delta = stats::runif(n, 0, 0.1)
  mu = matrix(stats::rexp(n * d1) + 0.2, n)
  b = matrix(stats::rexp(n * d1) + 0.5, n)
  cens = stats::rbinom(n, 1, 0.2)
  D = cbind(sample.int(d1, n, TRUE), sample.int(d1, n, TRUE))
  D[cens == 1, 2] = NA
  lg = function(t, a) statmod::dinvgauss(tau[t] - delta[t], b[t, a] / mu[t, a], shape = b[t, a]^2, log = TRUE)
  lS = function(t, a) statmod::pinvgauss(tau[t] - delta[t], b[t, a] / mu[t, a], shape = b[t, a]^2,
                                         lower.tail = FALSE, log.p = TRUE)
  manual = sapply(seq_len(n), function(t){
    if (cens[t] == 1) return(sum(sapply(seq_len(d1), function(a) lS(t, a))))
    won = D[t, 2]
    lg(t, won) + sum(sapply(setdiff(seq_len(d1), won), function(a) lS(t, a)))
  })
  expect_equal(log_likelihood_ind(tau, mu, b, delta, cens, D), manual, tolerance = 1e-8)
  expect_equal(log_likelihood(tau, mu, b, delta, cens, D, TRUE), sum(manual), tolerance = 1e-8)
  expect_equal(log_likelihood(tau[1:3], mu[1:3, ], b[1:3, ], delta[1:3], cens[1:3], D[1:3, ], FALSE),
               exp(sum(manual[1:3])), tolerance = 1e-8)
})

test_that("a response before the offset has zero likelihood", {
  ll = log_likelihood_ind(c(0.5, 0.05), matrix(1, 2, 2), matrix(1, 2, 2), c(0.1, 0.1),
                          c(0, 0), cbind(c(1, 1), c(1, 2)))
  expect_true(is.finite(ll[1]))
  expect_equal(ll[2], -Inf)
})
