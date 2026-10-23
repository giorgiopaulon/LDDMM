small_data = function(){
  data("data", package = "lddmm", envir = environment())
  data[data$subject %in% sort(unique(data$subject))[1:6], ]
}

test_that("every model runs and its caches stay consistent", {
  dat = small_data()
  for (boundaries in c("flexible", "fixed", "constant", "fixed-constant")){
    for (random_effects in c("correct-incorrect", "shared")){
      for (cluster in c(TRUE, FALSE)){
        set.seed(1)
        fit = LDDMM(dat, list(s_sigma_mu = 0.1, s_sigma_b = 0.1), cluster = cluster,
                    boundaries = boundaries, random_effects = random_effects,
                    Niter = 4, burnin = 2, thin = 1, verbose = FALSE, debug_checks = TRUE)
        expect_equal(dim(fit$post_ind_mu), c(5, 6, 2, 2, 2))
        expect_false(anyNA(fit$post_mean_b))
        expect_true(is.finite(compute_WAIC(fit)))
      }
    }
  }
})

test_that("proposal scales adapt during burn-in without changing parameter shapes", {
  dat = small_data()
  set.seed(1)
  fit = LDDMM(dat, list(s_sigma_mu = 0.1, s_sigma_b = 0.1), boundaries = "fixed",
              Niter = 104, burnin = 102, thin = 1, verbose = FALSE, debug_checks = TRUE)
  expect_length(fit$sigma2_mu_ua, 2)
  expect_null(dim(fit$sigma2_1_b))
  expect_true(all(vapply(fit$acceptance, function(a) all(a >= 0 & a <= 1), logical(1))))
})

test_that("summaries and plots work on a fit", {
  dat = small_data()
  set.seed(1)
  fit = LDDMM(dat, list(s_sigma_mu = 0.1, s_sigma_b = 0.1), Niter = 4, burnin = 2,
              thin = 1, verbose = FALSE, save_z = TRUE)
  expect_equal(dim(fit$Z), c(2, 2, 6, 2))
  means = extract_post_mean(dat, fit, "drift")
  expect_equal(nrow(means), 5 * 6)
  expect_setequal(unique(means$subject), unique(dat$subject))
  expect_equal(nrow(extract_post_draws(dat, fit, "boundary")), 5 * 6 * 2)
  for (p in list(plot_post_pars(dat, fit, "drift"), plot_post_pars(dat, fit, "boundary"),
                 plot_accuracy(dat), plot_RT(dat),
                 plot_predictive(dat, fit, "response"), plot_predictive(dat, fit, "rt"))){
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
})

test_that("invalid data are rejected with an informative error", {
  dat = small_data()
  expect_error(LDDMM(dat[, -6], Niter = 4, burnin = 2), "missing column")
  bad = dat
  bad$s[1] = 0
  expect_error(LDDMM(bad, Niter = 4, burnin = 2), "coded 1, ..., K")
  bad = dat
  bad$d[1] = 5
  expect_error(LDDMM(bad, Niter = 4, burnin = 2), "same codes")
  bad = dat
  bad$block[bad$block == 3] = 7
  expect_error(LDDMM(bad, Niter = 4, burnin = 2), "consecutive")
})
