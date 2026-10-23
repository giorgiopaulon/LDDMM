# Simulate-and-recover study in the style of supplement §S.9.
#
# Data are simulated from the model's likelihood with known parameters on the
# design of the tone-learning study: 20 participants, 10 blocks, 4 input tones,
# 40 trials per block. The truth has local clustering of the correct-response
# drifts: tones {1, 3} and {2, 4} form two groups that coincide in the first
# blocks, separate in the middle of training and merge again at the end. The
# fit is then compared with the truth.
#
# Usage:
#   Rscript tests/validation/recovery.R <seed> <Niter> <burnin> <thin> <random_effects> [mcmc_seed]
# <seed> sets the simulated data; the optional <mcmc_seed> reruns the sampler
# on the same data with a different random number stream. <random_effects> is
# the model used both to simulate and to fit (versions before 1.0.0 only
# support "shared"). Environment variables: LDDMM_LIB selects the library holding the
# lddmm version to test (default: the usual library path), LDDMM_OUT the output
# directory (default: validation_out).

args = commandArgs(trailingOnly = TRUE)
seed = as.integer(args[1])
Niter = as.integer(args[2])
burnin = as.integer(args[3])
thin = as.integer(args[4])
random_effects = args[5]
mcmc_seed = if (length(args) > 5) as.integer(args[6]) else NA
lib = Sys.getenv("LDDMM_LIB", unset = "")
lib = if (nzchar(lib)) lib else NULL
suppressMessages(library(lddmm, lib.loc = lib))
version = utils::packageVersion("lddmm", lib.loc = lib)

n_ind = 20; n_blocks = 10; d0 = 4; reps = 10
J = n_blocks + 1
knots = seq_len(n_blocks)

# True spline coefficients (log scale) of every (s, d) curve at the J knots.
k = seq_len(J)
drift_A = log(0.9 + 1.5 * (k - 1) / (J - 1))                    # tones 1 and 3
gap = ifelse(k >= 4 & k <= 9, 0.45 * sin(pi * (k - 3) / 7), 0)  # groups apart at knots 4-9
drift_B = drift_A - gap                                          # tones 2 and 4
drift_wrong = log(0.45 - 0.15 * (k - 1) / (J - 1))
beta_mu = array(NA_real_, c(J, d0, d0))
beta_b = array(NA_real_, c(J, d0, d0))
for (s in seq_len(d0)){
  for (d in seq_len(d0)){
    beta_mu[, s, d] = if (s != d) drift_wrong else if (s %in% c(1, 3)) drift_A else drift_B
    beta_b[, s, d] = if (s == d) rep(log(1.6), J) else rep(log(1.9), J)
  }
}
# True co-clustering of the correct-response curves at each knot (pairs of tones).
pairs = utils::combn(d0, 2)
true_coclust = sapply(k, function(kk) apply(pairs, 2, function(p) {
  isTRUE(all.equal(beta_mu[kk, p[1], p[1]], beta_mu[kk, p[2], p[2]]))
}))

# Random effects: smooth subject curves (article eq. 10), one for correct and
# one for incorrect responses (or a single shared curve).
R = if (random_effects == "shared") 1 else 2
P = lddmm::P_smooth1(J)
draw_curves = function(sigma2_a, sigma2_s){
  L = chol(diag(J) / sigma2_a + P / sigma2_s)
  array(t(backsolve(L, matrix(stats::rnorm(J * n_ind * R), J))), c(n_ind, J, R))
}
re_var = function(sigma2_a, sigma2_s){
  B = lddmm::B_basis(knots, knots)
  rowSums((B %*% solve(diag(J) / sigma2_a + P / sigma2_s)) * B)
}
truth = list(sigma2_mu_ua = 0.04, sigma2_mu_us = 0.01, sigma2_b_ua = 0.02, sigma2_b_us = 0.005)

set.seed(seed)
U_mu = draw_curves(truth$sigma2_mu_ua, truth$sigma2_mu_us)
U_b = draw_curves(truth$sigma2_b_ua, truth$sigma2_b_us)
delta = matrix(stats::runif(d0 * n_ind, 0.15, 0.35), d0, n_ind)

dat = expand.grid(rep = seq_len(reps), s = seq_len(d0), block = knots, subject = seq_len(n_ind))
n = nrow(dat)
B = lddmm::B_basis(dat$block, knots)
# value of random-effect curve r[t] of the subject of every trial t
re_value = function(U, r){
  idx = cbind(rep(dat$subject, J), rep(k, each = n), rep(r, J))
  rowSums(B * matrix(U[idx], n))
}
log_mu = log_b = matrix(NA_real_, n, d0)
for (d in seq_len(d0)){
  r = if (R == 1) rep(1, n) else ifelse(dat$s == d, 1, 2)
  log_mu[, d] = rowSums(B * t(beta_mu[, dat$s, d])) + re_value(U_mu, r)
  log_b[, d] = rowSums(B * t(beta_b[, dat$s, d])) + re_value(U_b, r)
}
off = delta[cbind(dat$s, dat$subject)]
times = off + matrix(statmod::rinvgauss(n * d0, mean = exp(log_b - log_mu), shape = exp(2 * log_b)), n)
dat$d = max.col(-times, ties.method = "first")
dat$r_time = apply(times, 1, min)
dat$cens = 0
dat = dat[, c("subject", "block", "s", "d", "r_time", "cens")]

# Population-level truth on the block grid: exp(f(t) + var(u(t)) / 2).
Bg = lddmm::B_basis(knots, knots)
pop = function(beta, v) exp(array(Bg %*% matrix(beta, J), c(n_blocks, d0, d0)) + v / 2)
true_mu = pop(beta_mu, re_var(truth$sigma2_mu_ua, truth$sigma2_mu_us))
true_b = pop(beta_b, re_var(truth$sigma2_b_ua, truth$sigma2_b_us))
ind_truth = function(beta, U){
  out = array(NA_real_, c(n_blocks, n_ind, d0, d0))
  for (s in seq_len(d0)) for (d in seq_len(d0)){
    r = if (R == 1 || s == d) 1 else 2
    out[, , s, d] = exp(as.numeric(Bg %*% beta[, s, d]) + Bg %*% t(U[, , r]))
  }
  out
}
true_ind_mu = ind_truth(beta_mu, U_mu)
true_ind_b = ind_truth(beta_b, U_b)

fit_args = list(data = dat, hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1),
                Niter = Niter, burnin = burnin, thin = thin)
if (version >= "1.0.0"){  # arguments added in 1.0.0
  fit_args = c(fit_args, list(random_effects = random_effects, save_z = TRUE,
                              save_loglik = FALSE, verbose = FALSE))
}
if (!is.na(mcmc_seed)) set.seed(mcmc_seed)
t0 = proc.time()[["elapsed"]]
fit = do.call(lddmm::LDDMM, fit_args)
elapsed = proc.time()[["elapsed"]] - t0

interval_summary = function(draws, truth, margin){
  lo = apply(draws, margin, stats::quantile, 0.05)
  hi = apply(draws, margin, stats::quantile, 0.95)
  est = apply(draws, margin, mean)
  c(coverage90 = mean(truth >= lo & truth <= hi),
    rmse_log = sqrt(mean((log(est) - log(truth))^2)),
    width_log = mean(log(hi) - log(lo)))
}
res = list(
  seed = seed, version = as.character(version), random_effects = random_effects, seconds = elapsed,
  seconds_per_iter = elapsed / Niter,
  pop_mu = interval_summary(fit$post_mean_mu, true_mu, 1:3),
  pop_b = interval_summary(fit$post_mean_b, true_b, 1:3),
  ind_mu = interval_summary(fit$post_ind_mu, true_ind_mu, 1:4),
  ind_b = interval_summary(fit$post_ind_b, true_ind_b, 1:4),
  delta = c(coverage90 = mean(delta >= apply(fit$post_ind_delta, 1:2, stats::quantile, 0.05) &
                                delta <= apply(fit$post_ind_delta, 1:2, stats::quantile, 0.95)),
            rmse = sqrt(mean((apply(fit$post_ind_delta, 1:2, mean) - delta)^2))),
  est_mu = apply(fit$post_mean_mu, 1:3, mean), true_mu = true_mu,
  est_b = apply(fit$post_mean_b, 1:3, mean), true_b = true_b,
  sigma2 = rbind(truth = unlist(truth),
                 posterior_mean = c(mean(fit$sigma2_mu_ua), mean(fit$sigma2_mu_us),
                                    mean(fit$sigma2_b_ua), mean(fit$sigma2_b_us))),
  true_coclust = true_coclust, mcmc_seed = mcmc_seed,
  trace_b_11 = fit$post_mean_b[1, 1, 1, ], trace_mu_11 = fit$post_mean_mu[1, 1, 1, ],
  trace_delta = apply(fit$post_ind_delta, 3, mean),
  trace_sigma2_mu_ua = fit$sigma2_mu_ua, true_delta_mean = mean(delta),
  est_delta = apply(fit$post_ind_delta, 1:2, mean), true_delta = delta)
if (!is.null(fit$Z)){
  res$post_coclust = sapply(k, function(kk) apply(pairs, 2, function(p) {
    mean(fit$Z[p[1], p[1], kk, ] == fit$Z[p[2], p[2], kk, ])
  }))
}
out_dir = Sys.getenv("LDDMM_OUT", unset = "validation_out")
dir.create(out_dir, showWarnings = FALSE)
tag = if (is.na(mcmc_seed)) "" else sprintf("_mcmc%d", mcmc_seed)
saveRDS(res, file.path(out_dir, sprintf("recovery_%s_%s_seed%03d%s.rds", as.character(version),
                                        random_effects, seed, tag)))
cat(sprintf("seed %d, lddmm %s, %s random effects: %.0fs (%.3f s/iter)\n", seed, as.character(version),
            random_effects, elapsed, elapsed / Niter))
print(rbind(pop_mu = res$pop_mu, pop_b = res$pop_b, ind_mu = res$ind_mu, ind_b = res$ind_b))
print(res$delta)
if (!is.null(res$post_coclust)){
  cat("posterior co-clustering of correct-response drifts (rows: tone pairs, columns: knots)\n")
  print(round(res$post_coclust, 2))
}
