# Geweke (2004) joint distribution test of the full sampler.
#
# The successive-conditional simulator alternates one sampler sweep with a
# fresh draw of the data from the likelihood. If every update leaves its full
# conditional invariant, the parameters it visits are distributed as the prior;
# any error in a conditional shows up as a mismatch with direct prior draws.
#
# Usage:
#   Rscript tests/validation/geweke.R <boundaries> <random_effects> <cluster> <iterations> [seed]
# The environment variable LDDMM_LIB selects the library holding lddmm
# (default: the usual library path).

args = commandArgs(trailingOnly = TRUE)
boundaries = if (length(args) > 0) args[1] else "flexible"
random_effects = if (length(args) > 1) args[2] else "correct-incorrect"
cluster = if (length(args) > 2) as.logical(args[3]) else TRUE
n_iter = if (length(args) > 3) as.integer(args[4]) else 20000
seed = if (length(args) > 4) as.integer(args[5]) else 2024

lib = Sys.getenv("LDDMM_LIB", unset = "")
suppressMessages(library(lddmm, lib.loc = if (nzchar(lib)) lib else NULL))
script_dir = dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "simulate.R"))

hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1, sigma2_1_max = 0.2, sigma2_us_max = 0.2,
              sigma2_ua_max = 0.3, mu_beta0_mu = 0, sigma2_beta0_mu = 0.1,
              mu_beta0_b = log(1.5), sigma2_beta0_b = 0.1)
delta_upper = 0.2

# Label-invariant summaries of the state.
summaries = function(st){
  D = st$design
  coef_mu = lddmm:::label_coefs(st$beta_mu_star, st$z)
  out = c(delta_11 = st$delta[1, 1], delta_22 = st$delta[2, 2],
          mu_coef_111 = coef_mu[1, 1, 1], mu_coef_221 = coef_mu[D$J, 2, 1],
          sigma2_1_mu = st$sigma2_1_mu,
          sigma2_mu_ua = st$sigma2_mu_ua, sigma2_mu_us = st$sigma2_mu_us,
          sigma2_b_ua = st$sigma2_b_ua,
          u_mu_12C = st$U_mu[1, 2, 1], u_mu_1JI = st$U_mu[1, D$J, D$R],
          u_b_11 = st$U_b[1, 1, D$R])
  if (D$spline_boundary) out["sigma2_b_us"] = st$sigma2_b_us
  if (D$flexible){
    coef_b = lddmm:::label_coefs(st$beta_b_star, st$z)
    out = c(out, b_coef_111 = coef_b[1, 1, 1], b_coef_212 = coef_b[2, 1, 2],
            sigma2_1_b = st$sigma2_1_b)
  } else if (D$boundaries == "fixed"){
    out = c(out, b_coef_11 = st$beta_b[1, 1], b_coef_J2 = st$beta_b[D$J, 2],
            sigma2_1_b = st$sigma2_1_b)
  } else {
    out = c(out, b_const_1 = st$b_const[1])
  }
  if (D$cluster){
    out = c(out, alpha_C = st$alpha_C, alpha_I = st$alpha_I,
            same_label_k2 = as.numeric(st$z[1, 1, 2] == st$z[2, 2, 2]),
            n_labels_k2 = length(unique(as.vector(st$z[, , 2]))))
  }
  out
}

set.seed(seed)
st = toy_state(n_ind = 3, n_blocks = 2, d0 = 2, reps = 3, boundaries, random_effects, cluster, hypers)

# Marginal-conditional simulator: independent prior draws.
prior = t(replicate(n_iter, summaries(prior_draw(st, delta_upper))))

# Successive-conditional simulator. Proposal scales are tuned in a first phase
# and then frozen.
st = simulate_responses(prior_draw(st, delta_upper), delta_upper)
n_tune = 2500
chain = matrix(NA_real_, n_iter, ncol(prior), dimnames = list(NULL, colnames(prior)))
t0 = proc.time()[["elapsed"]]
for (iter in seq_len(n_tune + n_iter)){
  st = lddmm:::sampler_sweep(st, post_burnin = FALSE)
  if (iter <= n_tune && iter %% 50 == 0) st = lddmm:::adapt_scales(st, iter / 50)
  st = simulate_responses(st, delta_upper)
  if (iter > n_tune) chain[iter - n_tune, ] = summaries(st)
}
elapsed = proc.time()[["elapsed"]] - t0

# Compare means and a central quantile, with batch-means standard errors for the
# chain (and the effective sample size they imply).
batch_se = function(x, n_batch = 50){
  b = colMeans(matrix(x[seq_len(n_batch * (length(x) %/% n_batch))], ncol = n_batch))
  stats::sd(b) / sqrt(n_batch)
}
rows = lapply(colnames(prior), function(nm){
  p = prior[, nm]
  ch = chain[, nm]
  q = stats::median(p)
  data.frame(summary = nm, prior_mean = mean(p), chain_mean = mean(ch),
             ess = round(stats::var(ch) / batch_se(ch)^2),
             z_mean = (mean(ch) - mean(p)) / sqrt(batch_se(ch)^2 + stats::var(p) / length(p)),
             z_below_median = (mean(ch <= q) - mean(p <= q)) /
               sqrt(batch_se(as.numeric(ch <= q))^2 + 0.25 / length(p)))
})
res = do.call(rbind, rows)
cat(sprintf("Geweke test: boundaries = %s, random_effects = %s, cluster = %s, %d iterations, seed %d (%.0fs)\n",
            boundaries, random_effects, cluster, n_iter, seed, elapsed))
print(format(res, digits = 3), row.names = FALSE)
worst = max(abs(c(res$z_mean, res$z_below_median)), na.rm = TRUE)
cat(sprintf("largest |z| = %.2f over %d statistics\n", worst, 2 * nrow(res)))
