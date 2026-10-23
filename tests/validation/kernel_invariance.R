# Invariance test of individual sampler updates.
#
# Draw (parameters, data) independently from the joint model, apply one update
# a few times, and check that the parameters still follow their prior. An
# update that leaves the posterior invariant passes for every data set, so the
# updated parameters are again prior draws. All draws are independent, so the
# tests (two-sample Kolmogorov-Smirnov against a fresh prior sample, and a z
# test of the means) have exact p-values, unlike the Geweke test.
#
# Usage:
#   Rscript tests/validation/kernel_invariance.R <boundaries> <random_effects> <cluster> <update> <draws> [repeats]
# where <update> is the name of one entry of SAMPLER_STEPS (e.g. offsets,
# scale_mu, coefficients, labels) or "sweep" for a full sweep. LDDMM_LIB selects
# the library holding lddmm.

args = commandArgs(trailingOnly = TRUE)
boundaries = args[1]
random_effects = args[2]
cluster = as.logical(args[3])
update = args[4]
n_draws = as.integer(args[5])
repeats = if (length(args) > 5) as.integer(args[6]) else 3

lib = Sys.getenv("LDDMM_LIB", unset = "")
suppressMessages(library(lddmm, lib.loc = if (nzchar(lib)) lib else NULL))
script_dir = dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "simulate.R"))

hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1, sigma2_1_max = 0.2, sigma2_us_max = 0.2,
              sigma2_ua_max = 0.3, mu_beta0_mu = 0, sigma2_beta0_mu = 0.1,
              mu_beta0_b = log(1.5), sigma2_beta0_b = 0.1)
delta_upper = 0.2

summaries = function(st){
  D = st$design
  coef_mu = lddmm:::label_coefs(st$beta_mu_star, st$z)
  out = c(delta_11 = st$delta[1, 1], delta_22 = st$delta[2, 2],
          mu_coef_111 = coef_mu[1, 1, 1], mu_coef_J21 = coef_mu[D$J, 2, 1],
          u_mu_12 = st$U_mu[1, 2, 1], u_b_11 = st$U_b[1, 1, D$R],
          sigma2_mu_ua = st$sigma2_mu_ua, sigma2_b_ua = st$sigma2_b_ua)
  if (D$flexible){
    coef_b = lddmm:::label_coefs(st$beta_b_star, st$z)
    out = c(out, b_coef_111 = coef_b[1, 1, 1], b_coef_212 = coef_b[2, 1, 2])
  } else if (D$boundaries == "fixed"){
    out = c(out, b_coef_11 = st$beta_b[1, 1], b_coef_J2 = st$beta_b[D$J, 2])
  } else {
    out = c(out, b_const_1 = st$b_const[1])
  }
  if (D$cluster) out = c(out, same_label_k2 = as.numeric(st$z[1, 1, 2] == st$z[2, 2, 2]))
  out
}

set.seed(1)
st0 = toy_state(n_ind = 3, n_blocks = 2, d0 = 2, reps = 3, boundaries, random_effects, cluster, hypers)
steps = if (update == "sweep") lddmm:::SAMPLER_STEPS else lddmm:::SAMPLER_STEPS[update]
if (any(vapply(steps, is.null, logical(1)))) stop("unknown update: ", update)

reference = t(replicate(n_draws, summaries(prior_draw(st0, delta_upper))))
after = t(replicate(n_draws, {
  st = simulate_responses(prior_draw(st0, delta_upper), delta_upper)
  for (r in seq_len(repeats)) for (f in steps) st = f(st, FALSE)
  summaries(st)
}))

res = do.call(rbind, lapply(colnames(reference), function(nm){
  x = after[, nm]
  y = reference[, nm]
  binary = all(c(x, y) %in% c(0, 1))
  data.frame(summary = nm, prior_mean = mean(y), after_mean = mean(x),
             z_mean = (mean(x) - mean(y)) / sqrt(stats::var(x) / length(x) + stats::var(y) / length(y)),
             ks_p = if (binary) NA else suppressWarnings(stats::ks.test(x, y)$p.value))
}))
cat(sprintf("Invariance of '%s' (x%d): boundaries = %s, random_effects = %s, cluster = %s, %d draws\n",
            update, repeats, boundaries, random_effects, cluster, n_draws))
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("largest |z| = %.2f; smallest KS p = %.3g\n", max(abs(res$z_mean)), min(res$ks_p, na.rm = TRUE)))
