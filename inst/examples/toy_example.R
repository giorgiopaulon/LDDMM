# Toy example: fit lddmm to the bundled data and look at the fit.
#
# Run it step by step in RStudio, or with Rscript. With the default 1,000
# iterations each fit takes about 1-2 minutes. The file is installed with the
# package: system.file("examples", "toy_example.R", package = "lddmm").

library(lddmm)
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw(base_size = 14))

# 1. Data ---------------------------------------------------------------------
# One row per trial: subject, block, stimulus s, decision d, response time in
# seconds, censoring indicator. 99 subjects, 5 blocks, 2 categories.
data("data")
head(data)
plot_accuracy(data)   # proportion of each decision by block and stimulus
plot_RT(data)         # mean response time by block, stimulus and decision

# 2. Fit ----------------------------------------------------------------------
# For a real analysis use something like Niter = 5000, burnin = 2000, thin = 5.
# save_z = TRUE keeps the cluster labels, needed for the co-clustering plot.
set.seed(123)
fit = LDDMM(data, hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1),
            Niter = 1000, burnin = 500, thin = 2, save_z = TRUE)

# 3. Population curves --------------------------------------------------------
# Posterior means and 90% credible bands of mu_{d,s}(t) and b_{d,s}(t). Panels
# are input stimuli, colours are decisions (the curve with d = s is the
# correct response).
plot_post_pars(data, fit, par = 'drift')
plot_post_pars(data, fit, par = 'boundary')

# 4. Subject-level curves -----------------------------------------------------
# Posterior means of each subject's drift for correct responses.
drift_ind = extract_post_mean(data, fit, par = 'drift') |>
  pivot_longer(starts_with('drift_'), names_to = 'stimulus', values_to = 'drift')
ggplot(drift_ind, aes(block, drift, group = subject)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~ stimulus) +
  labs(y = 'subject drift, correct responses')

# 5. Local clustering ---------------------------------------------------------
# Posterior probability that the correct-response curves (drift and boundary,
# which share cluster labels) of the two stimuli are equal at each block. At
# block t a curve depends on spline coefficients t and t + 1, so the curves are
# equal there when both coefficients share a cluster (article §3.2.1).
Z = fit$Z
same = Z[1, 1, , ] == Z[2, 2, , ]                    # coefficient x draw
equal_at_block = same[-nrow(same), ] & same[-1, ]   # block x draw
coclust = data.frame(block = sort(unique(data$block)), prob = rowMeans(equal_at_block))
ggplot(coclust, aes(block, prob)) +
  geom_line() + geom_point() + ylim(0, 1) +
  labs(y = 'P(curves equal), correct responses s = 1 vs s = 2')

# 6. Posterior predictive checks -----------------------------------------------
# For every saved draw, the fit simulates one response per subject, block and
# stimulus. Compare the simulated responses with the data: proportions of each
# decision, and quantiles of the response times.
plot_predictive(data, fit, type = 'response')
plot_predictive(data, fit, type = 'rt')

# 7. Convergence checks -------------------------------------------------------
# Trace plots of the smoothness and variance parameters, and of one drift and
# boundary: they should fluctuate around a stable level, without trends.
traces = data.frame(draw = seq_along(fit$sigma2_1_mu),
                    sigma2_1_mu = fit$sigma2_1_mu, sigma2_1_b = fit$sigma2_1_b,
                    sigma2_mu_ua = fit$sigma2_mu_ua, sigma2_b_ua = fit$sigma2_b_ua,
                    drift_block1_s1d1 = fit$post_mean_mu[1, 1, 1, ],
                    boundary_block1_s1d1 = fit$post_mean_b[1, 1, 1, ]) |>
  pivot_longer(-draw)
ggplot(traces, aes(draw, value)) +
  geom_line() +
  facet_wrap(~ name, scales = 'free_y')
round(sapply(fit$acceptance, mean), 2)   # Metropolis acceptance rates after burn-in

# 8. Model comparison ---------------------------------------------------------
# WAIC of this model vs constant boundaries over time (smaller is better).
set.seed(123)
fit_const = LDDMM(data, hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1),
                  boundaries = 'constant', Niter = 1000, burnin = 500, thin = 2)
c(flexible = compute_WAIC(fit), constant = compute_WAIC(fit_const))
