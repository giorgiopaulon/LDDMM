# Toy example with simulated data: simulate responses from the model's
# likelihood with known parameters, fit lddmm, and compare the posterior with
# the truth.
#
# Run it step by step in RStudio. With the settings below the fit takes a few
# minutes. The file is installed with the package:
# system.file("examples", "toy_simulation.R", package = "lddmm").

library(lddmm)
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw(base_size = 14))
set.seed(2024)

# 1. Design --------------------------------------------------------------------
n_subj = 15            # participants
n_blocks = 8           # training blocks
K = 3                  # input tones = possible decisions
trials_per_tone = 20   # trials per tone and block
J = n_blocks + 1       # spline coefficients (one per knot, plus one)
knot = seq_len(J)

# 2. True parameters -----------------------------------------------------------
# Fixed effects are spline coefficients on the log scale, indexed
# [knot, stimulus s, decision d]. All curves start out flat; from knot 4 on:
# * tone 2 is learned (its correct drift rises) and tone 3 is learned faster;
# * tone 1 is "anti-learned": its correct drift falls while the drift towards
#   response 2 rises above it, so stimulus 1 is answered 2 more and more often;
# * the other incorrect drifts decrease slowly;
# * boundaries are constant: 1.5 for correct responses and for the confusion
#   curve (s = 1, d = 2), 1.8 for the other incorrect responses.
# So the three correct curves share a cluster up to knot 3 and then separate.
# Drift and boundary share cluster labels, so the confusion curve, whose
# boundary differs from the other incorrect curves, is a cluster of its own.
after = pmax(knot - 3, 0) / (J - 3)             # 0 up to knot 3, then rising to 1
drift_wrong = log(0.5 - 0.2 * (knot - 1) / (J - 1))
beta_mu = beta_b = array(NA_real_, c(J, K, K))
for (s in 1:K){
  for (d in 1:K){
    beta_mu[, s, d] = drift_wrong
    beta_b[, s, d] = if (s == d) log(1.5) else log(1.8)
  }
}
beta_mu[, 1, 1] = log(0.85 - 0.3 * after)                         # tone 1, correct: falls
beta_mu[, 2, 2] = log(0.85 + 0.55 * after)                        # tone 2, correct: rises
beta_mu[, 3, 3] = log(0.85 + 0.95 * after)                        # tone 3, correct: rises faster
beta_mu[, 1, 2] = log(exp(drift_wrong) + 0.8 * after)             # tone 1 -> response 2: rises
beta_mu[knot <= 3, 1, 2] = drift_wrong[knot <= 3]                 # (exactly equal before)
beta_b[, 1, 2] = log(1.5)                                         # boundary of correct responses

# Random effects (article eq. 10): every participant has one smooth curve for
# correct responses and one for incorrect responses, for both drift and
# boundary. 'a' is the variance, 's' the smoothness parameter.
truth_var = list(mu = c(a = 0.05, s = 0.01), b = c(a = 0.02, s = 0.005))
P = P_smooth1(J)
draw_curves = function(v){
  L = chol(diag(J) / v[["a"]] + P / v[["s"]])
  replicate(2, t(backsolve(L, matrix(rnorm(J * n_subj), J))), simplify = "array")
}
U_mu = draw_curves(truth_var$mu)   # [subject, knot, curve: 1 correct, 2 incorrect]
U_b = draw_curves(truth_var$b)

# Offsets (non-decision times), one per tone and participant.
delta = matrix(runif(K * n_subj, 0.15, 0.3), K, n_subj)

# 3. Simulate the responses ----------------------------------------------------
# Inverse Gaussian draws (Michael, Schucany and Haas, 1976).
rinvgauss = function(n, mean, shape){
  y = rnorm(n)^2
  x = mean + mean^2 * y / (2 * shape) - mean / (2 * shape) * sqrt(4 * mean * shape * y + mean^2 * y^2)
  ifelse(runif(n) <= mean / (mean + x), x, mean^2 / x)
}

sim = expand.grid(trial = seq_len(trials_per_tone), s = 1:K, block = seq_len(n_blocks),
                  subject = seq_len(n_subj))
n = nrow(sim)
B = B_basis(sim$block, seq_len(n_blocks))
spline_at = function(coef) rowSums(B * coef)  # coef: one row of coefficients per trial

# Each trial is a race: accumulator d finishes after delta + an inverse Gaussian
# time with mean b / mu and shape b^2; the first one decides the response.
finish = matrix(NA_real_, n, K)
for (d in 1:K){
  correct = sim$s == d
  log_mu = spline_at(t(beta_mu[, sim$s, d])) +
    spline_at(U_mu[sim$subject, , 1] * correct + U_mu[sim$subject, , 2] * !correct)
  log_b = spline_at(t(beta_b[, sim$s, d])) +
    spline_at(U_b[sim$subject, , 1] * correct + U_b[sim$subject, , 2] * !correct)
  finish[, d] = delta[cbind(sim$s, sim$subject)] +
    rinvgauss(n, mean = exp(log_b - log_mu), shape = exp(2 * log_b))
}
data_sim = data.frame(subject = sim$subject, block = sim$block, s = sim$s,
                      d = max.col(-finish), r_time = apply(finish, 1, min), cens = 0)
summary(data_sim$r_time)
plot_accuracy(data_sim)
plot_RT(data_sim)

# 4. Fit -----------------------------------------------------------------------
fit = LDDMM(data_sim, hypers = list(s_sigma_mu = 0.1, s_sigma_b = 0.1),
            Niter = 3000, burnin = 1000, thin = 4, save_z = TRUE)

# Posterior predictive checks: decisions and response times simulated from the
# fit, compared with the data.
plot_predictive(data_sim, fit, type = 'response')
plot_predictive(data_sim, fit, type = 'rt')

# 5. True population curves ------------------------------------------------------
# The population curves of the article integrate out the random effects:
# theta(t) = exp{f(t) + var(u(t)) / 2}.
Bk = B_basis(seq_len(n_blocks), seq_len(n_blocks))
re_var = function(v) rowSums((Bk %*% solve(diag(J) / v[["a"]] + P / v[["s"]])) * Bk)
true_pop = function(beta, v) exp(array(Bk %*% matrix(beta, J), c(n_blocks, K, K)) + re_var(v) / 2)
true_mu = true_pop(beta_mu, truth_var$mu)
true_b = true_pop(beta_b, truth_var$b)

# 6. Population curves: posterior vs truth --------------------------------------
curves = function(draws, truth, parameter){
  data.frame(expand.grid(block = seq_len(n_blocks), s = 1:K, d = 1:K),
             mean = as.vector(apply(draws, 1:3, mean)),
             low = as.vector(apply(draws, 1:3, quantile, 0.05)),
             upp = as.vector(apply(draws, 1:3, quantile, 0.95)),
             truth = as.vector(truth), parameter = parameter)
}
pop = rbind(curves(fit$post_mean_mu, true_mu, 'drift'),
            curves(fit$post_mean_b, true_b, 'boundary'))
ggplot(pop, aes(block, group = d)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = factor(d)), alpha = 0.3) +
  geom_line(aes(y = mean, col = factor(d)), linewidth = 1) +
  geom_line(aes(y = truth), linetype = 'dashed') +
  facet_grid(parameter ~ s, scales = 'free_y', labeller = label_both) +
  labs(y = NULL, col = 'decision d', fill = 'decision d',
       caption = 'posterior mean and 90% band; dashed: truth')
# share of true population values inside the 90% bands
pop |> group_by(parameter) |> summarise(coverage = mean(truth >= low & truth <= upp))

# 7. Subject curves: posterior vs truth ------------------------------------------
# Drifts of correct responses for the participants with the largest and the
# smallest true random effect.
true_ind = function(beta, U, i, s) exp(Bk %*% beta[, s, s] + Bk %*% U[i, , 1])
speed = sapply(seq_len(n_subj), function(i) mean(U_mu[i, , 1]))
chosen = c(which.max(speed), which.min(speed))
ind = do.call(rbind, lapply(chosen, function(i) do.call(rbind, lapply(1:K, function(s){
  draws = fit$post_ind_mu[, i, s, s, ]
  data.frame(block = seq_len(n_blocks), s = s, subject = paste('subject', i),
             mean = rowMeans(draws), low = apply(draws, 1, quantile, 0.05),
             upp = apply(draws, 1, quantile, 0.95),
             truth = as.vector(true_ind(beta_mu, U_mu, i, s)))
}))))
ggplot(ind, aes(block, group = subject)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = subject), alpha = 0.3) +
  geom_line(aes(y = mean, col = subject), linewidth = 1) +
  geom_line(aes(y = truth), linetype = 'dashed') +
  facet_wrap(~ s, labeller = label_both) +
  labs(y = 'drift, correct responses', caption = 'dashed: truth')

# 8. Offsets: posterior vs truth --------------------------------------------------
offsets = data.frame(truth = as.vector(delta),
                     mean = as.vector(apply(fit$post_ind_delta, 1:2, mean)),
                     low = as.vector(apply(fit$post_ind_delta, 1:2, quantile, 0.05)),
                     upp = as.vector(apply(fit$post_ind_delta, 1:2, quantile, 0.95)))
ggplot(offsets, aes(truth, mean)) +
  geom_abline(linetype = 'dashed') +
  geom_errorbar(aes(ymin = low, ymax = upp), width = 0, alpha = 0.5) +
  geom_point() +
  labs(x = 'true offset', y = 'posterior mean and 90% interval')

# 9. Local clustering: posterior vs truth ----------------------------------------
# Probability that two curves are equal at each block (both adjacent spline
# coefficients in the same cluster), and the truth. Drift and boundary share the
# cluster labels. Curves are named (s, d) = (stimulus, decision).
compare = list(c(1, 1, 2, 2), c(1, 1, 3, 3), c(2, 2, 3, 3),   # correct responses
               c(1, 2, 1, 3), c(1, 3, 2, 1))                 # incorrect responses
clustering = do.call(rbind, lapply(compare, function(x){
  same = fit$Z[x[1], x[2], , ] == fit$Z[x[3], x[4], , ]          # coefficient x draw
  equal_coef = beta_mu[, x[1], x[2]] == beta_mu[, x[3], x[4]] &
    beta_b[, x[1], x[2]] == beta_b[, x[3], x[4]]
  data.frame(block = seq_len(n_blocks),
             pair = sprintf('(%d,%d) vs (%d,%d)', x[1], x[2], x[3], x[4]),
             posterior = rowMeans(same[-J, ] & same[-1, ]),
             truth = as.numeric(equal_coef[-J] & equal_coef[-1]))
}))
ggplot(clustering, aes(block)) +
  geom_line(aes(y = posterior)) + geom_point(aes(y = posterior)) +
  geom_point(aes(y = truth), shape = 4, size = 3, col = 'red') +
  facet_wrap(~ pair) + ylim(0, 1) +
  labs(y = 'P(curves equal)', caption = 'curves (stimulus, decision); red crosses: truth')

# 10. Random-effect variances: posterior vs truth --------------------------------
variances = data.frame(sigma2_mu_ua = fit$sigma2_mu_ua, sigma2_mu_us = fit$sigma2_mu_us,
                       sigma2_b_ua = fit$sigma2_b_ua, sigma2_b_us = fit$sigma2_b_us) |>
  pivot_longer(everything())
truth_lines = data.frame(name = c('sigma2_mu_ua', 'sigma2_mu_us', 'sigma2_b_ua', 'sigma2_b_us'),
                         value = c(truth_var$mu, truth_var$b))
ggplot(variances, aes(value)) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept = value), data = truth_lines, col = 'red') +
  facet_wrap(~ name, scales = 'free') +
  labs(x = NULL, caption = 'red: truth')

