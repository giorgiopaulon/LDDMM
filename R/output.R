# Storage of posterior draws, posterior predictive draws, and the functions that
# summarize a fit.

init_output = function(D, samp_size, save_z, save_loglik){
  G = nrow(D$Bgrid)
  d0 = D$d0
  out = list(
    post_mean_delta = array(NA_real_, c(samp_size, d0)),
    post_mean_mu = array(NA_real_, c(G, d0, d0, samp_size)),
    post_mean_b = array(NA_real_, c(G, d0, d0, samp_size)),
    post_ind_delta = array(NA_real_, c(d0, D$n_ind, samp_size)),
    post_ind_mu = array(NA_real_, c(G, D$n_ind, d0, d0, samp_size)),
    post_ind_b = array(NA_real_, c(G, D$n_ind, d0, d0, samp_size)),
    sigma2_mu_us = rep(NA_real_, samp_size), sigma2_mu_ua = rep(NA_real_, samp_size),
    sigma2_b_us = rep(NA_real_, samp_size), sigma2_b_ua = rep(NA_real_, samp_size),
    sigma2_1_mu = rep(NA_real_, samp_size), sigma2_1_b = rep(NA_real_, samp_size),
    pred_ans = array(NA_integer_, c(D$n_blocks, d0, samp_size)),
    pred_time = array(NA_real_, c(D$n_blocks, d0, samp_size)),
    pred_ans_ind = array(NA_integer_, c(D$n_ind, D$n_blocks, d0, samp_size)),
    pred_time_ind = array(NA_real_, c(D$n_ind, D$n_blocks, d0, samp_size)))
  if (save_loglik) out$loglik = array(NA_real_, c(samp_size, D$n))
  if (save_z) out$Z = array(NA_integer_, c(d0, d0, D$J, samp_size))
  out
}

# Curves of one parameter on the block grid, on the log scale:
#   fixed[g, s, d] (population fixed effects),
#   re[[r]][g, i]  (random-effect curve r of subject i),
#   var_re[g]      (prior variance of a random-effect curve, for the
#                   population-level correction exp(f + var / 2); article §3.2.2).
grid_curves = function(st, which){
  D = st$design
  G = nrow(D$Bgrid)
  d0 = D$d0
  U = st[[paste0("U_", which)]]
  Prec = re_precision(st, which)
  spline = which == "mu" || D$spline_boundary
  if (which == "mu" || D$flexible){
    coef = label_coefs(if (which == "mu") st$beta_mu_star else st$beta_b_star, st$z)
    fixed = array(D$Bgrid %*% matrix(coef, D$J), c(G, d0, d0))
  } else if (D$boundaries == "fixed"){
    fixed = aperm(array(D$Bgrid %*% st$beta_b, c(G, d0, d0)), c(1, 3, 2))
  } else if (D$boundaries == "constant"){
    fixed = aperm(array(st$b_const, c(d0, d0, G)), c(3, 1, 2))
  } else {
    fixed = aperm(array(st$b_const, c(d0, d0, G)), c(3, 2, 1))
  }
  Bu = if (spline) D$Bgrid else matrix(1, G, 1)
  re = lapply(seq_len(D$R), function(r) Bu %*% t(matrix(U[, , r], D$n_ind)))
  var_re = rowSums((Bu %*% solve(Prec)) * Bu)
  list(fixed = fixed, re = re, var_re = var_re)
}

# Accumulator d of stimulus s uses random-effect curve 1 (correct, d == s) or
# 2 (incorrect); with shared random effects there is a single curve.
re_index = function(D, s, d) if (D$R == 1 || s == d) 1 else 2

# Race of the d0 accumulators: offset plus inverse Gaussian first-passage
# times; returns the winning accumulator and its response time.
simulate_race = function(log_mu, log_b, delta){
  mu = exp(log_mu)
  b = exp(log_b)
  times = matrix(delta + rinvgauss_cpp(as.numeric(b / mu), as.numeric(b^2)), ncol = ncol(log_mu))
  list(ans = max.col(-times, ties.method = "first"), time = apply(times, 1, min))
}

save_draws = function(st, out, it){
  D = st$design
  d0 = D$d0
  G = nrow(D$Bgrid)
  mu = grid_curves(st, "mu")
  b = grid_curves(st, "b")

  out$post_mean_mu[, , , it] = exp(mu$fixed + 0.5 * mu$var_re)
  out$post_mean_b[, , , it] = exp(b$fixed + 0.5 * b$var_re)
  for (s in seq_len(d0)){
    for (d in seq_len(d0)){
      out$post_ind_mu[, , s, d, it] = exp(mu$fixed[, s, d] + mu$re[[re_index(D, s, d)]])
      out$post_ind_b[, , s, d, it] = exp(b$fixed[, s, d] + b$re[[re_index(D, s, d)]])
    }
  }
  out$post_ind_delta[, , it] = st$delta
  out$post_mean_delta[it, ] = rowMeans(st$delta)
  out$sigma2_mu_us[it] = st$sigma2_mu_us
  out$sigma2_mu_ua[it] = st$sigma2_mu_ua
  out$sigma2_b_ua[it] = st$sigma2_b_ua
  out$sigma2_1_mu[it] = st$sigma2_1_mu
  # Constant boundaries have neither smooth curves nor a smoothness parameter.
  if (D$spline_boundary){
    out$sigma2_b_us[it] = st$sigma2_b_us
    out$sigma2_1_b[it] = st$sigma2_1_b
  }

  # Posterior predictive draws at every block: population level (population
  # curves, average offset) and subject level.
  pop = simulate_race(matrix(log(out$post_mean_mu[, , , it]), G * d0),
                      matrix(log(out$post_mean_b[, , , it]), G * d0),
                      rep(rowMeans(st$delta), each = G))
  out$pred_ans[, , it] = pop$ans
  out$pred_time[, , it] = pop$time
  # rows ordered (subject, block, stimulus), columns are the accumulators
  dims = c(G, D$n_ind, d0, d0)
  ind_mu = aperm(array(log(out$post_ind_mu[, , , , it]), dims), c(2, 1, 3, 4))
  ind_b = aperm(array(log(out$post_ind_b[, , , , it]), dims), c(2, 1, 3, 4))
  ind = simulate_race(matrix(ind_mu, D$n_ind * G * d0), matrix(ind_b, D$n_ind * G * d0),
                      as.vector(apply(st$delta, 1, rep, times = G)))
  out$pred_ans_ind[, , , it] = ind$ans
  out$pred_time_ind[, , , it] = ind$time

  if (!is.null(out$loglik)) out$loglik[it, ] = rowSums(st$LL)
  if (!is.null(out$Z)) out$Z[, , , it] = st$z
  out
}


#' Parameter posterior means
#'
#' Posterior means of the subject-level drift or boundary parameters of the
#' correct responses (d = s) from an [LDDMM()] fit.
#'
#' @param data the data frame used to fit the model (see [LDDMM()]).
#' @param fit fit from the [LDDMM()] function
#' @param par parameter to output ('drift', or 'boundary')
#' @return Data frame with columns `block`, `subject` and one column
#'   `<par>_<s>_<s>` per stimulus with the posterior means.
#' @export
extract_post_mean = function(data, fit, par = c('drift', 'boundary')){
  block = d = s = subject = value = NULL
  par = match.arg(par)
  draws = if (par == 'drift') fit$post_ind_mu else fit$post_ind_b
  post_mean = apply(draws, 1:4, mean)
  dplyr::tibble(reshape2::melt(post_mean, varnames = c('block', 'subject', 's', 'd'))) |>
    dplyr::filter(s == d) |>
    dplyr::mutate(block = (min(data$block):max(data$block))[block],
                  subject = sort(unique(data$subject))[subject]) |>
    tidyr::pivot_wider(names_from = c(s, d), values_from = value,
                       names_prefix = paste0(par, '_'))
}


#' Parameter posterior draws
#'
#' Posterior draws of the subject-level drift or boundary parameters of the
#' correct responses (d = s) from an [LDDMM()] fit.
#'
#' @inheritParams extract_post_mean
#' @return Data frame with columns `block`, `subject`, `draw` and one column
#'   `<par>_<s>_<s>` per stimulus with the posterior draws.
#' @export
extract_post_draws = function(data, fit, par = c('drift', 'boundary')){
  block = d = s = subject = value = NULL
  par = match.arg(par)
  draws = if (par == 'drift') fit$post_ind_mu else fit$post_ind_b
  dplyr::tibble(reshape2::melt(draws, varnames = c('block', 'subject', 's', 'd', 'draw'))) |>
    dplyr::filter(s == d) |>
    dplyr::mutate(block = (min(data$block):max(data$block))[block],
                  subject = sort(unique(data$subject))[subject]) |>
    tidyr::pivot_wider(names_from = c(s, d), values_from = value,
                       names_prefix = paste0(par, '_'))
}


#' Plot posterior estimates
#'
#' Posterior means and 90% credible intervals of the population-level drift or
#' boundary curves from an [LDDMM()] fit. Each panel is an input stimulus `s`;
#' colours are the decisions `d`.
#'
#' @inheritParams extract_post_mean
#' @return A `ggplot` object.
#' @export
plot_post_pars = function(data, fit, par = c('drift', 'boundary')){
  block = d = low = upp = value = NULL
  par = match.arg(par)
  draws = if (par == 'drift') fit$post_mean_mu else fit$post_mean_b
  tex_string = if (par == 'drift') "$\\mu_{d,s}(t)$" else "$b_{d,s}(t)$"
  blocks = min(data$block):max(data$block)

  summ = reshape2::melt(apply(draws, 1:3, mean), varnames = c('block', 's', 'd'))
  summ$low = as.vector(apply(draws, 1:3, stats::quantile, probs = 0.05))
  summ$upp = as.vector(apply(draws, 1:3, stats::quantile, probs = 0.95))
  summ$block = blocks[summ$block]

  ggplot(summ) +
    geom_line(aes(x = block, y = value, col = factor(d)), linewidth = 1.5) +
    geom_ribbon(aes(x = block, ymin = low, ymax = upp, fill = factor(d)), alpha = 0.4) +
    facet_wrap( ~ s, nrow = 2) +
    labs(x = 'block', y = TeX(tex_string)) +
    scale_color_brewer(name = "", palette = "Set1") +
    scale_fill_brewer(name = "", palette = "Set1") +
    scale_x_continuous(breaks = blocks) +
    theme(legend.position = "none")
}


#' Posterior predictive checks
#'
#' Compares the data with responses simulated from the fitted model. For every
#' saved draw, [LDDMM()] simulates one response of every subject to every
#' stimulus in every block (`pred_ans_ind` and `pred_time_ind`).
#' * `type = "response"`: proportion of each decision, by stimulus and block.
#'   Lines are the posterior predictive mean of the proportion; bands are the
#'   90% interval of the proportion among the simulated responses of one draw
#'   (one per participant), so they also contain the sampling variability of
#'   that many responses. Points are the observed proportions.
#' * `type = "rt"`: 10%, 50% and 90% quantiles of the response times, by
#'   stimulus and block. Lines are the posterior predictive quantiles (pooled
#'   over subjects and draws); points are the observed quantiles.
#'
#' Censored trials are left out of the observed data.
#'
#' @inheritParams extract_post_mean
#' @param type which check to plot: `"response"` or `"rt"`.
#' @return A `ggplot` object.
#' @export
plot_predictive = function(data, fit, type = c('response', 'rt')){
  block = s = d = r_time = n = freq = value = low = upp = quantile = NULL
  type = match.arg(type)
  blocks = min(data$block):max(data$block)
  n_ind = dim(fit$pred_ans_ind)[1]
  K = dim(fit$pred_ans_ind)[3]
  observed = data[data$cens == 0, ]
  # predictions of block t and stimulus s: subjects x draws
  pred = function(x, t, s) matrix(x[, t, s, ], n_ind)

  if (type == 'response'){
    grid = expand.grid(t = seq_along(blocks), s = seq_len(K), d = seq_len(K))
    summ = t(mapply(function(t, s, d){
      p = colMeans(pred(fit$pred_ans_ind, t, s) == d)
      c(mean = mean(p), low = stats::quantile(p, 0.05, names = FALSE),
        upp = stats::quantile(p, 0.95, names = FALSE))
    }, grid$t, grid$s, grid$d))
    predicted = data.frame(block = blocks[grid$t], s = grid$s, d = grid$d, summ)
    # observed proportions, with zeros for decisions never given
    obs = merge(expand.grid(block = blocks, s = seq_len(K), d = seq_len(K)),
                dplyr::count(observed, block, s, d), all.x = TRUE)
    obs$n[is.na(obs$n)] = 0
    obs = obs |>
      dplyr::group_by(block, s) |>
      dplyr::mutate(freq = n / sum(n)) |>
      dplyr::ungroup()
    return(
      ggplot(predicted, aes(x = block)) +
        geom_ribbon(aes(ymin = low, ymax = upp, fill = factor(d)), alpha = 0.25) +
        geom_line(aes(y = mean, col = factor(d))) +
        geom_point(aes(y = freq, col = factor(d)), data = obs) +
        facet_wrap(~ s, labeller = label_both) +
        scale_color_brewer(name = 'decision', palette = 'Set1') +
        scale_fill_brewer(name = 'decision', palette = 'Set1') +
        scale_x_continuous(breaks = blocks) +
        labs(y = 'proportion of responses',
             caption = paste('lines: posterior predictive mean; bands: 90% interval with one',
                             'simulated response per participant; points: data'))
    )
  }

  probs = c(0.1, 0.5, 0.9)
  grid = expand.grid(t = seq_along(blocks), s = seq_len(K))
  q = t(mapply(function(t, s) stats::quantile(pred(fit$pred_time_ind, t, s), probs, names = FALSE),
               grid$t, grid$s))
  predicted = data.frame(block = blocks[rep(grid$t, length(probs))], s = rep(grid$s, length(probs)),
                         quantile = rep(probs, each = nrow(grid)), value = as.vector(q))
  obs = observed |>
    dplyr::group_by(block, s) |>
    dplyr::reframe(quantile = probs, value = stats::quantile(r_time, probs, names = FALSE))
  ggplot(predicted, aes(x = block, y = value, col = factor(quantile))) +
    geom_line() +
    geom_point(data = obs) +
    facet_wrap(~ s, labeller = label_both) +
    scale_color_brewer(name = 'quantile', palette = 'Dark2') +
    scale_x_continuous(breaks = blocks) +
    labs(y = 'response time', caption = 'lines: posterior predictive quantiles; points: data')
}


#' Calculate WAIC
#'
#' Watanabe-Akaike information criterion (Gelman, Hwang and Vehtari, 2014),
#' which estimates the expected out-of-sample prediction error using a
#' bias-corrected adjustment of the within-sample error.
#'
#' @param model_fit results of a model fit from the [LDDMM()] function, run with
#'   `save_loglik = TRUE`
#' @return A scalar indicating the WAIC (smaller WAIC denotes better fit)
#' @export
compute_WAIC = function(model_fit) {
  if (is.null(model_fit$loglik)) stop("The fit has no `loglik`; refit with save_loglik = TRUE.")
  loglik = model_fit$loglik
  # log pointwise predictive density, computed stably on the log scale
  lppd = sum(apply(loglik, 2, logsumexp)) - ncol(loglik) * log(nrow(loglik))
  p_waic = sum(apply(loglik, 2, stats::var))
  -2 * (lppd - p_waic)
}
