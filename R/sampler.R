# The MCMC loop shared by every boundary model.

# Random-walk scales are tuned towards a 44% acceptance rate in batches of 50
# iterations during burn-in only (Roberts and Rosenthal, 2009); after burn-in
# the kernel is fixed, so the chain is a standard Metropolis-within-Gibbs
# sampler.
ADAPT_BATCH = 50

adapt_scales = function(st, n_batch){
  step = min(0.25, 1 / sqrt(n_batch))
  for (name in names(st$mh)){
    rate = st$mh[[name]]$batch / ADAPT_BATCH
    st$mh[[name]]$sd = st$mh[[name]]$sd * exp(step * (2 * (rate > 0.44) - 1))
    st$mh[[name]]$batch = st$mh[[name]]$batch * 0
  }
  st
}

# Debugging aid: every cache must equal its value recomputed from the
# parameters (up to floating-point rounding).
check_caches = function(st, step){
  D = st$design
  diffs = c(
    lmu_fix = max(abs(st$lmu_fix - fixed_part(st, "mu"))),
    lb_fix = max(abs(st$lb_fix - fixed_part(st, "b"))),
    lmu_re = max(abs(st$lmu_re - random_part(st, "mu"))),
    lb_re = max(abs(st$lb_re - random_part(st, "b"))),
    delta_t = max(abs(st$delta_t - st$delta[cbind(D$stim, D$subj)])))
  LL = ll_matrix(D$tau, st$delta_t, D$dec, st$lmu_fix + st$lmu_re, st$lb_fix + st$lb_re)
  finite = is.finite(LL) & is.finite(st$LL)
  diffs["LL"] = max(c(0, abs(LL[finite] - st$LL[finite])))
  if (any(is.finite(LL) != is.finite(st$LL))) diffs["LL"] = Inf
  if (any(diffs > 1e-8)){
    stop(sprintf("cache mismatch after %s: %s", step,
                 paste(names(diffs), signif(diffs, 3), sep = " = ", collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

# One sweep of the sampler: the updates of supplement §S.5.2 (R/update_fixed.R,
# R/update_random.R) interleaved with the joint moves of R/joint_moves.R. Each
# step returns the updated state; `check` runs after every step.
SAMPLER_STEPS = list(
  delta = function(st, pb) update_delta(st, pb),
  coefficients = function(st, pb) update_cluster_coefs(st, pb),
  boundary = function(st, pb) update_boundary(st, pb),
  labels = function(st, pb) update_labels(st),
  transitions = function(st, pb) update_transitions(st, pb),
  translation_mu = function(st, pb) translate(st, "mu"),
  translation_b = function(st, pb) translate(st, "b"),
  smoothness = function(st, pb) update_smoothness(st, pb),
  random_effects_b = function(st, pb) update_random_effects(st, "b", pb),
  re_variances_b = function(st, pb) update_re_variances(st, "b", pb),
  scale_b = function(st, pb) scale_random_effects(st, "b", pb),
  random_effects_mu = function(st, pb) update_random_effects(st, "mu", pb),
  re_variances_mu = function(st, pb) update_re_variances(st, "mu", pb),
  scale_mu = function(st, pb) scale_random_effects(st, "mu", pb),
  offsets = function(st, pb) shift_offsets(st, pb))

# The caches are updated incrementally within a sweep; rebuilding the
# random-effect parts (and the boundary fixed part, which the label update does
# not rebuild in the unclustered boundary models) from the parameters at the
# end of every sweep keeps rounding errors from accumulating over long runs.
sampler_sweep = function(st, post_burnin, check = function(st, step) invisible(TRUE)){
  for (name in names(SAMPLER_STEPS)){
    st = SAMPLER_STEPS[[name]](st, post_burnin)
    check(st, name)
  }
  st$lmu_re = random_part(st, "mu")
  st$lb_re = random_part(st, "b")
  if (!st$design$flexible) st$lb_fix = fixed_part(st, "b")
  st
}

run_sampler = function(st, Niter, burnin, thin, save_z, save_loglik, verbose, debug_checks){
  samp_size = (Niter - burnin) %/% thin
  out = init_output(st$design, samp_size, save_z, save_loglik)
  check = if (debug_checks) check_caches else function(st, step) invisible(TRUE)
  check(st, "initialization")

  if (verbose) pb = utils::txtProgressBar(style = 3)
  it = 0
  for (iter in seq_len(Niter)){
    post_burnin = iter > burnin
    st = sampler_sweep(st, post_burnin, check)
    if (!post_burnin && iter %% ADAPT_BATCH == 0) st = adapt_scales(st, iter / ADAPT_BATCH)
    if (post_burnin && (iter - burnin) %% thin == 0){
      it = it + 1
      out = save_draws(st, out, it)
    }
    if (verbose) utils::setTxtProgressBar(pb, iter / Niter)
  }
  if (verbose) close(pb)

  n_post = Niter - burnin
  out$acceptance = lapply(st$mh, function(m) if (n_post > 0) m$total / n_post else m$total)
  out$settings = list(boundaries = st$design$boundaries,
                      random_effects = st$design$random_effects,
                      cluster = st$design$cluster, z_max = st$design$Z,
                      informed = st$design$informed, hypers = st$hypers,
                      Niter = Niter, burnin = burnin, thin = thin)
  out
}
