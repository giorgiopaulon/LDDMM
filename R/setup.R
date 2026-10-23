# Data checks, design quantities, hyperparameters and initial state of the
# sampler.

validate_data = function(data){
  data = as.data.frame(data)
  required = c("subject", "block", "s", "d", "r_time", "cens")
  missing = setdiff(required, colnames(data))
  if (length(missing)){
    stop("`data` is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyNA(data[, c("subject", "block", "s", "r_time", "cens")])){
    stop("`subject`, `block`, `s`, `r_time` and `cens` must not contain NA.", call. = FALSE)
  }
  if (!all(data$cens %in% c(0, 1))){
    stop("`cens` must be 0 (observed) or 1 (censored).", call. = FALSE)
  }
  if (any(data$r_time <= 0)){
    stop("`r_time` must be positive.", call. = FALSE)
  }
  d0 = max(data$s)
  if (any(data$s != round(data$s)) || !setequal(unique(data$s), seq_len(d0))){
    stop("`s` must be coded 1, ..., K with every category present.", call. = FALSE)
  }
  obs = data$cens == 0
  if (anyNA(data$d[obs]) || !all(data$d[obs] %in% seq_len(d0))){
    stop("`d` must use the same codes 1, ..., K as `s` for uncensored trials.", call. = FALSE)
  }
  blocks = sort(unique(data$block))
  if (length(blocks) < 2 || any(blocks != round(blocks)) || any(diff(blocks) != 1)){
    stop("`block` must take at least two consecutive integer values.", call. = FALSE)
  }
  seen = table(factor(data$subject), factor(data$s, levels = seq_len(d0)))
  if (any(seen == 0)){
    stop("Every subject must have at least one trial with every stimulus ",
         "(the offsets delta_s^(i) are otherwise not identified).", call. = FALSE)
  }
  data
}


# Index sets and basis matrices reused at every iteration.
build_design = function(data, boundaries, random_effects, cluster, z_max, informed){
  n = nrow(data)
  subjects = sort(unique(data$subject))
  subj = match(data$subject, subjects)
  n_ind = length(subjects)
  d0 = as.integer(max(data$s))
  knots = min(data$block):max(data$block)
  time = as.integer(data$block - min(data$block) + 1)
  n_blocks = length(knots)
  J = n_blocks + 1
  B = B_basis(data$block, knots)
  P = P_smooth1(J)
  stim = as.integer(data$s)

  if (cluster){
    Z = if (is.null(z_max)) min(d0^2, 2 * d0) else as.integer(z_max)
    if (Z < d0 + 1 || Z > min(d0^2, 32)){
      stop("`z_max` must be between K + 1 and min(K^2, 32), where K is the number of categories.",
           call. = FALSE)
    }
  } else {
    Z = d0^2
  }

  by_cell = function(...) split(seq_len(n), interaction(..., drop = FALSE, lex.order = TRUE))
  # rows_sb[[s]][[t]]: trials with stimulus s in block t
  sb = by_cell(factor(stim, seq_len(d0)), factor(time, seq_len(n_blocks)))
  rows_sb = lapply(seq_len(d0), function(s) unname(sb[(s - 1) * n_blocks + seq_len(n_blocks)]))
  # rows_si[[s]][[i]]: trials of subject i with stimulus s
  si = by_cell(factor(stim, seq_len(d0)), factor(subj, seq_len(n_ind)))
  rows_si = lapply(seq_len(d0), function(s) unname(si[(s - 1) * n_ind + seq_len(n_ind)]))
  # rows_it[[i]][[t]]: trials of subject i in block t
  it = by_cell(factor(subj, seq_len(n_ind)), factor(time, seq_len(n_blocks)))
  rows_it = lapply(seq_len(n_ind), function(i) unname(it[(i - 1) * n_blocks + seq_len(n_blocks)]))
  rows_b = unname(split(seq_len(n), factor(time, seq_len(n_blocks))))

  # Spline coefficient k is supported on blocks k - 1 and k.
  near_knot = function(rows_by_block, k){
    c(if (k > 1) rows_by_block[[k - 1]], if (k <= n_blocks) rows_by_block[[k]])
  }
  rows_sk = lapply(rows_sb, function(rs) lapply(seq_len(J), function(k) near_knot(rs, k)))
  rows_k = lapply(seq_len(J), function(k) near_knot(rows_b, k))
  rows_ik_spline = lapply(rows_it, function(ri) lapply(seq_len(J), function(k) near_knot(ri, k)))
  rows_ik_scalar = lapply(seq_len(n_ind), function(i) list(which(subj == i)))

  # Upper bound of the offset prior: the smallest uncensored response time of
  # subject i to stimulus s (all of them if every such trial is censored).
  delta_max = matrix(NA_real_, d0, n_ind)
  for (s in seq_len(d0)){
    for (i in seq_len(n_ind)){
      rows = rows_si[[s]][[i]]
      unc = rows[data$cens[rows] == 0]
      delta_max[s, i] = min(data$r_time[if (length(unc)) unc else rows])
    }
  }

  list(n = n, d0 = d0, n_ind = n_ind, n_blocks = n_blocks, J = J, Z = Z,
       subjects = subjects, knots = knots, tau = data$r_time,
       dec = ifelse(data$cens == 1, -1L, as.integer(data$d) - 1L),
       stim = stim, subj = subj, time = time, B = B,
       Bgrid = B_basis(knots, knots), P = P,
       P_eigen = pmax(eigen(P, symmetric = TRUE, only.values = TRUE)$values, 0),
       rows_sb = rows_sb, rows_si = rows_si, rows_sk = rows_sk, rows_k = rows_k,
       rows_s = lapply(seq_len(d0), function(s) which(stim == s)),
       rows_ik_spline = rows_ik_spline, rows_ik_scalar = rows_ik_scalar,
       delta_max = delta_max,
       boundaries = boundaries, flexible = boundaries == "flexible",
       spline_boundary = boundaries %in% c("flexible", "fixed"),
       random_effects = random_effects, re_mode = if (random_effects == "shared") 0L else 1L,
       R = if (random_effects == "shared") 1L else 2L,
       cluster = cluster, informed = informed)
}


# Log drift and log boundary matching the first two moments of tau - delta
# (mean b / mu, variance b / mu^3) in every (block, stimulus, decision) cell with
# at least three responses (supplement §S.5.1).
moment_estimates = function(design, delta_t){
  obs = design$dec >= 0
  cells = split(which(obs), list(design$time[obs], design$stim[obs], design$dec[obs] + 1L),
                drop = TRUE)
  est = do.call(rbind, lapply(cells, function(rows){
    x = design$tau[rows] - delta_t[rows]
    if (length(x) < 3 || stats::var(x) <= 0) return(NULL)
    m = mean(x)
    v = stats::var(x)
    data.frame(t = design$time[rows[1]], s = design$stim[rows[1]], d = design$dec[rows[1]] + 1L,
               n = length(x), log_mu = 0.5 * log(m) - 0.5 * log(v),
               log_b = 1.5 * log(m) - 0.5 * log(v))
  }))
  if (is.null(est)) stop("Too few responses to initialize the sampler.", call. = FALSE)
  est
}


# Starting values of the random effects: deviations of each subject's moment
# estimates (pooled over blocks and stimuli) from the population ones, for
# correct and for incorrect responses (or pooled when the random effects are
# shared). Returns an n_ind x 2 x R array with log drift and log boundary
# deviations, centered over subjects.
subject_deviations = function(D, delta_t){
  moments = function(rows){
    x = D$tau[rows] - delta_t[rows]
    if (length(x) < 5 || stats::var(x) <= 0) return(c(NA_real_, NA_real_))
    m = mean(x)
    v = stats::var(x)
    c(0.5 * log(m) - 0.5 * log(v), 1.5 * log(m) - 0.5 * log(v))
  }
  obs = D$dec >= 0
  type = if (D$R == 2) ifelse(D$dec + 1L == D$stim, 1L, 2L) else rep(1L, D$n)
  dev = array(0, c(D$n_ind, 2, D$R))
  for (r in seq_len(D$R)){
    pop = moments(which(obs & type == r))
    for (i in seq_len(D$n_ind)){
      e = moments(which(obs & type == r & D$subj == i)) - pop
      dev[i, , r] = ifelse(is.na(e), 0, e)
    }
    dev_r = matrix(dev[, , r], D$n_ind)
    dev[, , r] = dev_r - rep(colMeans(dev_r), each = D$n_ind)
  }
  dev
}


set_hypers = function(hypers, est){
  if (is.null(hypers)) hypers = list()
  if (is.null(hypers$s_sigma_mu) || is.null(hypers$s_sigma_b)){
    message("`hypers` did not contain s_sigma_mu or s_sigma_b; using 0.1 for both.")
  }
  defaults = list(
    s_sigma_mu = 0.1, s_sigma_b = 0.1,
    sigma2_1_max = 0.2, sigma2_us_max = 0.2, sigma2_ua_max = Inf, pseudo_sd = 0.1,
    mu_beta0_mu = mean(est$log_mu), sigma2_beta0_mu = stats::var(est$log_mu) + 0.25,
    mu_beta0_b = mean(est$log_b), sigma2_beta0_b = stats::var(est$log_b) + 0.25)
  for (nm in names(defaults)){
    if (is.null(hypers[[nm]]) || is.na(hypers[[nm]])) hypers[[nm]] = defaults[[nm]]
  }
  hypers
}


# Coefficients of every (s, d) curve: coef[k, s, d] = beta_star[k, z[s, d, k]].
label_coefs = function(beta_star, z){
  d0 = dim(z)[1]
  J = dim(z)[3]
  vals = beta_star[cbind(rep(seq_len(J), each = d0 * d0), as.vector(z))]
  aperm(array(vals, c(d0, d0, J)), c(3, 1, 2))
}

# Fixed-effect part of a parameter on the log scale for every (trial, accumulator).
fixed_part = function(st, which){
  D = st$design
  if (which == "mu" || D$flexible){
    coef = label_coefs(if (which == "mu") st$beta_mu_star else st$beta_b_star, st$z)
    out = matrix(0, D$n, D$d0)
    for (s in seq_len(D$d0)){
      rows = D$rows_s[[s]]
      out[rows, ] = D$B[rows, , drop = FALSE] %*% coef[, s, ]
    }
    return(out)
  }
  switch(D$boundaries,
         "fixed" = D$B %*% st$beta_b,
         "constant" = st$b_const[D$stim, , drop = FALSE],
         "fixed-constant" = matrix(st$b_const, D$n, D$d0, byrow = TRUE))
}

# Random-effect part on the log scale. Curve 1 acts on the accumulator of the
# correct decision and curve 2 on the others (or curve 1 on all of them when
# the random effects are shared).
random_part = function(st, which){
  D = st$design
  U = if (which == "mu") st$U_mu else st$U_b
  Bu = if (which == "mu" || D$spline_boundary) D$B else matrix(1, D$n, 1)
  per_curve = sapply(seq_len(D$R), function(r){
    rowSums(Bu * matrix(U[D$subj, , r], D$n, ncol(Bu)))
  })
  per_curve = matrix(per_curve, D$n, D$R)
  out = matrix(per_curve[, D$R], D$n, D$d0)
  if (D$R == 2) out[cbind(seq_len(D$n), D$stim)] = per_curve[, 1]
  out
}

rebuild_likelihood = function(st){
  D = st$design
  st$LL = ll_matrix(D$tau, st$delta_t, D$dec, st$lmu_fix + st$lmu_re, st$lb_fix + st$lb_re)
  st
}


init_state = function(design, hypers){
  D = design
  st = list(design = D)
  st$delta = D$delta_max / 2
  st$delta_t = st$delta[cbind(D$stim, D$subj)]

  est = moment_estimates(D, st$delta_t)
  st$hypers = set_hypers(hypers, est)
  hy = st$hypers
  # Pooled starting values for every (s, d) combination.
  start_mu = matrix(hy$mu_beta0_mu, D$d0, D$d0)
  start_b = matrix(hy$mu_beta0_b, D$d0, D$d0)
  pooled = stats::aggregate(cbind(log_mu, log_b) ~ s + d, data = est, FUN = mean)
  start_mu[cbind(pooled$s, pooled$d)] = pooled$log_mu
  start_b[cbind(pooled$s, pooled$d)] = pooled$log_b

  # Labels: supplement §S.5.1 puts every correct-response curve in its own
  # cluster and all incorrect-response curves in one cluster.
  if (D$cluster){
    lab = matrix(D$d0 + 1L, D$d0, D$d0)
    diag(lab) = seq_len(D$d0)
  } else {
    lab = matrix(seq_len(D$d0^2), D$d0, D$d0, byrow = TRUE)
  }
  st$z = array(lab, c(D$d0, D$d0, D$J))
  label_start = function(start, mu0){
    vals = sapply(seq_len(D$Z), function(h) if (any(lab == h)) mean(start[lab == h]) else mu0)
    matrix(vals, D$J, D$Z, byrow = TRUE)
  }
  st$beta_mu_star = label_start(start_mu, hy$mu_beta0_mu)
  st$beta_b_star = if (D$flexible) label_start(start_b, hy$mu_beta0_b) else matrix(0, D$J, D$Z)
  if (D$boundaries == "fixed") st$beta_b = matrix(colMeans(start_b), D$J, D$d0, byrow = TRUE)
  if (D$boundaries == "constant") st$b_const = start_b
  if (D$boundaries == "fixed-constant") st$b_const = colMeans(start_b)

  # Random effects start at constant curves given by the subjects' moment
  # estimates, with their spread as the starting variance.
  dev = subject_deviations(D, st$delta_t)
  K_b = if (D$spline_boundary) D$J else 1L
  st$U_mu = aperm(array(dev[, 1, ], c(D$n_ind, D$R, D$J)), c(1, 3, 2))
  st$U_b = aperm(array(dev[, 2, ], c(D$n_ind, D$R, K_b)), c(1, 3, 2))
  start_var = function(x) min(max(stats::var(as.vector(x)), 0.01), 0.5 * hy$sigma2_ua_max)
  st$sigma2_mu_ua = start_var(dev[, 1, ])
  st$sigma2_b_ua = start_var(dev[, 2, ])
  st$sigma2_mu_us = st$sigma2_b_us = 0.005
  st$sigma2_1_mu = st$sigma2_1_b = 0.005

  # Transition dynamics (article §3.2.1); alpha ~ Ga(1, 1).
  st$alpha_C = st$alpha_I = 1
  st = update_transition_matrices(st)

  st$lmu_fix = fixed_part(st, "mu")
  st$lb_fix = fixed_part(st, "b")
  st$lmu_re = random_part(st, "mu")
  st$lb_re = random_part(st, "b")
  st = rebuild_likelihood(st)

  # Direction of the offset moves (shift_offsets): with m the mean decision
  # time, delta - eps is matched by log mu + eps / (2 m) and log b + 3 eps / (2 m).
  m = mean((D$tau - st$delta_t)[D$dec >= 0])
  st$ridge = c(mu = 1 / (2 * m), b = 3 / (2 * m))

  # Random-walk scales of the updates this model uses, adapted during burn-in.
  mh = function(sd) list(sd = sd, batch = sd * 0, total = sd * 0)
  st$mh = list(
    delta = mh(matrix(0.3, D$d0, D$n_ind)),
    beta1_mu = mh(rep(0.1, D$Z)),
    U_mu = mh(array(0.4, dim(st$U_mu))), U_b = mh(array(0.4, dim(st$U_b))),
    sigma2_1_mu = mh(0.1), sigma2_mu_ua = mh(0.2), sigma2_mu_us = mh(0.1),
    sigma2_b_ua = mh(0.2), scale_mu = mh(0.05), scale_b = mh(0.05),
    ridge_pop = mh(0.005), ridge_subj = mh(rep(0.01, D$n_ind)))
  if (D$flexible) st$mh$beta1_b = mh(rep(0.1, D$Z))
  if (!D$flexible){
    st$mh$boundary = mh(switch(D$boundaries, "fixed" = matrix(0.05, D$J, D$d0),
                               "constant" = matrix(0.05, D$d0, D$d0), rep(0.05, D$d0)))
  }
  if (D$spline_boundary) st$mh[c("sigma2_1_b", "sigma2_b_us")] = list(mh(0.1), mh(0.1))
  if (D$cluster) st$mh[c("alpha_C", "alpha_I")] = list(mh(0.1), mh(0.1))
  st
}
