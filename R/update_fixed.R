# Updates of the offsets, fixed effects, cluster labels, transition dynamics
# and smoothness parameters (supplement §S.5.2 steps 1-6).
#
# Every update keeps the caches in sync: lmu_fix / lb_fix (fixed parts on the
# log scale), lmu_re / lb_re (random parts) and LL (per-trial, per-accumulator
# log-likelihood terms), all n x d0.

record_mh = function(st, name, accepted, post_burnin){
  st$mh[[name]]$batch = st$mh[[name]]$batch + accepted
  if (post_burnin) st$mh[[name]]$total = st$mh[[name]]$total + accepted
  st
}


# (1) Offsets delta_s^{(i)}.
update_delta = function(st, post_burnin){
  D = st$design
  res = update_delta_cpp(D$tau, D$dec, st$lmu_fix + st$lmu_re, st$lb_fix + st$lb_re,
                         st$LL, st$delta, D$rows_si, D$delta_max, st$mh$delta$sd)
  st$delta = res$delta
  st$LL = res$LL
  st$delta_t = st$delta[cbind(D$stim, D$subj)]
  record_mh(st, "delta", res$accepted, post_burnin)
}


# Cache entries touched by coefficient k of the curves (s, d) in `members`
# (linear indices into a d0 x d0 matrix): trials of stimulus s at blocks k - 1
# and k, accumulator d, with spline weight B[t, k].
coef_entries = function(D, members, k){
  s = (members - 1) %% D$d0 + 1
  d = (members - 1) %/% D$d0 + 1
  rows = lapply(s, function(x) D$rows_sk[[x]][[k]])
  len = lengths(rows)
  rows = unlist(rows, use.names = FALSE)
  cols = rep(d, len)
  list(rows = rows, cols = cols, lin = rows + (cols - 1) * D$n, w = D$B[cbind(rows, k)])
}

# log p(beta_{k+1, h'} | beta_{k, .}) for the children h' of label h at knot k,
# with beta_{k, h} set to `value` (article eq. 9, normalized reading: a label
# with parent set P is N(mean(beta_{k, P}), sigma2 / |P|)).
log_children = function(beta, zk, znext, k, h, children, value, sigma2){
  beta_k = beta[k, ]
  beta_k[h] = value
  out = 0
  for (h2 in children){
    parents = unique(zk[znext == h2])
    out = out + log_dnorm(beta[k + 1, h2], mean(beta_k[parents]), sigma2 / length(parents))
  }
  out
}

# Pseudo-prior of unused labels (they integrate out, so it does not affect the
# posterior; Carlin and Chib, 1995): a mixture over the labels in use at the
# same knot of N(beta*_{k,j}, tau2), jointly for drift and boundary when they
# share labels. Returns the summed log density of the unused labels' values.
log_pseudo = function(bm_k, bb_k, used, unused, tau2, flex){
  out = 0
  for (u in unused){
    v = log_dnorm(bm_k[u], bm_k[used], tau2)
    if (flex) v = v + log_dnorm(bb_k[u], bb_k[used], tau2)
    out = out + logsumexp(v) - log(length(used))
  }
  out
}

# (2) Cluster-specific spline coefficients beta*_{k,h} (drift; drift and
# boundary jointly when boundaries = "flexible", since they share labels).
#   * unused label: exact draw from the pseudo-prior;
#   * k >= 2: independence proposal from the conditional prior
#     N(mean of the parent coefficients, sigma2_1 / #parents);
#   * k == 1: random walk under the N(mu0, sigma2_0) prior.
# Moves of a used label also change the pseudo-prior density of the unused
# labels at the same knot, which enters the acceptance ratio.
update_cluster_coefs = function(st, post_burnin){
  D = st$design
  hy = st$hypers
  flex = D$flexible
  bm = st$beta_mu_star
  bb = st$beta_b_star
  tau2 = hy$pseudo_sd^2
  acc_mu = acc_b = numeric(D$Z)
  for (k in seq_len(D$J)){
    zk = st$z[, , k]
    zprev = if (k > 1) st$z[, , k - 1]
    znext = if (k < D$J) st$z[, , k + 1]
    used = sort(unique(as.vector(zk)))
    unused = setdiff(seq_len(D$Z), used)
    for (h in seq_len(D$Z)){
      members = which(zk == h)
      if (!length(members)){
        j = used[sample.int(length(used), 1)]
        bm[k, h] = stats::rnorm(1, bm[k, j], sqrt(tau2))
        if (flex) bb[k, h] = stats::rnorm(1, bb[k, j], sqrt(tau2))
        next
      }
      prop_b = bb[k, h]
      if (k > 1){
        parents = unique(zprev[members])
        prop_mu = stats::rnorm(1, mean(bm[k - 1, parents]), sqrt(st$sigma2_1_mu / length(parents)))
        if (flex) prop_b = stats::rnorm(1, mean(bb[k - 1, parents]), sqrt(st$sigma2_1_b / length(parents)))
        log_ratio = 0
      } else {
        prop_mu = bm[1, h] + st$mh$beta1_mu$sd[h] * stats::rnorm(1)
        log_ratio = log_dnorm(prop_mu, hy$mu_beta0_mu, hy$sigma2_beta0_mu) -
          log_dnorm(bm[1, h], hy$mu_beta0_mu, hy$sigma2_beta0_mu)
        if (flex){
          prop_b = bb[1, h] + st$mh$beta1_b$sd[h] * stats::rnorm(1)
          log_ratio = log_ratio + log_dnorm(prop_b, hy$mu_beta0_b, hy$sigma2_beta0_b) -
            log_dnorm(bb[1, h], hy$mu_beta0_b, hy$sigma2_beta0_b)
        }
      }
      if (k < D$J){
        children = unique(znext[members])
        log_ratio = log_ratio +
          log_children(bm, zk, znext, k, h, children, prop_mu, st$sigma2_1_mu) -
          log_children(bm, zk, znext, k, h, children, bm[k, h], st$sigma2_1_mu)
        if (flex){
          log_ratio = log_ratio +
            log_children(bb, zk, znext, k, h, children, prop_b, st$sigma2_1_b) -
            log_children(bb, zk, znext, k, h, children, bb[k, h], st$sigma2_1_b)
        }
      }
      if (length(unused)){
        bm_k = bm[k, ]
        bb_k = bb[k, ]
        bm_k[h] = prop_mu
        bb_k[h] = prop_b
        log_ratio = log_ratio + log_pseudo(bm_k, bb_k, used, unused, tau2, flex) -
          log_pseudo(bm[k, ], bb[k, ], used, unused, tau2, flex)
      }
      e = coef_entries(D, members, k)
      shift_mu = e$w * (prop_mu - bm[k, h])
      shift_b = e$w * (prop_b - bb[k, h])
      ll_new = ll_entries(D$tau, st$delta_t, D$dec, e$rows, e$cols,
                          st$lmu_fix[e$lin] + shift_mu + st$lmu_re[e$lin],
                          st$lb_fix[e$lin] + shift_b + st$lb_re[e$lin])
      log_ratio = log_ratio + sum(ll_new) - sum(st$LL[e$lin])
      if (log(stats::runif(1)) < log_ratio){
        bm[k, h] = prop_mu
        st$lmu_fix[e$lin] = st$lmu_fix[e$lin] + shift_mu
        if (flex){
          bb[k, h] = prop_b
          st$lb_fix[e$lin] = st$lb_fix[e$lin] + shift_b
        }
        st$LL[e$lin] = ll_new
        if (k == 1){
          acc_mu[h] = 1
          acc_b[h] = 1
        }
      }
    }
  }
  st$beta_mu_star = bm
  st$beta_b_star = bb
  st = record_mh(st, "beta1_mu", acc_mu, post_burnin)
  if (flex) st = record_mh(st, "beta1_b", acc_b, post_burnin)
  st
}


# (2') Boundaries that are not clustered: random-walk updates of
#   * boundaries = "fixed":          spline b_d(t), first-difference prior;
#   * boundaries = "constant":       scalar b_{d,s};
#   * boundaries = "fixed-constant": scalar b_d.
# Starting coefficients and the scalars have a N(mu0_b, sigma2_0_b) prior.
update_boundary = function(st, post_burnin){
  D = st$design
  hy = st$hypers
  if (D$flexible) return(st)
  sd = st$mh$boundary$sd
  acc = sd * 0
  # Metropolis step for one boundary coefficient acting on trials `rows` of
  # accumulator `col` with spline weights `w`.
  try_move = function(st, rows, col, w, cur, prop, log_prior_diff){
    lin = rows + (col - 1) * D$n
    shift = w * (prop - cur)
    ll_new = ll_entries(D$tau, st$delta_t, D$dec, rows, rep(col, length(rows)),
                        st$lmu_fix[lin] + st$lmu_re[lin], st$lb_fix[lin] + shift + st$lb_re[lin])
    accepted = log(stats::runif(1)) < log_prior_diff + sum(ll_new) - sum(st$LL[lin])
    if (accepted){
      st$lb_fix[lin] = st$lb_fix[lin] + shift
      st$LL[lin] = ll_new
    }
    list(st = st, accepted = accepted)
  }

  if (D$boundaries == "fixed"){
    for (k in seq_len(D$J)){
      rows = D$rows_k[[k]]
      w = D$B[cbind(rows, k)]
      for (d in seq_len(D$d0)){
        cur = st$beta_b[k, d]
        prop = cur + sd[k, d] * stats::rnorm(1)
        lp = function(x){
          out = if (k == 1) log_dnorm(x, hy$mu_beta0_b, hy$sigma2_beta0_b) else
            log_dnorm(x, st$beta_b[k - 1, d], st$sigma2_1_b)
          if (k < D$J) out = out + log_dnorm(st$beta_b[k + 1, d], x, st$sigma2_1_b)
          out
        }
        move = try_move(st, rows, d, w, cur, prop, lp(prop) - lp(cur))
        st = move$st
        if (move$accepted){
          st$beta_b[k, d] = prop
          acc[k, d] = 1
        }
      }
    }
  } else {
    cells = if (D$boundaries == "constant") as.matrix(expand.grid(s = seq_len(D$d0), d = seq_len(D$d0))) else
      cbind(s = NA, d = seq_len(D$d0))
    for (j in seq_len(nrow(cells))){
      s = cells[j, "s"]
      d = cells[j, "d"]
      idx = if (is.na(s)) d else cbind(s, d)
      rows = if (is.na(s)) seq_len(D$n) else D$rows_s[[s]]
      cur = st$b_const[idx]
      prop = cur + sd[idx] * stats::rnorm(1)
      move = try_move(st, rows, d, 1, cur, prop,
                      log_dnorm(prop, hy$mu_beta0_b, hy$sigma2_beta0_b) -
                        log_dnorm(cur, hy$mu_beta0_b, hy$sigma2_beta0_b))
      st = move$st
      if (move$accepted){
        st$b_const[idx] = prop
        acc[idx] = 1
      }
    }
  }
  record_mh(st, "boundary", acc, post_burnin)
}


# (3) Cluster labels, one input stimulus at a time (exact locally informed
# Hamming ball sampler in src/labels.cpp), followed by a rebuild of the caches
# that depend on the labels.
update_labels = function(st){
  D = st$design
  if (!D$cluster) return(st)
  hy = st$hypers
  d0 = D$d0
  z_all = matrix(aperm(st$z, c(2, 1, 3)), d0 * d0, D$J)  # row (s - 1) * d0 + d
  storage.mode(z_all) = "integer"
  for (s in seq_len(d0)){
    z_all[(s - 1) * d0 + seq_len(d0), ] = hb_update_stimulus(
      s, z_all, D$Z, st$beta_mu_star, st$beta_b_star, D$flexible, D$rows_sb[[s]],
      D$tau, st$delta_t, D$dec, st$lmu_re, st$lb_fix, st$lb_re,
      st$logQ_C, st$logQ_I, st$logpi0_C, st$logpi0_I,
      st$sigma2_1_mu, st$sigma2_1_b, hy$mu_beta0_mu, hy$sigma2_beta0_mu,
      hy$mu_beta0_b, hy$sigma2_beta0_b, hy$pseudo_sd^2, D$informed)
  }
  st$z = aperm(array(z_all, c(d0, d0, D$J)), c(2, 1, 3))
  st$lmu_fix = fixed_part(st, "mu")
  if (D$flexible) st$lb_fix = fixed_part(st, "b")
  rebuild_likelihood(st)
}


# Label paths of the correct (d = s) and incorrect (d != s) chains.
label_paths = function(st){
  d0 = st$design$d0
  correct = diag(d0) == 1
  paths = matrix(st$z, d0 * d0, st$design$J)  # row (d - 1) * d0 + s
  list(C = paths[as.vector(correct), , drop = FALSE], I = paths[!as.vector(correct), , drop = FALSE])
}

# (5) Transition matrices and initial distributions: rows ~ Dirichlet(alpha / z_max + counts).
update_transition_matrices = function(st){
  Z = st$design$Z
  paths = label_paths(st)
  for (type in c("C", "I")){
    p = paths[[type]]
    alpha = st[[paste0("alpha_", type)]]
    counts = matrix(tabulate((p[, -ncol(p)] - 1) * Z + p[, -1], Z * Z), Z, Z, byrow = TRUE)
    st[[paste0("logQ_", type)]] = t(apply(counts, 1, function(x) log_rdirichlet(alpha / Z + x)))
    st[[paste0("logpi0_", type)]] = log_rdirichlet(alpha / Z + tabulate(p[, 1], Z))
  }
  st
}

# Concentrations alpha^(C), alpha^(I) ~ Ga(1, 1), restricted to [0.01, 10].
update_alpha = function(st, post_burnin){
  Z = st$design$Z
  for (type in c("C", "I")){
    logQ = st[[paste0("logQ_", type)]]
    logpi0 = st[[paste0("logpi0_", type)]]
    sum_log = sum(logQ) + sum(logpi0)
    target = function(a){
      stats::dgamma(a, 1, 1, log = TRUE) + (Z + 1) * (lgamma(a) - Z * lgamma(a / Z)) +
        (a / Z - 1) * sum_log
    }
    name = paste0("alpha_", type)
    res = mh_positive(st[[name]], target, st$mh[[name]]$sd, upper = 10, lower = 0.01)
    st[[name]] = res$value
    st = record_mh(st, name, res$accepted, post_burnin)
  }
  st
}

update_transitions = function(st, post_burnin){
  if (!st$design$cluster) return(st)
  st = update_transition_matrices(st)
  update_alpha(st, post_burnin)
}


# Sufficient statistics of the smoothness parameter of the clustered
# coefficients: sum over labels in use at knots k >= 2 of
# m (beta_{k,h} - mean of its m parents)^2, and the number of such labels.
smoothness_stats = function(beta, z){
  J = dim(z)[3]
  ss = 0
  count = 0
  for (k in 2:J){
    zk = z[, , k]
    zp = z[, , k - 1]
    for (h in unique(as.vector(zk))){
      parents = unique(zp[zk == h])
      ss = ss + length(parents) * (beta[k, h] - mean(beta[k - 1, parents]))^2
      count = count + 1
    }
  }
  c(ss = ss, count = count)
}

# (6) Smoothness parameters sigma2_1 of the fixed-effect curves: sqrt(sigma2_1)
# ~ half-Cauchy(0, s_sigma), truncated at sigma2_1_max.
update_smoothness = function(st, post_burnin){
  D = st$design
  hy = st$hypers
  draw = function(st, name, stats, scale){
    target = function(v) -0.5 * stats[["count"]] * log(v) - 0.5 * stats[["ss"]] / v +
      log_prior_sd_halfcauchy(v, scale)
    res = mh_positive(st[[name]], target, st$mh[[name]]$sd, upper = hy$sigma2_1_max)
    st[[name]] = res$value
    record_mh(st, name, res$accepted, post_burnin)
  }
  st = draw(st, "sigma2_1_mu", smoothness_stats(st$beta_mu_star, st$z), hy$s_sigma_mu)
  if (D$flexible){
    st = draw(st, "sigma2_1_b", smoothness_stats(st$beta_b_star, st$z), hy$s_sigma_b)
  } else if (D$boundaries == "fixed"){
    st = draw(st, "sigma2_1_b", c(ss = sum(diff(st$beta_b)^2), count = (D$J - 1) * D$d0),
              hy$s_sigma_b)
  }
  st
}
