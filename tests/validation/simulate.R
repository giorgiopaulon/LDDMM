# Simulation helpers shared by the validation scripts: draws from the prior of
# every parameter, and response data from the likelihood given the parameters.
# They work on the sampler's internal state (see R/setup.R), so the simulated
# model is exactly the one the sampler targets.

# Truncated inverse-CDF draws.
r_halfcauchy_var = function(upper) tan(stats::runif(1) * atan(upper))  # v ~ C+(0, 1), v <= upper
r_halfcauchy_sd_var = function(scale, upper){                          # sqrt(v) ~ C+(0, scale), v <= upper
  (scale * tan(stats::runif(1) * atan(sqrt(upper) / scale)))^2
}
r_alpha = function(){                                                  # Ga(1, 1) on [0.01, 10]
  stats::qgamma(stats::runif(1, stats::pgamma(0.01, 1), stats::pgamma(10, 1)), 1)
}

# Directed prior of the clustered coefficients given the labels (article eq. 9).
# `pars` holds one list(sigma2, mu0, s20) per clustered parameter (drift, and
# boundary when it shares the labels). Unused labels follow the pseudo-prior:
# a mixture over the used labels at the same knot of N(beta*_{k,j}, tau2),
# drawn jointly for all parameters.
r_cluster_coefs = function(z, Z, tau2, pars){
  J = dim(z)[3]
  beta = lapply(pars, function(p) matrix(NA_real_, J, Z))
  for (k in seq_len(J)){
    zk = z[, , k]
    used = sort(unique(as.vector(zk)))
    for (h in used){
      parents = if (k > 1) unique(z[, , k - 1][zk == h])
      for (nm in names(pars)){
        p = pars[[nm]]
        beta[[nm]][k, h] = if (k == 1) stats::rnorm(1, p$mu0, sqrt(p$s20)) else
          stats::rnorm(1, mean(beta[[nm]][k - 1, parents]), sqrt(p$sigma2 / length(parents)))
      }
    }
    for (u in setdiff(seq_len(Z), used)){
      j = used[sample.int(length(used), 1)]
      for (nm in names(pars)) beta[[nm]][k, u] = stats::rnorm(1, beta[[nm]][k, j], sqrt(tau2))
    }
  }
  beta
}

r_curves = function(n_ind, K, R, Prec){
  L = chol(Prec)
  U = array(0, c(n_ind, K, R))
  for (r in seq_len(R)) U[, , r] = t(backsolve(L, matrix(stats::rnorm(K * n_ind), K)))
  U
}

# Replace every parameter of `st` by a draw from the prior and rebuild the
# parameter caches. `delta_upper` is the upper end of the uniform offset prior.
prior_draw = function(st, delta_upper){
  D = st$design
  hy = st$hypers
  d0 = D$d0
  Z = D$Z
  J = D$J
  if (D$cluster){
    st$alpha_C = r_alpha()
    st$alpha_I = r_alpha()
    for (type in c("C", "I")){
      a = st[[paste0("alpha_", type)]] / Z
      st[[paste0("logQ_", type)]] = t(replicate(Z, lddmm:::log_rdirichlet(rep(a, Z))))
      st[[paste0("logpi0_", type)]] = lddmm:::log_rdirichlet(rep(a, Z))
    }
    z = array(NA_integer_, c(d0, d0, J))
    for (s in seq_len(d0)){
      for (d in seq_len(d0)){
        type = if (s == d) "C" else "I"
        logpi0 = st[[paste0("logpi0_", type)]]
        logQ = st[[paste0("logQ_", type)]]
        z[s, d, 1] = sample.int(Z, 1, prob = exp(logpi0))
        for (k in 2:J) z[s, d, k] = sample.int(Z, 1, prob = exp(logQ[z[s, d, k - 1], ]))
      }
    }
    st$z = z
  }
  st$sigma2_1_mu = r_halfcauchy_sd_var(hy$s_sigma_mu, hy$sigma2_1_max)
  pars = list(mu = list(sigma2 = st$sigma2_1_mu, mu0 = hy$mu_beta0_mu, s20 = hy$sigma2_beta0_mu))
  if (D$flexible){
    st$sigma2_1_b = r_halfcauchy_sd_var(hy$s_sigma_b, hy$sigma2_1_max)
    pars$b = list(sigma2 = st$sigma2_1_b, mu0 = hy$mu_beta0_b, s20 = hy$sigma2_beta0_b)
  }
  coefs = r_cluster_coefs(st$z, Z, hy$pseudo_sd^2, pars)
  st$beta_mu_star = coefs$mu
  if (D$flexible){
    st$beta_b_star = coefs$b
  } else if (D$boundaries == "fixed"){
    st$sigma2_1_b = r_halfcauchy_sd_var(hy$s_sigma_b, hy$sigma2_1_max)
    beta = matrix(NA_real_, J, d0)
    beta[1, ] = stats::rnorm(d0, hy$mu_beta0_b, sqrt(hy$sigma2_beta0_b))
    for (k in 2:J) beta[k, ] = stats::rnorm(d0, beta[k - 1, ], sqrt(st$sigma2_1_b))
    st$beta_b = beta
  } else {
    st$b_const[] = stats::rnorm(length(st$b_const), hy$mu_beta0_b, sqrt(hy$sigma2_beta0_b))
  }
  for (which in c("mu", "b")){
    st[[paste0("sigma2_", which, "_ua")]] = r_halfcauchy_var(hy$sigma2_ua_max)
    st[[paste0("sigma2_", which, "_us")]] = r_halfcauchy_var(hy$sigma2_us_max)
    U = st[[paste0("U_", which)]]
    st[[paste0("U_", which)]] = r_curves(D$n_ind, dim(U)[2], D$R, lddmm:::re_precision(st, which))
  }
  st$delta[] = stats::runif(length(st$delta), 0, delta_upper)
  st$delta_t = st$delta[cbind(D$stim, D$subj)]
  st$lmu_fix = lddmm:::fixed_part(st, "mu")
  st$lb_fix = lddmm:::fixed_part(st, "b")
  st$lmu_re = lddmm:::random_part(st, "mu")
  st$lb_re = lddmm:::random_part(st, "b")
  st
}

# Replace the responses in `st` by a draw from the likelihood (no censoring)
# and rebuild everything that depends on them.
simulate_responses = function(st, delta_upper){
  D = st$design
  mu = exp(st$lmu_fix + st$lmu_re)
  b = exp(st$lb_fix + st$lb_re)
  times = st$delta_t + matrix(lddmm:::rinvgauss_cpp(as.numeric(b / mu), as.numeric(b^2)), D$n)
  st$design$tau = apply(times, 1, min)
  st$design$dec = as.integer(max.col(-times, ties.method = "first") - 1L)
  for (s in seq_len(D$d0)){
    for (i in seq_len(D$n_ind)){
      st$design$delta_max[s, i] = min(delta_upper, st$design$tau[D$rows_si[[s]][[i]]])
    }
  }
  lddmm:::rebuild_likelihood(st)
}

# A sampler state for a balanced design, ready to be overwritten by prior_draw().
toy_state = function(n_ind, n_blocks, d0, reps, boundaries, random_effects, cluster, hypers){
  dat = expand.grid(rep = seq_len(reps), s = seq_len(d0), block = seq_len(n_blocks),
                    subject = seq_len(n_ind))
  dat$d = dat$s
  dat$r_time = stats::rexp(nrow(dat)) + 0.5
  dat$cens = 0
  design = lddmm:::build_design(dat, boundaries, random_effects, cluster, NULL, TRUE)
  lddmm:::init_state(design, hypers)
}
