
require(Rcpp)
require(RcppArmadillo)
require(mvtnorm)
require(gtools)
require(ggplot2)
require(latex2exp)
require(RColorBrewer)
require(LaplacesDemon)
require(plyr)

`%notin%` = Negate(`%in%`)


#' Descriptive plots 
#'
#' Plot the accuracy of the raw data.
#' 
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @return Individual and population level raw accuracies
plot_accuracy = function(data){
  
  block = d = s = cens = subject = tally = n = freq = NULL
  
  data_aggr = data |>
    dplyr::filter(cens == 0) |>
    dplyr::group_by(block, s, d) |> 
    tally() |> 
    dplyr::mutate(freq = n / sum(n), 
                  s = factor(s), 
                  d = factor(d)) |> 
    dplyr::ungroup() |> 
    dplyr::select(block, s, d, freq)
  
  data |>
    dplyr::filter(cens == 0) |>
    dplyr::group_by(subject, block, s, d) |> 
    tally() |> 
    dplyr::mutate(freq = n / sum(n), 
                  subject = factor(subject), 
                  s = factor(s), 
                  d = factor(d)) |> 
    dplyr::ungroup() |> 
    ggplot() +
    geom_point(aes(x = block, y = freq, group = subject, col = d), alpha = 0.2) +
    geom_line(aes(x = block, y = freq, group = interaction(subject, d), col = d), alpha = 0.2) +
    geom_line(aes(x = block, y = freq, col = d), size = 1.5, data = data_aggr) +
    facet_wrap( ~ s, nrow = 2, ncol = 2) +
    labs(x = 'block', y = 'probability') +
    scale_color_brewer(name = "", palette = 'Set1') +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks=seq(1, max(data$block), by = 1)) +
    theme(legend.position = "top")
}


#' Descriptive plots 
#'
#' Plot the mean response times of the raw data.
#' 
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @return Population level raw response times
plot_RT = function(data){
  
  mean_r_time = block = d = s = cens = r_time = NULL
  
  data |>
    dplyr::filter(cens == 0) |>
    dplyr::group_by(block, s, d) |>
    dplyr::summarise(mean_r_time = mean(r_time)) |>
    dplyr::ungroup() |> 
    ggplot() +
    geom_line(aes(x = block, y = mean_r_time, col = factor(d)), size = 1.5) +
    facet_wrap( ~ s) +
    labs(x = 'block', y = 'response time') +
    scale_color_brewer(name = "", palette = "Set1") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks=seq(1, max(data$block), by = 1)) +
    theme(legend.position = "top")
}




# The Half Cauchy Distribution
# 
# Probability density function for the half Cauchy distribution
# 
# @param x vector of points where to evaluate the pdf
# @param sigma vector of scale parameters
# @param log should the result be returned on the log scale?
# @return Probability density function for the half Cauchy distribution with location 0 and scale sigma
dhalfcauhy = function(x, sigma, log = T){
  
  out = log(2) + log(sigma) - log(pi) - log(sigma^2 + x^2)
  if (!log){
    out = exp(out)
  }
  
  return(out)
}


#' Spline Penalty Matrix
#' 
#' Construct the covariance matrix P of the smoothness inducing prior for the
#' spline coefficients
#' 
#' @param K Number of spline knots
#' @return Covariance of the smoothness inducing prior (penalizing first
#' differences in the spline coefficients)
P_smooth1 = function(K){
  D = diag(rep(1,K))
  D = diff(D)
  P = crossprod(D)
  
  return(P)
}


# Locally Informed Moves
# 
# Computes the g() function necessary for locally informed moves.
# 
# @param x input of the function
# @param log should the result be computed on a log scale?
# @return g(x)
g_HB = function(x, log = T){
  # In the case we do not want to use the locally informed version, 
  # set out = 0
  out = x/2
  
  if (log){
    return(out)
  }
  else{
    return(exp(out))
  }
}


#' Spline Basis Functions
#'
#' Construct the J basis functions for the splines evaluated on a grid.
#' 
#' @param xgrid grid where we want to evaluate the spline functions (vector of length n)
#' @param knots vector of knots for the splines (vector of length K)
#' @return n x (K+1) - matrix representing the value of each basis function evaluated on xgrid
B_basis = function(xgrid, knots){
  n = length(xgrid) # number of grid points where to evaluate the spline
  K = length(knots) # number of knots
  delta = knots[2] - knots[1]
  
  B = array(0, dim = c(n, K + 1))
  for (j in 1:(K-1)){
    act_idx = (1:n)[(xgrid >= knots[j]) & (xgrid <= knots[j+1])]
    act_x = xgrid[act_idx]
    resc_x = (act_x - knots[j]) / (knots[j+1] - knots[j])
    
    B[act_idx,j] = (1/2) * (1 - resc_x)^2
    B[act_idx,j+1] = -(resc_x^2) + resc_x + 1/2
    B[act_idx,j+2] = (resc_x^2)/2
  }
  
  return(B)
}


#' Hamming Ball
#'
#' Computes the Hamming Ball centered at x with radius r.
#' 
#' @param x center of the Hamming Ball
#' @param S number of states 
#' @param r radius of the Hamming Ball
#' @return Hamming Ball
H_ball = function(x, S, r){
  K = length(x) # number of chains
  card_m = (S - 1)^(0:r) * choose(K, 0:r) 
  # sum(card_m) is the cardinality of each HB with radius r
  
  # First, put the current vector in its HB
  HB_temp = matrix(x, K, 1)
  
  for (m in 1:r){ # loop the possible HB radius
    HB_diff = matrix(rep(x, card_m[m+1]), K, card_m[m+1])
    index = utils::combn(K, m) # what elements of the vector to change?
    
    for (j in 1:ncol(index)){
      vec = NULL
      for (i in 1:nrow(index)){ 
        vec = c(vec, (1:S)[-x[index[i,j]]])
      }
      prop_col = t(gtools::permutations(n = length(unique(vec)), r = m, 
                                        v = vec, repeats.allowed = TRUE))
      keep_col = prop_col[,which(colSums(prop_col != x[index[,j]]) == m)]
      HB_diff[index[,j], ((j-1)*(card_m[m+1]/ncol(index))+1):(j*card_m[m+1]/ncol(index))] = keep_col
    }
    HB_temp = cbind(HB_temp,HB_diff) # save these in the Hamming Ball
  }
  
  return (HB_temp)
}


# Hamming Ball Sampler
#
# Samples a configuration uniformly within the HB.
# 
# @param x center of the Hamming Ball
# @param S number of states 
# @param r radius of the Hamming Ball
# @return sampled state
H_ball_unif = function(x, S, r){
  HB_temp = H_ball(x, S, r)
  HB_samp = HB_temp[,sample(1:ncol(HB_temp), 1)]
  
  return (HB_samp)
}


# Random effects hyperparameters update
#
# Metropolis Hastings step to update the variance parameter and the smoothness
# parameter for the random effects.
#
# @param sigma2_ua_old variance parameter at the previous iteration
# @param sigma2_us_old smoothness parameter at the previous iteration
# @param beta_u_old random effects at the current iteration
# @param P_smooth matrix penalizing the first order differences of the coefficients
# @param n_ind number of participants (random effects)
# @return sampled variance and smoothness parameters
sample_smooth_var = function(sigma2_ua_old, sigma2_us_old, 
                             beta_u_old, P_smooth, n_ind){
  
  J = ncol(beta_u_old)
  # Update \sigma_{u,a}^{2}
  sigma2_ua_prop = exp(rnorm(1, mean = log(sigma2_ua_old), sd = 0.2))
  lu = log(runif(1))
  log_prop = 0.5 * n_ind * log(det(diag(J)/sigma2_ua_prop + P_smooth/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(diag(J)/sigma2_ua_prop, beta_u_old))) -
    log(1 + sigma2_ua_prop^2)
  log_old = 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + P_smooth/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(diag(J)/sigma2_ua_old, beta_u_old))) -
    log(1 + sigma2_ua_old^2)
  alpha = min(c(0, log_prop + log(sigma2_ua_prop) -
                  log_old - log(sigma2_ua_old)))
  if (lu < alpha){
    sigma2_ua_old = sigma2_ua_prop
  }
  
  # Update \sigma_{u,s}^{2}
  sigma2_us_prop = exp(rnorm(1, mean = log(sigma2_us_old), sd = 0.1))
  while (sigma2_us_prop > 0.2){
    sigma2_us_prop = exp(rnorm(1, mean = log(sigma2_us_old), sd = 0.1))
  }
  lu = log(runif(1))
  log_prop = 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + P_smooth/sigma2_us_prop)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(P_smooth/sigma2_us_prop, beta_u_old))) -
    log(1 + sigma2_us_prop^2)
  log_old = 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + P_smooth/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(P_smooth/sigma2_us_old, beta_u_old))) -
    log(1 + sigma2_us_old^2)
  alpha = min(c(0, log_prop + log(sigma2_us_prop) -
                  log_old - log(sigma2_us_old)))
  if (lu < alpha){
    sigma2_us_old = sigma2_us_prop
  }
  
  return(list('sigma2_us_old' = sigma2_us_old, 
              'sigma2_ua_old' = sigma2_ua_old))
  
}


#' Drift Diffusion Model Fit
#'
#' Main function for the Gibbs sampler for the drift-diffusion model. Note that 
#' priors are noninformative and calibrated so that, for the most stable 
#' performance, the response times (variable `r_time` in the `data` dataframe) 
#' should lie between 0 and 10.
#'
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times. To avoid numerical 
#'           issues, the unit of measurement should be such that the numerical 
#'           values of most response times should lie between 0 and 10
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @param hypers hyperparameters of the MCMC: list containing "s_sigma_mu" and "s_sigma_b", 
#'               which are the smoothness parameters for drifts and boundaries, respectively)
#' @param boundaries whether to fit the unrestricted model (flexible), assume constant 
#'                   boundaries over time (constant) or fix the boundaries to the same level 
#'                   across predictors (fixed)
#' @param Niter total number of iterations
#' @param burnin burnin of the chain
#' @param thin thinning factor
#' @return List with the following MCMC posterior samples: 
#' * post_mean_delta: posterior samples for the population offset parameters
#' * post_mean_mu: posterior samples for the population drift parameters
#' * post_mean_b: posterior samples for the population boundary parameters
#' * post_ind_delta: posterior samples for the individual offset parameters
#' * post_ind_mu: posterior samples for the individual drift parameters
#' * post_ind_b: posterior samples for the individual boundary parameters
#' * sigma2_mu_us: posterior samples for the random effects drift smoothness parameters
#' * sigma2_mu_ua: posterior samples for the random effects drift variance parameters
#' * sigma2_b_us: posterior samples for the random effects boundary smoothness parameters
#' * sigma2_b_ua: posterior samples for the random effects boundary variance parameters
#' * sigma2_1_mu: posterior samples for the drift smoothness parameters
#' * sigma2_1_b: posterior samples for the boundary smoothness parameters
#' * pred_ans: predicted population-level categories
#' * pred_time: predicted population-level response times
#' * pred_ans_ind: predicted individual-level categories
#' * pred_time_ind: predicted individual-level response times
LDDMM = function(data, hypers, boundaries = 'flexible', 
                 Niter = 5000, burnin = 2000, thin = 5){
  
  # Check for data issues
  if ( (is.null(hypers[["s_sigma_mu"]])) | (is.null(hypers[["s_sigma_b"]]))){
    cat('The hyperparameters did not contain s_sigma_mu or s_sigma_b, running with their default values.\n\n')
    hypers$s_sigma_mu =  hypers$s_sigma_b = 0.1
  }
  if ( ("subject" %notin% colnames(data)) | 
       ("block" %notin% colnames(data)) | 
       ("s" %notin% colnames(data)) | 
       ("d" %notin% colnames(data)) | 
       ("r_time" %notin% colnames(data)) | 
       ("cens" %notin% colnames(data))){
    stop('The data is not in the correct format!')
  }
  
  # Call one of the main functions
  if (boundaries == 'flexible'){
    fit = LDDMM_full(data, hypers, Niter, burnin, thin)
  }
  else if (boundaries == 'constant') {
    fit = LDDMM_const_bound(data, hypers, Niter, burnin, thin)
  }
  else if (boundaries == 'fixed') {
    fit = LDDMM_fix_bound(data, hypers, Niter, burnin, thin)
  }
  else if (boundaries == 'fixed-constant') {
    fit = LDDMM_fixandconst_bound(data, hypers, Niter, burnin, thin)
  }
  else {
    stop("The argument boundaries can be only one of the following: flexible, constant, fixed, or fixed-constant")
  }
  
  return (fit)
}


LDDMM_full = function(data, hypers, cluster = TRUE, Niter = 5000, burnin = 2000, thin = 5){
  
  # Choose the number of knots (default is between the beginning and the end of 
  # the study, at every block)
  T_min = min(data$block)
  T_max = max(data$block)
  knots = T_min:T_max
  K = length(knots)
  P_smooth = P_smooth1(K + 1)
  
  # Choose a fine grid in order to plot the curves resulting from the spline basis
  # expansion
  xgrid = seq(T_min, T_max, by = 1)
  
  # Extract quantities of interest
  tau = data$r_time
  ind = data$subject
  time = data$block
  cens = data$cens
  D = cbind(data$s, data$d)
  
  n_ind = length(unique(ind))
  ind = as.numeric(plyr::mapvalues(factor(ind),
                                   from = levels(factor(ind)),
                                   to = 1:n_ind))
  B = B_basis(data$block, knots)
  Bgrid = B_basis(xgrid, knots)
  
  samp_size = (Niter - burnin)/thin # sample size
  p = ncol(D) # number of covariates
  d_j = rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] = length(unique(D[!is.na(D[,j]),j]))
  }
  J = ncol(B) # number of locations
  n = nrow(B) # total number of observations
  n_ind = length(unique(ind)) # number of individuals
  
  
  # Rescale the time steps in \{1, ..., T_max\}
  T_max = max(time)
  T_min = min(time)
  time = time - T_min + 1
  T_max = max(time)
  T_min = min(time)
  idx_xy = t(apply(expand.grid(y = 1:d_j[2], x = 1:d_j[1]), 1, rev))
  colnames(idx_xy) = NULL
  
  
  # Set HB structures
  r_HB = 1
  Z_max = if_else(cluster, min(prod(d_j), 6), prod(d_j))
  dim_HB = sum((Z_max - 1)^(0:r_HB) * choose(d_j[2], 0:r_HB))
  
  # Set MCMC objects
  z = array(NA, dim = c(prod(d_j), J, samp_size))
  post_mean_delta = array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_mean_b = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_ind_delta = array(NA, dim = c(d_j[1], n_ind, samp_size))
  post_ind_mu = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  post_ind_b = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  sigma2_mu_us = array(NA, dim = samp_size)
  sigma2_mu_ua = array(NA, dim = samp_size)
  sigma2_b_us = array(NA, dim = samp_size)
  sigma2_b_ua = array(NA, dim = samp_size)
  sigma2_1_mu = array(NA, dim = samp_size)
  sigma2_1_b = array(NA, dim = samp_size)
  pred_time = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  loglik_chain = array(NA, dim = c(samp_size, n))
  
  # Set initial values
  delta_old = array(NA, dim = c(d_j[1], n_ind))
  beta_mu_old = array(NA, dim = c(J, d_j))
  beta_b_old = array(NA, dim = c(J, d_j))
  delta_dat = array(NA, dim = n)
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
      if (length(idx_temp) > 0){
        delta_old[s_temp,i_temp] = min(tau[idx_temp])/2
        delta_dat[idx_temp] = delta_old[s_temp,i_temp]
      }
      else{
        delta_old[s_temp,i_temp] = 1E-3
      }
    }
    for(j in 1:d_j[1]){
      beta_mu_old[,s_temp,j] = rep(0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
      beta_b_old[,s_temp,j] = rep(1.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
    }
  }
  
  
  low_bound_mu = min(beta_mu_old) - 1.5
  upp_bound_mu = max(beta_mu_old) + 1
  low_bound_b = min(beta_b_old) - 1.5
  upp_bound_b = max(beta_b_old) + 1
  
  beta_mu_star_prop = array(NA, dim = c(J, Z_max))
  beta_mu_u_old = array(0, dim = c(n_ind, J))
  beta_b_star_prop = array(NA, dim = c(J, Z_max))
  beta_b_u_old = array(0, dim = c(n_ind, J))
  sigma2_1_mu_old = 0.005
  sigma2_1_b_old = 0.005
  sigma2_b_us_old = 0.005
  sigma2_mu_us_old = 0.005
  sigma2_b_ua_old = 0.005
  sigma2_mu_ua_old = 0.005
  
  
  # Message passing structures
  beta_mess = array(NA, dim = c(J, dim_HB))
  z_old = list()
  
  for (j in 1:d_j[1]){
    z_old[[j]] = array(NA, dim = c(d_j[2], J))
    for (jj in 1:d_j[2]){
      if (j == jj){
        z_old[[j]][jj,] = if_else(cluster, j, jj + (j - 1) * d_j[1])
      }
      else{
        z_old[[j]][jj,] = if_else(cluster, sample((d_j[1] + 1):Z_max, 1), jj + (j - 1) * d_j[1])
      }
    }
  }
  z_temp = do.call(rbind, z_old)
  
  
  beta_mu_star_old = array(NA, dim = c(J, Z_max))
  beta_b_star_old = array(NA, dim = c(J, Z_max))
  for (i in 1:Z_max){
    beta_mu_star_old[,i] = beta_mu_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                       idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]] 
    beta_b_star_old[,i] = beta_b_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                     idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]] 
  }
  
  rand_mat = array(rnorm(prod(d_j)), dim = d_j)
  idx_succ = which(rand_mat == diag(rand_mat))
  
  v_old = array(NA, dim = c(d_j[2], J))
  z_prop = list()
  
  # Transition dynamics objects
  alpha_S_old = alpha_F_old = 1
  Q_S_old = array(NA, dim = c(Z_max, Z_max))
  Q_F_old = array(NA, dim = c(Z_max, Z_max))
  pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                        table_int(z_temp[idx_succ,1], Z_max))
  pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                        table_int(z_temp[-idx_succ,1], Z_max))
  tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
  tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
  for (h in 1:Z_max){
    Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + tr_count_S[h,])
    Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + tr_count_F[h,])
  }
  
  
  # Auxiliary variables
  prob_mat = array(NA, dim = c(dim_HB, dim_HB))
  prob_vec = array(NA, dim = dim_HB)
  
  # MH proposal parameters
  sd_MH_delta = array(0.3, dim = c(d_j[1], n_ind))
  sd_MH_beta_mu = 0.4
  sd_MH_beta_b = array(0.25, dim = c(J, Z_max))
  sd_beta_mu_u =  array(0.4, dim = c(n_ind, J))
  sd_beta_b_u =  array(0.4, dim = c(n_ind, J))
  acc_delta = array(0, dim = c(d_j[1], n_ind))
  acc_beta_mu = array(0, dim = c(J, Z_max))
  acc_beta_b = array(0, dim = c(J, Z_max))
  acc_beta_mu_u = array(0, dim = c(n_ind, J))
  acc_beta_b_u = array(0, dim = c(n_ind, J))
  n_batch = 0
  
  # Auxiliary variables
  B_beta_mu_dat = array(0, dim = c(n, d_j[1]))
  B_beta_mu_u_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_u_dat = array(0, dim = c(n, d_j[1]))
  mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
  b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
  
  
  # Gibbs Sampler
  it = 1
  pb = txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (0) Adaptively tune the MH variance for the proposals of \delta_{s, i}, 
    # beta_u_mu, beta_u_b
    if (iter %% 20 == 0){
      n_batch = n_batch + 1
      delta_n = min(0.01, n_batch^(-0.5))
      
      for (i in 1:n_ind){
        for (x_temp in 1:d_j[1]){
          if (acc_delta[x_temp,i]/iter > 0.44){
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) + delta_n)
          }
          else{
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) - delta_n)
          }
        }
        for (k in 1:J){
          if (acc_beta_mu_u[i,k]/iter > 0.44){
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) + delta_n)
          }
          else{
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) - delta_n)
          }
          if (acc_beta_b_u[i,k]/iter > 0.44){
            sd_beta_b_u[i,k] = exp(log(sd_beta_b_u[i,k]) + delta_n)
          }
          else{
            sd_beta_b_u[i,k] = exp(log(sd_beta_b_u[i,k]) - delta_n)
          }
        }
      }
    }
    
    
    # (1) Update of the delta parameter: \delta_{s,i}: MH with log normal 
    #     proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp = tau[idx_temp]
        cens_temp = cens[idx_temp]
        D_temp = D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        while (delta_prop > min(tau[idx_temp])){
          delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        }
        loglik_prop = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                     b_dat[idx_temp,], 
                                     rep(delta_prop, length(idx_temp)), 
                                     cens_temp, D_temp, TRUE)
        loglik_old = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                    b_dat[idx_temp,], 
                                    rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                    cens_temp, D_temp, TRUE)
        
        alpha_acc = min(0, loglik_prop + log(delta_prop) -
                          loglik_old - log(delta_old[s_temp,i_temp]))
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] = delta_prop
          delta_dat[idx_temp] = delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] = acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    # (2) Joint update of b, mu parameters: \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time = T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time = T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time = (k-1):k
      }
      
      for (h in 1:Z_max){ # loop over possible latent values
        # tuples (combinations of covariates) that are clustered together via 
        # the latent h
        idx_cov = matrix(idx_xy[which(z_temp[,k] == h),], length(which(z_temp[,k] == h)), p)
        
        X_1k = unique(idx_cov[,1]) # all possible values of x in this cluster
        X_2k = unique(idx_cov[,2]) # all possible values of y in this cluster
        
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{j,k}: posterior update
          # Pick data with covariate levels of x clustered in group h and 
          # at the correct locations
          idx_i = which( (D[,1] %in% X_1k) & (time %in% idx_time) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[which((D[,1] %in% X_1k) & (time %in% idx_time)),]
          
          
          # Normal proposal distribution centered on the current value
          if (k == 1){
            beta_b_star_prop[,h] = beta_b_star_old[,h]
            beta_b_star_prop[k,h] = rnorm(1, beta_b_star_old[k,h], sd_MH_beta_b[k,h])
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_mu_star_old[k,h], sd_MH_beta_mu)
          }
          # Normal proposal from the prior
          else if (k == J){
            beta_pre = beta_b_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_b_star_prop[,h] = beta_b_star_old[,h]
            beta_b_star_prop[k,h] = rnorm(1, beta_pre, sqrt(sigma2_1_b_old))
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          # Normal proposal from the prior
          else {
            beta_pre = beta_b_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_b_star_prop[,h] = beta_b_star_old[,h]
            beta_b_star_prop[k,h] = rnorm(1, beta_pre, sqrt(sigma2_1_b_old))
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          B_beta_b_prop_dat = B_beta_b_dat[idx_i,]
          B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
          
          # Modify the proposed values in the corresponding positions
          for (hh in 1:nrow(idx_cov)){
            B_beta_b_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] = 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_b_star_prop[,h]
            B_beta_mu_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] = 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_mu_star_prop[,h]
          }
          
          # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
          b_prop_dat = exp(B_beta_b_prop_dat + B_beta_b_u_dat[idx_i,])
          mu_prop_dat = exp(B_beta_mu_prop_dat + B_beta_mu_u_dat[idx_i,])
          
          
          if (k == 1){
            beta_post = beta_b_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_prop_dat, delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_prop[k,h] - beta_post)^2) - 
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_old[k,h] - beta_post)^2) - 
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          else if (k == J){
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_prop_dat, delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE)
          }
          else {
            beta_post = beta_b_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_prop_dat, delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_prop[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_old[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          
          alpha_acc = min(0, logpost_prop - logpost_old)
          l_u = log(runif(1))
          
          if (l_u < alpha_acc){
            beta_b_star_old[k,h] = beta_b_star_prop[k,h]
            B_beta_b_dat[idx_i,] = B_beta_b_prop_dat
            b_dat[idx_i,] = b_prop_dat
            acc_beta_b[k,h] = acc_beta_b[k,h] + 1
            
            beta_mu_star_old[k,h] = beta_mu_star_prop[k,h]
            B_beta_mu_dat[idx_i,] = B_beta_mu_prop_dat
            mu_dat[idx_i,] = mu_prop_dat
            acc_beta_mu[k,h] = acc_beta_mu[k,h] + 1
          }
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_b_star_old[k,h] = runif(1, low_bound_b, upp_bound_b)
          beta_mu_star_old[k,h] = runif(1, low_bound_mu, upp_bound_mu)
        }
      }
    }
    
    
    # (3) Update the cluster assignments
    for (x_temp in 1:d_j[1]){ # loop over possible latent values
      if (cluster) {
        beta_mess = array(-Inf, dim = c(J, dim_HB))
        beta_mess[J,] = 1/dim_HB
        
        v_old[,J] = H_ball_unif(z_old[[x_temp]][,J], S = Z_max, r = r_HB)
        z_prop[[J]] = H_ball(v_old[,J], S = Z_max, r = r_HB)
        for (k in (J - 1):1){
          idx_i = which( (time == k) & (D[,1] == x_temp) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[idx_i,]
          
          # (i) Sample the auxiliary variables
          v_temp = H_ball(z_old[[x_temp]][,k], S = Z_max, r = r_HB)
          
          probs = rep(-Inf, dim_HB)
          for (h in 1:dim_HB){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
            B_beta_b_prop_dat = B_beta_b_dat[idx_i,]
            
            B_beta_mu_prop_dat = 0.5 * (beta_mu_star_old[k,v_temp[,h]] +
                                          beta_mu_star_old[k+1,z_old[[x_temp]][,k+1]])
            B_beta_b_prop_dat = 0.5 * (beta_b_star_old[k,v_temp[,h]] +
                                         beta_b_star_old[k+1,z_old[[x_temp]][,k+1]])
            
            mu_dat_prop = exp(t(B_beta_mu_prop_dat + t(B_beta_mu_u_dat[idx_i,])))
            b_dat_prop = exp(t(B_beta_b_prop_dat + t(B_beta_b_u_dat[idx_i,])))
            
            probs[h] = g_HB(log_likelihood(tau_temp, mu_dat_prop, b_dat_prop,
                                           delta_dat[idx_i], cens_temp, D_temp, TRUE))
          }
          probs = as.numeric(normalise_log(probs))
          
          v_old[,k] = v_temp[,sample(1:dim_HB, 1, prob = probs)]
          z_prop[[k]] = H_ball(v_old[,k], S = Z_max, r = r_HB)
          
          
          # (ii) Pass messages backwards only in the restricted state space given
          #      by the slice
          z_kp1_temp = which(beta_mess[k+1,] > 0)
          prob_mat = array(-Inf, dim = c(dim_HB, dim_HB))
          
          for (h1 in z_kp1_temp){
            for (h2 in 1:dim_HB){
              
              B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
              B_beta_b_prop_dat = B_beta_b_dat[idx_i,]
              B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k,z_prop[[k]][,h2]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_b_prop_dat_1 = 0.5 * (beta_b_star_old[k,z_prop[[k]][,h2]] +
                                             beta_b_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k,v_old[,k]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_b_prop_dat_2 = 0.5 * (beta_b_star_old[k,v_old[,k]] +
                                             beta_b_star_old[k+1,z_prop[[k+1]][,h1]])
              
              mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_i,])))
              b_dat_prop_1 = exp(t(B_beta_b_prop_dat_1 + t(B_beta_b_u_dat[idx_i,])))
              mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_i,])))
              b_dat_prop_2 = exp(t(B_beta_b_prop_dat_2 + t(B_beta_b_u_dat[idx_i,])))
              
              prob_mat[h2,h1] = log(beta_mess[k+1,h1]) -
                0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k,z_prop[[k]][,h2]] -
                                               beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])^2) -
                0.5 / sigma2_1_b_old * sum((beta_b_star_old[k,z_prop[[k]][,h2]] -
                                              beta_b_star_old[k+1,z_prop[[k+1]][,h1]])^2) +
                log_likelihood(tau_temp, mu_dat_prop_1, b_dat_prop_1,
                               delta_dat[idx_i], cens_temp, D_temp, TRUE) +
                g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat_prop_2,
                                    delta_dat[idx_i], cens_temp, D_temp, TRUE)) +
                sum(log(Q_F_old[cbind(z_prop[[k]][-x_temp,h2],z_prop[[k+1]][-x_temp,h1])])) +
                log(Q_S_old[z_prop[[k]][x_temp,h2],z_prop[[k+1]][x_temp,h1]])
            }
          }
          if ( sum(is.infinite(sum_rows_log(prob_mat))) == dim_HB){
            beta_mess[k,] = 1/dim_HB
          }
          else{
            beta_mess[k,] = as.numeric(sum_rows_log(prob_mat))
            beta_mess[k,] = as.numeric(normalise_log(beta_mess[k,]))
          }
        }
        
        
        # (iii) Sample states forward (only on allowed states)
        idx_fail = (1:d_j[2])[-x_temp]
        # Sample z_1
        prob_vec = log(beta_mess[1,]) + log(pi_S_0[z_prop[[1]][x_temp,]]) +
          colSums(matrix(log(pi_F_0[z_prop[[1]][-x_temp,]]), d_j[2] - 1, dim_HB))
        prob_vec = as.numeric(normalise_log(prob_vec))
        
        idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
        z_old[[x_temp]][,1] = z_prop[[1]][,idx_samp]
        
        # Sample z_k
        for (k in 2:J){
          idx_km1 = which( (time == k - 1) & (D[,1] == x_temp) )
          tau_temp = tau[idx_km1]
          cens_temp = cens[idx_km1]
          D_temp = D[idx_km1,]
          
          prob_vec = log(beta_mess[k,]) + 
            log(Q_S_old[cbind(z_old[[x_temp]][x_temp,k-1], z_prop[[k]][x_temp,])])
          for (kkk in idx_fail){
            prob_vec = prob_vec + log(Q_F_old[cbind(z_old[[x_temp]][kkk,k-1], z_prop[[k]][kkk,])])
          }
          
          for (z_k_temp in which(is.finite(prob_vec))){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_km1,]
            B_beta_b_prop_dat = B_beta_b_dat[idx_km1,]
            
            B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_b_prop_dat_1 = 0.5 * (beta_b_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                           beta_b_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k-1,v_old[,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_b_prop_dat_2 = 0.5 * (beta_b_star_old[k-1,v_old[,k-1]] +
                                           beta_b_star_old[k,z_prop[[k]][,z_k_temp]])
            
            mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_km1,])))
            b_dat_prop_1 = exp(t(B_beta_b_prop_dat_1 + t(B_beta_b_u_dat[idx_km1,])))
            mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_km1,])))
            b_dat_prop_2 = exp(t(B_beta_b_prop_dat_2 + t(B_beta_b_u_dat[idx_km1,])))
            
            prob_vec[z_k_temp] = prob_vec[z_k_temp] -
              0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                             beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])^2) -
              0.5 / sigma2_1_b_old * sum((beta_b_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                            beta_b_star_old[k,z_prop[[k]][,z_k_temp]])^2) +
              log_likelihood(tau_temp, mu_dat_prop_1, b_dat_prop_1,
                             delta_dat[idx_km1], cens_temp, D_temp, TRUE) +
              g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat_prop_2,
                                  delta_dat[idx_km1], cens_temp, D_temp, TRUE))
          }
          prob_vec = as.numeric(normalise_log(prob_vec))
          
          idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
          z_old[[x_temp]][,k] = z_prop[[k]][,idx_samp]
        }
      }
      
      # (4) Assign the cluster specific curves f_{\mu}
      for (y_temp in 1:d_j[2]){
        beta_mu_old[,x_temp,y_temp] = beta_mu_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
        beta_b_old[,x_temp,y_temp] = beta_b_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
      }
      B_beta_mu_dat[(D[,1] == x_temp),] = B[D[,1] == x_temp,] %*% beta_mu_old[,x_temp,]
      B_beta_b_dat[(D[,1] == x_temp),] = B[D[,1] == x_temp,] %*% beta_b_old[,x_temp,]
    }
    z_temp = do.call(rbind, z_old)
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
    
    
    
    # (5) Update the transition probabilities
    pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + table_int(z_temp[idx_succ,1], Z_max))
    pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + table_int(z_temp[-idx_succ,1], Z_max))
    tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
    tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
    for (h in 1:Z_max){
      Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                                 tr_count_S[h,])
      Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                                 tr_count_F[h,])
    }
    
    
    alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    while ((alpha_S_prop < 0.01) | (alpha_S_prop > 10)){
      alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    }
    alpha_acc = dgamma(alpha_S_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_S_prop) -
      Z_max^2 * lgamma(alpha_S_prop/Z_max) +
      (alpha_S_prop/Z_max - 1) * sum(log(Q_S_old)) - (
        dgamma(alpha_S_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_S_old) -
          Z_max^2 * lgamma(alpha_S_old/Z_max) +
          (alpha_S_old/Z_max - 1) * sum(log(Q_S_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_S_old = alpha_S_prop
      }
    }
    
    alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    while ((alpha_F_prop < 0.01) | (alpha_F_prop > 10)){
      alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    }
    alpha_acc = dgamma(alpha_F_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_F_prop) -
      Z_max^2 * lgamma(alpha_F_prop/Z_max) +
      (alpha_F_prop/Z_max - 1) * sum(log(Q_F_old)) - (
        dgamma(alpha_F_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_F_old) -
          Z_max^2 * lgamma(alpha_F_old/Z_max) +
          (alpha_F_old/Z_max - 1) * sum(log(Q_F_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_F_old = alpha_F_prop
      }
    }
    
    
    # (6) Correction term for random effects.
    corr_mu = colMeans(beta_mu_u_old)
    corr_b = colMeans(beta_b_u_old)
    beta_mu_u_old = t(t(beta_mu_u_old) - corr_mu)
    beta_b_u_old = t(t(beta_b_u_old) - corr_b)
    
    for (k in 1:J){
      if (k == 1){
        idx_time = which(time == T_min)
      }
      else if (k == J){
        idx_time = which(time == T_max)
      }
      else {
        idx_time = which(time %in% (k-1):k)
      }
      beta_mu_old[k,,] = beta_mu_old[k,,] + corr_mu[k]
      beta_b_old[k,,] = beta_b_old[k,,] + corr_b[k]
      beta_mu_star_old[k,] = beta_mu_star_old[k,] + corr_mu[k]
      beta_b_star_old[k,] = beta_b_star_old[k,] + corr_b[k]
      B_beta_mu_dat[idx_time,] = B_beta_mu_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_mu_u_dat[idx_time,] = B_beta_mu_u_dat[idx_time,] - 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_b_dat[idx_time,] = B_beta_b_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_b)
      B_beta_b_u_dat[idx_time,] = B_beta_b_u_dat[idx_time,] -
        as.numeric(B[idx_time,] %*% corr_b)
    }
    
    
    # (7) Update the cluster specific smoothness parameters
    RSS_mu = RSS_b = 0
    for (h1 in 1:d_j[1]){
      for (h2 in 1:d_j[2]){
        RSS_mu = RSS_mu + as.numeric(crossprod(beta_mu_old[,h1,h2], P_smooth) %*% 
                                       beta_mu_old[,h1,h2])
        RSS_b = RSS_b + as.numeric(crossprod(beta_b_old[,h1,h2], P_smooth) %*% 
                                     beta_b_old[,h1,h2])
      }
    }
    
    sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    while(sigma2_1_mu_temp > 0.2){
      sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_mu_temp) - 
                      0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_temp) - 
                      0.5/sigma2_1_mu_temp * RSS_mu + 
                      dhalfcauhy(sqrt(sigma2_1_mu_temp), hypers$s_sigma_mu, T) -
                      (log(sigma2_1_mu_old) - 
                         0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_old) -
                         0.5/sigma2_1_mu_old * RSS_mu +
                         dhalfcauhy(sqrt(sigma2_1_mu_old), hypers$s_sigma_mu, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_mu_old = sigma2_1_mu_temp
    }
    
    sigma2_1_b_temp = exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    while(sigma2_1_b_temp > 0.2){
      sigma2_1_b_temp = exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_b_temp) - 
                      0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_b_temp) - 
                      0.5/sigma2_1_b_temp * RSS_b + 
                      dhalfcauhy(sqrt(sigma2_1_b_temp), hypers$s_sigma_b, T) -
                      (log(sigma2_1_b_old) - 
                         0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_b_old) -
                         0.5/sigma2_1_b_old * RSS_b +
                         dhalfcauhy(sqrt(sigma2_1_b_old), hypers$s_sigma_b, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_b_old = sigma2_1_b_temp
    }
    
    
    # (8a) Update random effects for b
    reff = sample_reff_b(tau, D, cens, beta_b_u_old, delta_dat, B_beta_b_dat,
                         b_dat, mu_dat, B, P_smooth, ind, time,
                         sigma2_b_us_old, sigma2_b_ua_old, sd_beta_b_u, 
                         acc_beta_b_u)
    beta_b_u_old = reff$beta_u_old
    B_beta_b_u_dat = reff$B_beta_u_dat
    b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
    acc_beta_b_u = reff$acc_beta_u
    
    
    # (8b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_b_ua_old, sigma2_b_us_old,
                               beta_b_u_old, P_smooth, n_ind)
    sigma2_b_ua_old = ls_var$sigma2_ua_old
    sigma2_b_us_old = ls_var$sigma2_us_old
    
    
    # (9a) Update random effects for mu
    reff = sample_reff_mu(tau, D, cens, beta_mu_u_old, delta_dat, b_dat,
                          B_beta_mu_dat, mu_dat, B, P_smooth, ind, time,
                          sigma2_mu_us_old, sigma2_mu_ua_old, sd_beta_mu_u,
                          acc_beta_mu_u)
    beta_mu_u_old = reff$beta_u_old
    B_beta_mu_u_dat = reff$B_beta_u_dat
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    acc_beta_mu_u = reff$acc_beta_u
    
    
    # (9b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_mu_ua_old, sigma2_mu_us_old,
                               beta_mu_u_old, P_smooth, n_ind)
    sigma2_mu_ua_old = ls_var$sigma2_ua_old
    sigma2_mu_us_old = ls_var$sigma2_us_old
    
    
    
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # This is the correction for the random effects integration: we need to 
      # compute the variance of the random effects as detailed in the Supplementary
      # Materials
      cov_reff_mu = solve(diag(J) / sigma2_mu_ua_old + P_smooth / sigma2_mu_us_old)
      corr_term_gr_mu = rep(0, nrow(Bgrid))
      corr_term_mu = rep(0, T_max)
      cov_reff_b = solve(diag(J) / sigma2_b_ua_old + P_smooth / sigma2_b_us_old)
      corr_term_gr_b = rep(0, nrow(Bgrid))
      corr_term_b = rep(0, T_max)
      for (k in 1:J){
        corr_term_gr_mu = corr_term_gr_mu + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_mu[,k]))
        corr_term_mu = corr_term_mu + rowSums(B_basis(1:T_max, knots)[,k] * 
                                                t(t(B_basis(1:T_max, knots)) * cov_reff_mu[,k]))
        corr_term_gr_b = corr_term_gr_b + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_b[,k]))
        corr_term_b = corr_term_b + rowSums(B_basis(1:T_max, knots)[,k] * 
                                              t(t(B_basis(1:T_max, knots)) * cov_reff_b[,k]))
      }
      
      
      
      for (h1 in 1:d_j[1]){
        for (h2 in 1:d_j[2]){
          post_mean_mu[,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + 0.5 * corr_term_gr_mu)
          post_mean_b[,h1,h2,it] = exp(Bgrid %*% beta_b_old[,h1,h2] + 0.5 * corr_term_gr_b)
          
          for (i in 1:n_ind){
            post_ind_mu[,i,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + Bgrid %*% beta_mu_u_old[i,])
            post_ind_b[,i,h1,h2,it] = exp(Bgrid %*% beta_b_old[,h1,h2] + Bgrid %*% beta_b_u_old[i,])
          }
        }
      }
      post_ind_delta[,,it] = delta_old
      
      
      # Sample from the predictive distributions of response category and 
      # response times (so we can use predictive checks to evaluate goodness of fit)
      for (t in 1:T_max){
        for (x_temp in 1:d_j[1]){
          mu_temp = as.numeric(exp(B_basis(t, knots) %*% beta_mu_old[,x_temp,] + 
                                     0.5 * corr_term_mu[t]))
          b_temp = as.numeric(exp(B_basis(t, knots) %*% beta_b_old[,x_temp,] + 
                                    0.5 * corr_term_b[t]))
          delta_temp = mean(delta_old[x_temp,])
          pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                               b_temp^2)
          pred_ans[t,x_temp,it] = which.min(pred_temp)
          pred_time[t,x_temp,it] = min(pred_temp)
          
          for (i in 1:n_ind){
            mu_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_mu_old[,x_temp,] + 
                                                              beta_mu_u_old[i,])))
            b_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_b_old[,x_temp,] + 
                                                             beta_b_u_old[i,])))
            delta_temp = delta_old[x_temp,i]
            pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                                 b_temp^2)
            pred_ans_ind[i,t,x_temp,it] = which.min(pred_temp)
            pred_time_ind[i,t,x_temp,it] = min(pred_temp)
          }
        }
      }
      
      
      # Save MCMC objects
      post_mean_delta[it,] = rowMeans(delta_old)
      sigma2_mu_us[it] = sigma2_mu_us_old
      sigma2_mu_ua[it] = sigma2_mu_ua_old
      sigma2_b_us[it] = sigma2_b_us_old
      sigma2_b_ua[it] = sigma2_b_ua_old
      sigma2_1_mu[it] = sigma2_1_mu_old
      sigma2_1_b[it] = sigma2_1_b_old
      loglik_chain[it,] = log_likelihood_ind(tau, mu_dat, b_dat, delta_dat, cens, D)
      # z[,,it] = z_temp		
      
      it = it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list(
    # 'Z' = z, 
    'post_mean_delta' = post_mean_delta, 
    'post_mean_mu' = post_mean_mu,
    'post_mean_b' = post_mean_b,
    'post_ind_delta' = post_ind_delta,
    'post_ind_mu' = post_ind_mu,
    'post_ind_b' = post_ind_b,
    'sigma2_mu_us' = sigma2_mu_us, 
    'sigma2_mu_ua' = sigma2_mu_ua,
    'sigma2_b_us' = sigma2_b_us, 
    'sigma2_b_ua' = sigma2_b_ua,
    'sigma2_1_mu' = sigma2_1_mu, 
    'sigma2_1_b' = sigma2_1_b, 
    'pred_ans' = pred_ans, 
    'pred_time' = pred_time,
    'pred_ans_ind' = pred_ans_ind, 
    'pred_time_ind' = pred_time_ind,
    'loglik' = loglik_chain
  ))
}


#' Parameter posterior means
#'
#' Function to extract the posterior means of the parameters of interest from a lddmm fit object.
#'
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @param fit fit from the lddmm function
#' @param par parameter to output ('drift', or 'boundary')
#' @return Matrix with the following columns: 
#' * subject: participant labels
#' * block: training blocks
#' * par_s_d, ...: posterior means for the requested parameters
extract_post_mean = function(data, fit, par = c('drift', 'boundary')){
  
  block = d = s = subject = value = NULL
  
  if (par == 'drift'){
    post_mean = apply(fit$post_ind_mu, 1:4, mean)
  }
  else if (par == 'boundary'){
    post_mean = apply(fit$post_ind_b, 1:4, mean)
  }
  
  T_min = min(data$block)
  T_max = max(data$block)
  
  estimated_pars = dplyr::tibble(reshape2::melt(post_mean, varnames = c('block', 'subject', 's', 'd'))) |> 
    dplyr::filter(s == d) |> 
    dplyr::mutate(block = (T_min:T_max)[block], 
                  subject = plyr::mapvalues(factor(subject), from = levels(factor(subject)), to = levels(factor(data$subject)))) |>  
    pivot_wider(names_from = c(s, d), 
                values_from = value, 
                names_prefix = paste0(par, '_'))
  
  return (estimated_pars)
}



#' Parameter posterior draws
#'
#' Function to extract the posterior draws of the parameters of interest from a lddmm fit object.
#'
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @param fit fit from the lddmm function
#' @param par parameter to output ('drift', or 'boundary')
#' @return Matrix with the following columns: 
#' * subject: participant labels
#' * block: training blocks
#' * draw: iteration of the MCMC estimates
#' * par_s_d, ...: posterior draws for the requested parameters
extract_post_draws = function(data, fit, par = c('drift', 'boundary')){
  
  block = d = s = subject = value = NULL
  
  if (par == 'drift'){
    post_draw = fit$post_ind_mu
  }
  else if (par == 'boundary'){
    post_draw = fit$post_ind_b
  }
  d1 = dim(post_draw)[3]
  samp_size = dim(post_draw)[5]
  
  T_min = min(data$block)
  T_max = max(data$block)
  
  draw_pars = dplyr::tibble(reshape2::melt(post_draw, varnames = c('block', 'subject', 's', 'd', 'draw'))) |> 
    dplyr::filter(s == d) |> 
    dplyr::mutate(block = (T_min:T_max)[block], 
                  subject = plyr::mapvalues(factor(subject), from = levels(factor(subject)), to = levels(factor(data$subject)))) |>  
    pivot_wider(names_from = c(s, d), 
                values_from = value, 
                names_prefix = paste0(par, '_'))
  
  return (draw_pars)
}


#' Plot posterior estimates
#'
#' Function to plot the posterior mean and credible intervals of the parameters of interest from a lddmm fit object.
#'
#' @param data dataframe with the following columns:
#' * subject: vector of size n containing the participant labels
#' * block: vector of size n containing the training blocks (longitudinal units)
#' * s: vector of size n containing the stimuli
#' * d: vector of size n containing the decisions
#' * r_time: vector of size n containing the response times
#' * cens: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#' @param fit fit from the lddmm function
#' @param par parameter to output ('drift', or 'boundary')
#' @return Posterior mean and 95% CI
plot_post_pars = function(data, fit, par = c('drift', 'boundary')){
  
  Var1 = Var3 = low = mean_r_time = upp = value = NULL 
  
  if (par == 'drift'){
    post_draw = fit$post_mean_mu
    tex_string = "$\\mu_{d,s}(t)$"
  }
  else if (par == 'boundary'){
    post_draw = fit$post_mean_b
    tex_string = "$\\b_{d,s}(t)$"
  }
  
  K = max(data$block) - min(data$block) + 1
  xgrid = seq(min(data$block), max(data$block), by = 1)
  
  post_mean = reshape2::melt(apply(post_draw, 1:3, mean))
  post_quant = reshape2::melt(apply(post_draw, 1:3, quantile, probs = 0.05))
  post_mean$Var1 = rep(xgrid, nrow(post_mean)/length(xgrid))
  post_mean$low = post_quant$value
  post_quant = reshape2::melt(apply(post_draw, 1:3, quantile, probs = 0.95))
  post_mean$upp = post_quant$value
  
  
  ggplot(post_mean) +
    geom_line(aes(x = Var1, y = value, col = factor(Var3)), size = 1.5) +
    geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3)), alpha = 0.4) +
    facet_wrap( ~ Var2, nrow = 2, ncol = 2) +
    labs(x = 'block', y = TeX(tex_string)) +
    scale_color_brewer(name = "", palette = "Set1") +
    scale_fill_brewer(name = "", palette = "Set1") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks=seq(1, K, by = 1)) +
    theme(legend.position = "none")
}


#' Calculate WAIC
#'
#' Function to compute the Watanabe-Akaike information criterion (Gelman, Hwang, 
#' Vehtari, 2014), which estimates the expected out-of-sample-prediction error 
#' using a bias-corrected adjustment of within-sample error.
#'
#' @param model_fit results of a model fit from the lddmm function
#' @return A scalar indicating the WAIC (smaller WAIC denotes better fit)
compute_WAIC = function(model_fit) {
  
  # log pointwise predictive density 
  LPPD = sum(log(apply(exp(model_fit$loglik), 2, mean))) 
  
  # Penalty term
  p_WAIC = sum(apply(model_fit$loglik, 2, var))   
  
  WAIC = -2 * (LPPD - p_WAIC)
  
  return (WAIC) 
}



LDDMM_fix_bound = function(data, hypers, cluster = TRUE, Niter = 5000, burnin = 2000, thin = 5){
  
  # Choose the number of knots (default is between the beginning and the end of 
  # the study, at every block)
  T_min = min(data$block)
  T_max = max(data$block)
  knots = T_min:T_max
  K = length(knots)
  P_smooth = P_smooth1(K + 1)
  
  # Choose a fine grid in order to plot the curves resulting from the spline basis
  # expansion
  xgrid = seq(T_min, T_max, by = 1)
  
  # Extract quantities of interest
  tau = data$r_time
  ind = data$subject
  time = data$block
  cens = data$cens
  D = cbind(data$s, data$d)
  
  n_ind = length(unique(ind))
  ind = as.numeric(plyr::mapvalues(factor(ind),
                                   from = levels(factor(ind)),
                                   to = 1:n_ind))
  B = B_basis(data$block, knots)
  Bgrid = B_basis(xgrid, knots)
  
  samp_size = (Niter - burnin)/thin # sample size
  p = ncol(D) # number of covariates
  d_j = rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] = length(unique(D[!is.na(D[,j]),j]))
  }
  J = ncol(B) # number of locations
  n = nrow(B) # total number of observations
  n_ind = length(unique(ind)) # number of individuals
  
  
  # Rescale the time steps in \{1, ..., T_max\}
  T_max = max(time)
  T_min = min(time)
  time = time - T_min + 1
  T_max = max(time)
  T_min = min(time)
  idx_xy = t(apply(expand.grid(y = 1:d_j[2], x = 1:d_j[1]), 1, rev))
  colnames(idx_xy) = NULL
  
  
  # Set HB structures
  r_HB = 1
  Z_max = if_else(cluster, min(prod(d_j), 6), prod(d_j))
  dim_HB = sum((Z_max - 1)^(0:r_HB) * choose(d_j[2], 0:r_HB))
  
  # Set MCMC objects
  z = array(NA, dim = c(prod(d_j), J, samp_size))
  post_mean_delta = array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_mean_b = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_ind_delta = array(NA, dim = c(d_j[1], n_ind, samp_size))
  post_ind_mu = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  post_ind_b = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  sigma2_mu_us = array(NA, dim = samp_size)
  sigma2_mu_ua = array(NA, dim = samp_size)
  sigma2_b_us = array(NA, dim = samp_size)
  sigma2_b_ua = array(NA, dim = samp_size)
  sigma2_1_mu = array(NA, dim = samp_size)
  sigma2_1_b = array(NA, dim = samp_size)
  pred_time = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  loglik_chain = array(NA, dim = c(samp_size, n))
  
  # Set initial values
  delta_old = array(NA, dim = c(d_j[1], n_ind))
  beta_mu_old = array(NA, dim = c(J, d_j))
  beta_b_old = beta_b_prop = array(NA, dim = c(J, d_j[2]))
  delta_dat = array(NA, dim = n)
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
      if (length(idx_temp) > 0){
        delta_old[s_temp,i_temp] = min(tau[idx_temp])/2
        delta_dat[idx_temp] = delta_old[s_temp,i_temp]
      }
      else{
        delta_old[s_temp,i_temp] = 1E-3
      }
    }
    for(j in 1:d_j[2]){
      beta_mu_old[,s_temp,j] = rep(0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
      beta_b_old[,j] = rep(1.5*log(mean(tau)) - log(sd(tau)), J)
    }
  }
  
  
  low_bound_mu = min(beta_mu_old) - 1.5
  upp_bound_mu = max(beta_mu_old) + 1
  low_bound_b = min(beta_b_old) - 1.5
  upp_bound_b = max(beta_b_old) + 1
  
  beta_mu_star_prop = array(NA, dim = c(J, Z_max))
  beta_mu_u_old = array(0, dim = c(n_ind, J))
  beta_b_u_old = array(0, dim = c(n_ind, J))
  sigma2_1_mu_old = 0.005
  sigma2_1_b_old = 0.005
  sigma2_b_us_old = 0.005
  sigma2_mu_us_old = 0.005
  sigma2_b_ua_old = 0.005
  sigma2_mu_ua_old = 0.005
  
  
  # Message passing structures
  beta_mess = array(NA, dim = c(J, dim_HB))
  z_old = list()
  
  for (j in 1:d_j[1]){
    z_old[[j]] = array(NA, dim = c(d_j[2], J))
    for (jj in 1:d_j[2]){
      if (j == jj){
        z_old[[j]][jj,] = if_else(cluster, j, jj + (j - 1) * d_j[1])
      }
      else{
        z_old[[j]][jj,] = if_else(cluster, sample((d_j[1] + 1):Z_max, 1), jj + (j - 1) * d_j[1])
      }
    }
  }
  z_temp = do.call(rbind, z_old)
  
  
  beta_mu_star_old = array(NA, dim = c(J, Z_max))
  for (i in 1:Z_max){
    beta_mu_star_old[,i] = beta_mu_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                       idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]] 
  }
  
  rand_mat = array(rnorm(prod(d_j)), dim = d_j)
  idx_succ = which(rand_mat == diag(rand_mat))
  
  v_old = array(NA, dim = c(d_j[2], J))
  z_prop = list()
  
  # Transition dynamics objects
  alpha_S_old = alpha_F_old = 1
  Q_S_old = array(NA, dim = c(Z_max, Z_max))
  Q_F_old = array(NA, dim = c(Z_max, Z_max))
  pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                        table_int(z_temp[idx_succ,1], Z_max))
  pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                        table_int(z_temp[-idx_succ,1], Z_max))
  tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
  tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
  for (h in 1:Z_max){
    Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + tr_count_S[h,])
    Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + tr_count_F[h,])
  }
  
  
  # Auxiliary variables
  prob_mat = array(NA, dim = c(dim_HB, dim_HB))
  prob_vec = array(NA, dim = dim_HB)
  
  # MH proposal parameters
  sd_MH_delta = array(0.3, dim = c(d_j[1], n_ind))
  sd_MH_beta_mu = 0.4
  sd_MH_beta_b = array(0.05, dim = c(J, d_j[2]))
  sd_beta_mu_u =  array(0.4, dim = c(n_ind, J))
  sd_beta_b_u =  array(0.4, dim = c(n_ind, J))
  acc_delta = array(0, dim = c(d_j[1], n_ind))
  acc_beta_mu = array(0, dim = c(J, Z_max))
  acc_beta_mu_u = array(0, dim = c(n_ind, J))
  acc_beta_b = array(0, dim = c(J, d_j[2]))
  acc_beta_b_u = array(0, dim = c(n_ind, J))
  n_batch = 0
  
  # Auxiliary variables
  B_beta_mu_dat = array(0, dim = c(n, d_j[1]))
  B_beta_mu_u_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_u_dat = array(0, dim = c(n, d_j[1]))
  mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
  b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
  beta_b_prop = beta_b_old
  
  # Gibbs Sampler
  it = 1
  pb = txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (0) Adaptively tune the MH variance for the proposals of \delta_{s, i}, 
    # beta_u_mu, beta_u_b
    if (iter %% 20 == 0){
      n_batch = n_batch + 1
      delta_n = min(0.01, n_batch^(-0.5))
      
      for (k in 1:J){
        for (d_temp in 1:d_j[2]){	
          if (acc_beta_b[k,d_temp]/iter > 0.44){	
            sd_MH_beta_b[k,d_temp] = exp(log(sd_MH_beta_b[k,d_temp]) + delta_n)	
          }	
          else{	
            sd_MH_beta_b[k,d_temp] = exp(log(sd_MH_beta_b[k,d_temp]) - delta_n)	
          }	
        }	
      }
      
      for (i in 1:n_ind){
        for (x_temp in 1:d_j[1]){
          if (acc_delta[x_temp,i]/iter > 0.44){
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) + delta_n)
          }
          else{
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) - delta_n)
          }
        }
        for (k in 1:J){
          if (acc_beta_mu_u[i,k]/iter > 0.44){
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) + delta_n)
          }
          else{
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) - delta_n)
          }
          if (acc_beta_b_u[i,k]/iter > 0.44){
            sd_beta_b_u[i,k] = exp(log(sd_beta_b_u[i,k]) + delta_n)
          }
          else{
            sd_beta_b_u[i,k] = exp(log(sd_beta_b_u[i,k]) - delta_n)
          }
        }
      }
    }
    
    
    # (1) Update of the delta parameter: \delta_{s,i}: MH with log normal 
    #     proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp = tau[idx_temp]
        cens_temp = cens[idx_temp]
        D_temp = D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        while (delta_prop > min(tau[idx_temp])){
          delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        }
        loglik_prop = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                     b_dat[idx_temp,], 
                                     rep(delta_prop, length(idx_temp)), 
                                     cens_temp, D_temp, TRUE)
        loglik_old = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                    b_dat[idx_temp,], 
                                    rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                    cens_temp, D_temp, TRUE)
        
        alpha_acc = min(0, loglik_prop + log(delta_prop) -
                          loglik_old - log(delta_old[s_temp,i_temp]))
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] = delta_prop
          delta_dat[idx_temp] = delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] = acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    # (2) Update of mu parameters: \mu_{x,y}^{(i)}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time = T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time = T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time = (k-1):k
      }
      
      for (h in 1:Z_max){ # loop over possible latent values
        # tuples (combinations of covariates) that are clustered together via 
        # the latent h
        idx_cov = matrix(idx_xy[which(z_temp[,k] == h),], length(which(z_temp[,k] == h)), p)
        
        X_1k = unique(idx_cov[,1]) # all possible values of x in this cluster
        X_2k = unique(idx_cov[,2]) # all possible values of y in this cluster
        
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{j,k}: posterior update
          # Pick data with covariate levels of x clustered in group h and 
          # at the correct locations
          idx_i = which( (D[,1] %in% X_1k) & (time %in% idx_time) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[which((D[,1] %in% X_1k) & (time %in% idx_time)),]
          
          
          # Normal proposal distribution centered on the current value
          if (k == 1){
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_mu_star_old[k,h], sd_MH_beta_mu)
          }
          # Normal proposal from the prior
          else if (k == J){
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          # Normal proposal from the prior
          else {
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
          
          # Modify the proposed values in the corresponding positions
          for (hh in 1:nrow(idx_cov)){
            B_beta_mu_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] = 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_mu_star_prop[,h]
          }
          
          # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
          mu_prop_dat = exp(B_beta_mu_prop_dat + B_beta_mu_u_dat[idx_i,])
          
          
          if (k == 1){
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          else if (k == J){
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE)
          }
          else {
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          
          alpha_acc = min(0, logpost_prop - logpost_old)
          l_u = log(runif(1))
          
          if (l_u < alpha_acc){
            beta_mu_star_old[k,h] = beta_mu_star_prop[k,h]
            B_beta_mu_dat[idx_i,] = B_beta_mu_prop_dat
            mu_dat[idx_i,] = mu_prop_dat
            acc_beta_mu[k,h] = acc_beta_mu[k,h] + 1
          }
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_mu_star_old[k,h] = runif(1, low_bound_mu, upp_bound_mu)
        }
      }
    }
    
    
    # (2) Update of b parameters: \b_{y}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time = T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time = T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time = (k-1):k
      }
      
      # Pick data with covariate levels of x clustered in group h and 
      # at the correct locations
      
      for (d in 1:d_j[2]){
        idx_i = which( (time %in% idx_time) )
        tau_temp = tau[idx_i]
        cens_temp = cens[idx_i]
        D_temp = D[which( (time %in% idx_time) ),]
        
        # Normal proposal distribution centered on the current value
        if (k == 1){
          beta_b_prop[k,] = beta_b_old[k,]
          beta_b_prop[k,d] = rnorm(1, beta_b_old[k,d], sd_MH_beta_b[k,d])
        }
        # Normal proposal from the prior
        else if (k == J){
          beta_b_prop[k,] = beta_b_old[k,]
          beta_b_prop[k,d] = rnorm(1, beta_b_old[k-1,d], sd_MH_beta_b[k,d])
        }
        # Normal proposal from the prior
        else {
          beta_b_prop[k,] = beta_b_old[k,]
          beta_b_prop[k,d] = rnorm(1, beta_b_old[k-1,d], sd_MH_beta_b[k,d])
        }
        B_beta_b_prop_dat = B_beta_b_dat[idx_i,]
        
        # Modify the proposed values in the corresponding positions
        B_beta_b_prop_dat[,d] = B[idx_i,] %*% beta_b_prop[,d]
        
        
        # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
        b_prop_dat = exp(B_beta_b_prop_dat + B_beta_b_u_dat[idx_i,])
        
        if (k == 1){
          beta_post1 = beta_b_old[k+1,d]
          
          logpost_prop = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                        b_prop_dat, delta_dat[idx_i], 
                                        cens_temp, D_temp, TRUE) -
            0.5/sigma2_1_b_old * sum((beta_b_prop[k,d] - beta_post1)^2)
          logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                       b_dat[idx_i,], delta_dat[idx_i], 
                                       cens_temp, D_temp, TRUE) -
            0.5/sigma2_1_b_old * sum((beta_b_old[k,d] - beta_post1)^2)
        }
        else if (k == J){
          logpost_prop = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                        b_prop_dat, delta_dat[idx_i], 
                                        cens_temp, D_temp, TRUE)
          logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                       b_dat[idx_i,], delta_dat[idx_i], 
                                       cens_temp, D_temp, TRUE)
        }
        else {
          beta_post1 = beta_b_old[k+1,d]
          
          logpost_prop = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                        b_prop_dat, delta_dat[idx_i], 
                                        cens_temp, D_temp, TRUE) -
            0.5/sigma2_1_b_old * sum((beta_b_prop[k,d] - beta_post1)^2)
          logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                       b_dat[idx_i,], delta_dat[idx_i], 
                                       cens_temp, D_temp, TRUE) -
            0.5/sigma2_1_b_old * sum((beta_b_old[k,d] - beta_post1)^2)
        }
        
        alpha_acc = min(0, logpost_prop - logpost_old)
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          beta_b_old[k,d] = beta_b_prop[k,d]
          B_beta_b_dat[idx_i,] = B_beta_b_prop_dat
          b_dat[idx_i,] = b_prop_dat
          acc_beta_b[k,d] = acc_beta_b[k,d] + 1
        }
      }
    }
    
    
    # (3) Update the cluster assignments
    for (x_temp in 1:d_j[1]){ # loop over possible latent values
      if (cluster) {
        beta_mess = array(-Inf, dim = c(J, dim_HB))
        beta_mess[J,] = 1/dim_HB
        
        v_old[,J] = H_ball_unif(z_old[[x_temp]][,J], S = Z_max, r = r_HB)
        z_prop[[J]] = H_ball(v_old[,J], S = Z_max, r = r_HB)
        for (k in (J - 1):1){
          idx_i = which( (time == k) & (D[,1] == x_temp) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[idx_i,]
          
          # (i) Sample the auxiliary variables
          v_temp = H_ball(z_old[[x_temp]][,k], S = Z_max, r = r_HB)
          
          probs = rep(-Inf, dim_HB)
          for (h in 1:dim_HB){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
            
            B_beta_mu_prop_dat = 0.5 * (beta_mu_star_old[k,v_temp[,h]] +
                                          beta_mu_star_old[k+1,z_old[[x_temp]][,k+1]])
            
            mu_dat_prop = exp(t(B_beta_mu_prop_dat + t(B_beta_mu_u_dat[idx_i,])))
            
            probs[h] = g_HB(log_likelihood(tau_temp, mu_dat_prop, b_dat[idx_i,],
                                           delta_dat[idx_i], cens_temp, D_temp, TRUE))
          }
          probs = as.numeric(normalise_log(probs))
          
          v_old[,k] = v_temp[,sample(1:dim_HB, 1, prob = probs)]
          z_prop[[k]] = H_ball(v_old[,k], S = Z_max, r = r_HB)
          
          
          # (ii) Pass messages backwards only in the restricted state space given
          #      by the slice
          z_kp1_temp = which(beta_mess[k+1,] > 0)
          prob_mat = array(-Inf, dim = c(dim_HB, dim_HB))
          
          for (h1 in z_kp1_temp){
            for (h2 in 1:dim_HB){
              
              B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
              B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k,z_prop[[k]][,h2]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k,v_old[,k]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              
              mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_i,])))
              mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_i,])))
              
              prob_mat[h2,h1] = log(beta_mess[k+1,h1]) -
                0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k,z_prop[[k]][,h2]] -
                                               beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])^2) +
                log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_i,],
                               delta_dat[idx_i], cens_temp, D_temp, TRUE) +
                g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_i,],
                                    delta_dat[idx_i], cens_temp, D_temp, TRUE)) +
                sum(log(Q_F_old[cbind(z_prop[[k]][-x_temp,h2],z_prop[[k+1]][-x_temp,h1])])) +
                log(Q_S_old[z_prop[[k]][x_temp,h2],z_prop[[k+1]][x_temp,h1]])
            }
          }
          if ( sum(is.infinite(sum_rows_log(prob_mat))) == dim_HB){
            beta_mess[k,] = 1/dim_HB
          }
          else{
            beta_mess[k,] = as.numeric(sum_rows_log(prob_mat))
            beta_mess[k,] = as.numeric(normalise_log(beta_mess[k,]))
          }
        }
        
        
        # (iii) Sample states forward (only on allowed states)
        idx_fail = (1:d_j[2])[-x_temp]
        # Sample z_1
        prob_vec = log(beta_mess[1,]) + log(pi_S_0[z_prop[[1]][x_temp,]]) +
          colSums(matrix(log(pi_F_0[z_prop[[1]][-x_temp,]]), d_j[2] - 1, dim_HB))
        prob_vec = as.numeric(normalise_log(prob_vec))
        
        idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
        z_old[[x_temp]][,1] = z_prop[[1]][,idx_samp]
        
        # Sample z_k
        for (k in 2:J){
          idx_km1 = which( (time == k - 1) & (D[,1] == x_temp) )
          tau_temp = tau[idx_km1]
          cens_temp = cens[idx_km1]
          D_temp = D[idx_km1,]
          
          prob_vec = log(beta_mess[k,]) + 
            log(Q_S_old[cbind(z_old[[x_temp]][x_temp,k-1], z_prop[[k]][x_temp,])])
          for (kkk in idx_fail){
            prob_vec = prob_vec + log(Q_F_old[cbind(z_old[[x_temp]][kkk,k-1], z_prop[[k]][kkk,])])
          }
          
          for (z_k_temp in which(is.finite(prob_vec))){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_km1,]
            
            B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k-1,v_old[,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            
            mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_km1,])))
            mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_km1,])))
            
            prob_vec[z_k_temp] = prob_vec[z_k_temp] -
              0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                             beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])^2) +
              log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_km1,],
                             delta_dat[idx_km1], cens_temp, D_temp, TRUE) +
              g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_km1,],
                                  delta_dat[idx_km1], cens_temp, D_temp, TRUE))
          }
          prob_vec = as.numeric(normalise_log(prob_vec))
          
          idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
          z_old[[x_temp]][,k] = z_prop[[k]][,idx_samp]
        }
      }
      
      # (4) Assign the cluster specific curves f_{\mu}
      for (y_temp in 1:d_j[2]){
        beta_mu_old[,x_temp,y_temp] = beta_mu_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
      }
      B_beta_mu_dat[(D[,1] == x_temp),] = B[D[,1] == x_temp,] %*% beta_mu_old[,x_temp,]
    }
    z_temp = do.call(rbind, z_old)
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    
    
    
    # (5) Update the transition probabilities
    pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + table_int(z_temp[idx_succ,1], Z_max))
    pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + table_int(z_temp[-idx_succ,1], Z_max))
    tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
    tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
    for (h in 1:Z_max){
      Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                                 tr_count_S[h,])
      Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                                 tr_count_F[h,])
    }
    
    
    alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    while ((alpha_S_prop < 0.01) | (alpha_S_prop > 10)){
      alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    }
    alpha_acc = dgamma(alpha_S_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_S_prop) -
      Z_max^2 * lgamma(alpha_S_prop/Z_max) +
      (alpha_S_prop/Z_max - 1) * sum(log(Q_S_old)) - (
        dgamma(alpha_S_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_S_old) -
          Z_max^2 * lgamma(alpha_S_old/Z_max) +
          (alpha_S_old/Z_max - 1) * sum(log(Q_S_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_S_old = alpha_S_prop
      }
    }
    
    alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    while ((alpha_F_prop < 0.01) | (alpha_F_prop > 10)){
      alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    }
    alpha_acc = dgamma(alpha_F_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_F_prop) -
      Z_max^2 * lgamma(alpha_F_prop/Z_max) +
      (alpha_F_prop/Z_max - 1) * sum(log(Q_F_old)) - (
        dgamma(alpha_F_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_F_old) -
          Z_max^2 * lgamma(alpha_F_old/Z_max) +
          (alpha_F_old/Z_max - 1) * sum(log(Q_F_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_F_old = alpha_F_prop
      }
    }
    
    
    # (6) Correction term for random effects.
    corr_mu = colMeans(beta_mu_u_old)
    corr_b = colMeans(beta_b_u_old)
    beta_mu_u_old = t(t(beta_mu_u_old) - corr_mu)
    beta_b_u_old = t(t(beta_b_u_old) - corr_b)
    
    for (k in 1:J){
      if (k == 1){
        idx_time = which(time == T_min)
      }
      else if (k == J){
        idx_time = which(time == T_max)
      }
      else {
        idx_time = which(time %in% (k-1):k)
      }
      beta_mu_old[k,,] = beta_mu_old[k,,] + corr_mu[k]
      beta_b_old[k,] = beta_b_old[k,] + corr_b[k]
      beta_mu_star_old[k,] = beta_mu_star_old[k,] + corr_mu[k]
      B_beta_mu_dat[idx_time,] = B_beta_mu_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_mu_u_dat[idx_time,] = B_beta_mu_u_dat[idx_time,] - 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_b_dat[idx_time,] = B_beta_b_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_b)
      B_beta_b_u_dat[idx_time,] = B_beta_b_u_dat[idx_time,] -
        as.numeric(B[idx_time,] %*% corr_b)
    }
    
    
    # (7) Update the cluster specific smoothness parameters
    RSS_mu = RSS_b = 0
    for (h2 in 1:d_j[2]){
      RSS_b = RSS_b + as.numeric(crossprod(beta_b_old[,h2], P_smooth) %*% 
                                   beta_b_old[,h2])
      for (h1 in 1:d_j[1]){
        RSS_mu = RSS_mu + as.numeric(crossprod(beta_mu_old[,h1,h2], P_smooth) %*% 
                                       beta_mu_old[,h1,h2])
      }
    }
    
    sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    while(sigma2_1_mu_temp > 0.2){
      sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_mu_temp) - 
                      0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_temp) - 
                      0.5/sigma2_1_mu_temp * RSS_mu + 
                      dhalfcauhy(sqrt(sigma2_1_mu_temp), hypers$s_sigma_mu, T) -
                      (log(sigma2_1_mu_old) - 
                         0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_old) -
                         0.5/sigma2_1_mu_old * RSS_mu +
                         dhalfcauhy(sqrt(sigma2_1_mu_old), hypers$s_sigma_mu, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_mu_old = sigma2_1_mu_temp
    }
    
    sigma2_1_b_temp = exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    while(sigma2_1_b_temp > 0.2){
      sigma2_1_b_temp = exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_b_temp) - 
                      0.5 * (d_j[2]*(J-1)) * log(sigma2_1_b_temp) - 
                      0.5/sigma2_1_b_temp * RSS_b + 
                      dhalfcauhy(sqrt(sigma2_1_b_temp), hypers$s_sigma_b, T) -
                      (log(sigma2_1_b_old) - 
                         0.5 * (d_j[2]*(J-1)) * log(sigma2_1_b_old) -
                         0.5/sigma2_1_b_old * RSS_b +
                         dhalfcauhy(sqrt(sigma2_1_b_old), hypers$s_sigma_b, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_b_old = sigma2_1_b_temp
    }
    
    
    # (8a) Update random effects for b
    reff = sample_reff_b(tau, D, cens, beta_b_u_old, delta_dat, B_beta_b_dat,
                         b_dat, mu_dat, B, P_smooth, ind, time,
                         sigma2_b_us_old, sigma2_b_ua_old, sd_beta_b_u, 
                         acc_beta_b_u)
    beta_b_u_old = reff$beta_u_old
    B_beta_b_u_dat = reff$B_beta_u_dat
    b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
    acc_beta_b_u = reff$acc_beta_u
    
    
    # (8b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_b_ua_old, sigma2_b_us_old,
                               beta_b_u_old, P_smooth, n_ind)
    sigma2_b_ua_old = ls_var$sigma2_ua_old
    sigma2_b_us_old = ls_var$sigma2_us_old
    
    
    # (9a) Update random effects for mu
    reff = sample_reff_mu(tau, D, cens, beta_mu_u_old, delta_dat, b_dat,
                          B_beta_mu_dat, mu_dat, B, P_smooth, ind, time,
                          sigma2_mu_us_old, sigma2_mu_ua_old, sd_beta_mu_u,
                          acc_beta_mu_u)
    beta_mu_u_old = reff$beta_u_old
    B_beta_mu_u_dat = reff$B_beta_u_dat
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    acc_beta_mu_u = reff$acc_beta_u
    
    
    # (9b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_mu_ua_old, sigma2_mu_us_old,
                               beta_mu_u_old, P_smooth, n_ind)
    sigma2_mu_ua_old = ls_var$sigma2_ua_old
    sigma2_mu_us_old = ls_var$sigma2_us_old
    
    
    
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # This is the correction for the random effects integration: we need to 
      # compute the variance of the random effects as detailed in the Supplementary
      # Materials
      cov_reff_mu = solve(diag(J) / sigma2_mu_ua_old + P_smooth / sigma2_mu_us_old)
      corr_term_gr_mu = rep(0, nrow(Bgrid))
      corr_term_mu = rep(0, T_max)
      cov_reff_b = solve(diag(J) / sigma2_b_ua_old + P_smooth / sigma2_b_us_old)
      corr_term_gr_b = rep(0, nrow(Bgrid))
      corr_term_b = rep(0, T_max)
      for (k in 1:J){
        corr_term_gr_mu = corr_term_gr_mu + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_mu[,k]))
        corr_term_mu = corr_term_mu + rowSums(B_basis(1:T_max, knots)[,k] * 
                                                t(t(B_basis(1:T_max, knots)) * cov_reff_mu[,k]))
        corr_term_gr_b = corr_term_gr_b + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_b[,k]))
        corr_term_b = corr_term_b + rowSums(B_basis(1:T_max, knots)[,k] * 
                                              t(t(B_basis(1:T_max, knots)) * cov_reff_b[,k]))
      }
      
      
      
      for (h1 in 1:d_j[1]){
        for (h2 in 1:d_j[2]){
          post_mean_mu[,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + 0.5 * corr_term_gr_mu)
          post_mean_b[,h1,h2,it] = exp(Bgrid %*% beta_b_old[,h2] + 0.5 * corr_term_gr_b)
          
          for (i in 1:n_ind){
            post_ind_mu[,i,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + Bgrid %*% beta_mu_u_old[i,])
            post_ind_b[,i,h1,h2,it] = exp(Bgrid %*% beta_b_old[,h2] + Bgrid %*% beta_b_u_old[i,])
          }
        }
      }
      post_ind_delta[,,it] = delta_old
      
      
      # Sample from the predictive distributions of response category and 
      # response times (so we can use predictive checks to evaluate goodness of fit)
      for (t in 1:T_max){
        for (x_temp in 1:d_j[1]){
          mu_temp = as.numeric(exp(B_basis(t, knots) %*% beta_mu_old[,x_temp,] + 
                                     0.5 * corr_term_mu[t]))
          b_temp = as.numeric(exp(B_basis(t, knots) %*% beta_b_old + 
                                    0.5 * corr_term_b[t]))
          delta_temp = mean(delta_old[x_temp,])
          pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                               b_temp^2)
          pred_ans[t,x_temp,it] = which.min(pred_temp)
          pred_time[t,x_temp,it] = min(pred_temp)
          
          for (i in 1:n_ind){
            mu_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_mu_old[,x_temp,] + 
                                                              beta_mu_u_old[i,])))
            b_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_b_old + 
                                                             beta_b_u_old[i,])))
            delta_temp = delta_old[x_temp,i]
            pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                                 b_temp^2)
            pred_ans_ind[i,t,x_temp,it] = which.min(pred_temp)
            pred_time_ind[i,t,x_temp,it] = min(pred_temp)
          }
        }
      }
      
      
      # Save MCMC objects
      post_mean_delta[it,] = rowMeans(delta_old)
      sigma2_mu_us[it] = sigma2_mu_us_old
      sigma2_mu_ua[it] = sigma2_mu_ua_old
      sigma2_b_us[it] = sigma2_b_us_old
      sigma2_b_ua[it] = sigma2_b_ua_old
      sigma2_1_mu[it] = sigma2_1_mu_old
      sigma2_1_b[it] = sigma2_1_b_old
      loglik_chain[it,] = log_likelihood_ind(tau, mu_dat, b_dat, delta_dat, cens, D)
      # z[,,it] = z_temp		
      
      it = it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list(
    # 'Z' = z, 
    'post_mean_delta' = post_mean_delta, 
    'post_mean_mu' = post_mean_mu,
    'post_mean_b' = post_mean_b,
    'post_ind_delta' = post_ind_delta,
    'post_ind_mu' = post_ind_mu,
    'post_ind_b' = post_ind_b,
    'sigma2_mu_us' = sigma2_mu_us, 
    'sigma2_mu_ua' = sigma2_mu_ua,
    'sigma2_b_us' = sigma2_b_us, 
    'sigma2_b_ua' = sigma2_b_ua,
    'sigma2_1_mu' = sigma2_1_mu, 
    'sigma2_1_b' = sigma2_1_b, 
    'pred_ans' = pred_ans, 
    'pred_time' = pred_time,
    'pred_ans_ind' = pred_ans_ind, 
    'pred_time_ind' = pred_time_ind,
    'loglik' = loglik_chain
  ))
}


LDDMM_const_bound = function(data, hypers, cluster = TRUE, Niter = 5000, burnin = 2000, thin = 5){
  
  # Choose the number of knots (default is between the beginning and the end of 
  # the study, at every block)
  T_min = min(data$block)
  T_max = max(data$block)
  knots = T_min:T_max
  K = length(knots)
  P_smooth = P_smooth1(K + 1)
  
  # Choose a fine grid in order to plot the curves resulting from the spline basis
  # expansion
  xgrid = seq(T_min, T_max, by = 1)
  
  # Extract quantities of interest
  tau = data$r_time
  ind = data$subject
  time = data$block
  cens = data$cens
  D = cbind(data$s, data$d)
  
  n_ind = length(unique(ind))
  ind = as.numeric(plyr::mapvalues(factor(ind),
                                   from = levels(factor(ind)),
                                   to = 1:n_ind))
  B = B_basis(data$block, knots)
  Bgrid = B_basis(xgrid, knots)
  
  samp_size = (Niter - burnin)/thin # sample size
  p = ncol(D) # number of covariates
  d_j = rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] = length(unique(D[!is.na(D[,j]),j]))
  }
  J = ncol(B) # number of locations
  n = nrow(B) # total number of observations
  n_ind = length(unique(ind)) # number of individuals
  
  
  # Rescale the time steps in \{1, ..., T_max\}
  T_max = max(time)
  T_min = min(time)
  time = time - T_min + 1
  T_max = max(time)
  T_min = min(time)
  idx_xy = t(apply(expand.grid(y = 1:d_j[2], x = 1:d_j[1]), 1, rev))
  colnames(idx_xy) = NULL
  
  
  # Set HB structures
  r_HB = 1
  Z_max = if_else(cluster, min(prod(d_j), 6), prod(d_j))
  dim_HB = sum((Z_max - 1)^(0:r_HB) * choose(d_j[2], 0:r_HB))
  
  # Set MCMC objects
  z = array(NA, dim = c(prod(d_j), J, samp_size))
  post_mean_delta = array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_mean_b = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_ind_delta = array(NA, dim = c(d_j[1], n_ind, samp_size))
  post_ind_mu = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  post_ind_b = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  sigma2_mu_us = array(NA, dim = samp_size)
  sigma2_mu_ua = array(NA, dim = samp_size)
  sigma2_b_us = array(NA, dim = samp_size)
  sigma2_b_ua = array(NA, dim = samp_size)
  sigma2_1_mu = array(NA, dim = samp_size)
  sigma2_1_b = array(NA, dim = samp_size)
  pred_time = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  loglik_chain = array(NA, dim = c(samp_size, n))
  
  # Set initial values
  delta_old = array(NA, dim = c(d_j[1], n_ind))
  beta_mu_old = array(NA, dim = c(J, d_j))
  delta_dat = array(NA, dim = n)
  b_old = array(NA, dim = d_j)
  B_beta_b_dat = array(0, dim = c(n, d_j[1]))
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
      if (length(idx_temp) > 0){
        delta_old[s_temp,i_temp] = min(tau[idx_temp])/2
        delta_dat[idx_temp] = delta_old[s_temp,i_temp]
      }
      else{
        delta_old[s_temp,i_temp] = 1E-3
      }
    }
    for(j in 1:d_j[2]){
      beta_mu_old[,s_temp,j] = rep(0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
      b_old[s_temp,j] = 0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp]))
      B_beta_b_dat[which((D[,1] == s_temp)),j] = b_old[s_temp,j]
    }
  }
  # b_old = 1.5*log(mean(tau)) - log(sd(tau))
  # B_beta_b_dat = array(b_old, dim = c(n, d_j[1]))
  
  
  low_bound_mu = min(beta_mu_old) - 1.5
  upp_bound_mu = max(beta_mu_old) + 1
  
  beta_mu_star_prop = array(NA, dim = c(J, Z_max))
  beta_mu_u_old = array(0, dim = c(n_ind, J))
  sigma2_1_mu_old = 0.005
  sigma2_mu_us_old = 0.005
  sigma2_mu_ua_old = 0.005
  
  
  # Message passing structures
  beta_mess = array(NA, dim = c(J, dim_HB))
  z_old = list()
  
  for (j in 1:d_j[1]){
    z_old[[j]] = array(NA, dim = c(d_j[2], J))
    for (jj in 1:d_j[2]){
      if (j == jj){
        z_old[[j]][jj,] = if_else(cluster, j, jj + (j - 1) * d_j[1])
      }
      else{
        z_old[[j]][jj,] = if_else(cluster, sample((d_j[1] + 1):Z_max, 1), jj + (j - 1) * d_j[1])
      }
    }
  }
  z_temp = do.call(rbind, z_old)
  
  
  beta_mu_star_old = array(NA, dim = c(J, Z_max))
  for (i in 1:Z_max){
    beta_mu_star_old[,i] = beta_mu_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                       idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]]
  }
  
  rand_mat = array(rnorm(prod(d_j)), dim = d_j)
  idx_succ = which(rand_mat == diag(rand_mat))
  
  v_old = array(NA, dim = c(d_j[2], J))
  z_prop = list()
  sigma2_b_ua_old = 0.005
  # Transition dynamics objects
  alpha_S_old = alpha_F_old = 1
  Q_S_old = array(NA, dim = c(Z_max, Z_max))
  Q_F_old = array(NA, dim = c(Z_max, Z_max))
  pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                        table_int(z_temp[idx_succ,1], Z_max))
  pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                        table_int(z_temp[-idx_succ,1], Z_max))
  tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
  tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
  for (h in 1:Z_max){
    Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + tr_count_S[h,])
    Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + tr_count_F[h,])
  }
  
  
  # Auxiliary variables
  prob_mat = array(NA, dim = c(dim_HB, dim_HB))
  prob_vec = array(NA, dim = dim_HB)
  beta_b_u_chain = array(NA, dim = c(samp_size, n_ind))
  
  # MH proposal parameters
  sd_MH_delta = array(0.3, dim = c(d_j[1], n_ind))
  sd_MH_beta_mu = 0.4
  sd_MH_beta_b = array(0.05, dim = d_j)
  acc_b = array(0, dim = d_j)
  sd_beta_mu_u =  array(0.4, dim = c(n_ind, J))
  sd_beta_b_u =  array(0.4, dim = n_ind)
  acc_delta = array(0, dim = c(d_j[1], n_ind))
  acc_beta_mu_u = array(0, dim = c(n_ind, J))
  acc_beta_b_u = array(0, dim = n_ind)
  n_batch = 0
  
  # Auxiliary variables
  B_beta_mu_dat = array(0, dim = c(n, d_j[1]))
  B_beta_mu_u_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_u_dat = array(0, dim = c(n, d_j[1]))
  mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
  b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
  beta_b_u_old = array(0, dim = n_ind)
  
  
  # Gibbs Sampler
  it = 1
  pb = txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (0) Adaptively tune the MH variance for the proposals of \delta_{s, i}, 
    # beta_u_mu, beta_u_b
    if (iter %% 20 == 0){
      n_batch = n_batch + 1
      delta_n = min(0.01, n_batch^(-0.5))
      
      for (i in 1:n_ind){
        if (acc_beta_b_u[i]/iter > 0.44){
          sd_beta_b_u[i] = exp(log(sd_beta_b_u[i]) + delta_n)
        }
        else{
          sd_beta_b_u[i] = exp(log(sd_beta_b_u[i]) - delta_n)
        }
        
        for (x_temp in 1:d_j[1]){
          if (acc_delta[x_temp,i]/iter > 0.44){
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) + delta_n)
          }
          else{
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) - delta_n)
          }
        }
        for (k in 1:J){
          if (acc_beta_mu_u[i,k]/iter > 0.44){
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) + delta_n)
          }
          else{
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) - delta_n)
          }
        }
      }
      for (d1 in 1:d_j[1]){
        for (d2 in 1:d_j[2]){
          if (acc_b[d1,d2]/iter > 0.44){
            sd_MH_beta_b[d1,d2] = exp(log(sd_MH_beta_b[d1,d2]) + delta_n)
          }
          else{
            sd_MH_beta_b[d1,d2] = exp(log(sd_MH_beta_b[d1,d2]) - delta_n)
          }
          
        }
      }
    }
    
    
    # (1) Update of the delta parameter: \delta_{s,i}: MH with log normal 
    #     proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp = tau[idx_temp]
        cens_temp = cens[idx_temp]
        D_temp = D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        while (delta_prop > min(tau[idx_temp])){
          delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        }
        loglik_prop = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                     b_dat[idx_temp,], 
                                     rep(delta_prop, length(idx_temp)), 
                                     cens_temp, D_temp, TRUE)
        loglik_old = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                    b_dat[idx_temp,], 
                                    rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                    cens_temp, D_temp, TRUE)
        
        alpha_acc = min(0, loglik_prop + log(delta_prop) -
                          loglik_old - log(delta_old[s_temp,i_temp]))
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] = delta_prop
          delta_dat[idx_temp] = delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] = acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    # (2) Update of mu parameters: \mu_{x,y}^{(i)}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time = T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time = T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time = (k-1):k
      }
      
      for (h in 1:Z_max){ # loop over possible latent values
        # tuples (combinations of covariates) that are clustered together via 
        # the latent h
        idx_cov = matrix(idx_xy[which(z_temp[,k] == h),], length(which(z_temp[,k] == h)), p)
        
        X_1k = unique(idx_cov[,1]) # all possible values of x in this cluster
        X_2k = unique(idx_cov[,2]) # all possible values of y in this cluster
        
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{j,k}: posterior update
          # Pick data with covariate levels of x clustered in group h and 
          # at the correct locations
          idx_i = which( (D[,1] %in% X_1k) & (time %in% idx_time) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[which((D[,1] %in% X_1k) & (time %in% idx_time)),]
          
          
          # Normal proposal distribution centered on the current value
          if (k == 1){
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_mu_star_old[k,h], sd_MH_beta_mu)
          }
          # Normal proposal from the prior
          else if (k == J){
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          # Normal proposal from the prior
          else {
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
          
          # Modify the proposed values in the corresponding positions
          for (hh in 1:nrow(idx_cov)){
            B_beta_mu_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] = 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_mu_star_prop[,h]
          }
          
          # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
          mu_prop_dat = exp(B_beta_mu_prop_dat + B_beta_mu_u_dat[idx_i,])
          
          
          if (k == 1){
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          else if (k == J){
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE)
          }
          else {
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          
          alpha_acc = min(0, logpost_prop - logpost_old)
          l_u = log(runif(1))
          
          if (l_u < alpha_acc){
            
            beta_mu_star_old[k,h] = beta_mu_star_prop[k,h]
            B_beta_mu_dat[idx_i,] = B_beta_mu_prop_dat
            mu_dat[idx_i,] = mu_prop_dat
          }
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_mu_star_old[k,h] = runif(1, low_bound_mu, upp_bound_mu)
        }
      }
    }
    
    # 2(a) Update the constant boundary parameters
    for (d1 in 1:d_j[1]){
      idx_i = which( (D[,1] == d1) )
      tau_temp = tau[idx_i]
      cens_temp = cens[idx_i]
      D_temp = D[which( (D[,1] == d1) ),]
      
      for (d2 in 1:d_j[2]){
        b_prop = rnorm(1, b_old[d1,d2], sd_MH_beta_b[d1,d2])
        
        # Modify the proposed values in the corresponding positions
        B_beta_b_prop_dat = B_beta_b_dat[idx_i,]
        B_beta_b_prop_dat[,d2] = b_prop
        
        # This is the proposed value for b
        b_prop_dat = exp(B_beta_b_prop_dat)
        
        logpost_prop = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                      b_prop_dat, delta_dat[idx_i], 
                                      cens_temp, D_temp, TRUE)
        logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                     b_dat[idx_i,], delta_dat[idx_i], 
                                     cens_temp, D_temp, TRUE)
        
        alpha_acc = min(0, logpost_prop - logpost_old)
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          b_old[d1,d2] = b_prop
          # b_dat = b_prop_dat
          B_beta_b_dat[idx_i,] = B_beta_b_prop_dat
          b_dat[idx_i,] = b_prop_dat
          acc_b[d1,d2] = acc_b[d1,d2] + 1
        }
      }
    }
    
    
    # (3) Update the cluster assignments
    for (x_temp in 1:d_j[1]){ # loop over possible latent values
      if (cluster) {
        
        beta_mess = array(-Inf, dim = c(J, dim_HB))
        beta_mess[J,] = 1/dim_HB
        
        v_old[,J] = H_ball_unif(z_old[[x_temp]][,J], S = Z_max, r = r_HB)
        z_prop[[J]] = H_ball(v_old[,J], S = Z_max, r = r_HB)
        for (k in (J - 1):1){
          idx_i = which( (time == k) & (D[,1] == x_temp) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[idx_i,]
          
          # (i) Sample the auxiliary variables
          v_temp = H_ball(z_old[[x_temp]][,k], S = Z_max, r = r_HB)
          
          probs = rep(-Inf, dim_HB)
          for (h in 1:dim_HB){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
            
            B_beta_mu_prop_dat = 0.5 * (beta_mu_star_old[k,v_temp[,h]] +
                                          beta_mu_star_old[k+1,z_old[[x_temp]][,k+1]])
            
            mu_dat_prop = exp(t(B_beta_mu_prop_dat + t(B_beta_mu_u_dat[idx_i,])))
            
            probs[h] = g_HB(log_likelihood(tau_temp, mu_dat_prop, b_dat[idx_i,],
                                           delta_dat[idx_i], cens_temp, D_temp, TRUE))
          }
          probs = as.numeric(normalise_log(probs))
          
          v_old[,k] = v_temp[,sample(1:dim_HB, 1, prob = probs)]
          z_prop[[k]] = H_ball(v_old[,k], S = Z_max, r = r_HB)
          
          
          # (ii) Pass messages backwards only in the restricted state space given
          #      by the slice
          z_kp1_temp = which(beta_mess[k+1,] > 0)
          prob_mat = array(-Inf, dim = c(dim_HB, dim_HB))
          
          for (h1 in z_kp1_temp){
            for (h2 in 1:dim_HB){
              
              B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
              B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k,z_prop[[k]][,h2]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k,v_old[,k]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              
              mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_i,])))
              mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_i,])))
              
              prob_mat[h2,h1] = log(beta_mess[k+1,h1]) -
                0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k,z_prop[[k]][,h2]] -
                                               beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])^2) +
                log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_i,],
                               delta_dat[idx_i], cens_temp, D_temp, TRUE) +
                g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_i,],
                                    delta_dat[idx_i], cens_temp, D_temp, TRUE)) +
                sum(log(Q_F_old[cbind(z_prop[[k]][-x_temp,h2],z_prop[[k+1]][-x_temp,h1])])) +
                log(Q_S_old[z_prop[[k]][x_temp,h2],z_prop[[k+1]][x_temp,h1]])
            }
          }
          if ( sum(is.infinite(sum_rows_log(prob_mat))) == dim_HB){
            beta_mess[k,] = 1/dim_HB
          }
          else{
            beta_mess[k,] = as.numeric(sum_rows_log(prob_mat))
            beta_mess[k,] = as.numeric(normalise_log(beta_mess[k,]))
          }
        }
        
        
        # (iii) Sample states forward (only on allowed states)
        idx_fail = (1:d_j[2])[-x_temp]
        # Sample z_1
        prob_vec = log(beta_mess[1,]) + log(pi_S_0[z_prop[[1]][x_temp,]]) +
          colSums(matrix(log(pi_F_0[z_prop[[1]][-x_temp,]]), d_j[2] - 1, dim_HB))
        prob_vec = as.numeric(normalise_log(prob_vec))
        
        idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
        z_old[[x_temp]][,1] = z_prop[[1]][,idx_samp]
        
        # Sample z_k
        for (k in 2:J){
          idx_km1 = which( (time == k - 1) & (D[,1] == x_temp) )
          tau_temp = tau[idx_km1]
          cens_temp = cens[idx_km1]
          D_temp = D[idx_km1,]
          
          prob_vec = log(beta_mess[k,]) + 
            log(Q_S_old[cbind(z_old[[x_temp]][x_temp,k-1], z_prop[[k]][x_temp,])])
          for (kkk in idx_fail){
            prob_vec = prob_vec + log(Q_F_old[cbind(z_old[[x_temp]][kkk,k-1], z_prop[[k]][kkk,])])
          }
          
          for (z_k_temp in which(is.finite(prob_vec))){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_km1,]
            
            B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k-1,v_old[,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            
            mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_km1,])))
            mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_km1,])))
            
            prob_vec[z_k_temp] = prob_vec[z_k_temp] -
              0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                             beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])^2) +
              log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_km1,],
                             delta_dat[idx_km1], cens_temp, D_temp, TRUE) +
              g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_km1,],
                                  delta_dat[idx_km1], cens_temp, D_temp, TRUE))
          }
          prob_vec = as.numeric(normalise_log(prob_vec))
          
          idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
          z_old[[x_temp]][,k] = z_prop[[k]][,idx_samp]
        }
      }
      
      # (4) Assign the cluster specific curves f_{\mu}
      for (y_temp in 1:d_j[2]){
        beta_mu_old[,x_temp,y_temp] = beta_mu_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
      }
      B_beta_mu_dat[(D[,1] == x_temp),] = B[D[,1] == x_temp,] %*% beta_mu_old[,x_temp,]
    }
    z_temp = do.call(rbind, z_old)
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    
    
    # (5) Update the transition probabilities
    pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + table_int(z_temp[idx_succ,1], Z_max))
    pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + table_int(z_temp[-idx_succ,1], Z_max))
    tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
    tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
    for (h in 1:Z_max){
      Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                                 tr_count_S[h,])
      Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                                 tr_count_F[h,])
    }
    
    
    alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    while ((alpha_S_prop < 0.01) | (alpha_S_prop > 10)){
      alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    }
    alpha_acc = dgamma(alpha_S_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_S_prop) -
      Z_max^2 * lgamma(alpha_S_prop/Z_max) +
      (alpha_S_prop/Z_max - 1) * sum(log(Q_S_old)) - (
        dgamma(alpha_S_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_S_old) -
          Z_max^2 * lgamma(alpha_S_old/Z_max) +
          (alpha_S_old/Z_max - 1) * sum(log(Q_S_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_S_old = alpha_S_prop
      }
    }
    
    alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    while ((alpha_F_prop < 0.01) | (alpha_F_prop > 10)){
      alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    }
    alpha_acc = dgamma(alpha_F_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_F_prop) -
      Z_max^2 * lgamma(alpha_F_prop/Z_max) +
      (alpha_F_prop/Z_max - 1) * sum(log(Q_F_old)) - (
        dgamma(alpha_F_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_F_old) -
          Z_max^2 * lgamma(alpha_F_old/Z_max) +
          (alpha_F_old/Z_max - 1) * sum(log(Q_F_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_F_old = alpha_F_prop
      }
    }
    
    
    # (6) Correction term for random effects.
    corr_mu = colMeans(beta_mu_u_old)
    beta_mu_u_old = t(t(beta_mu_u_old) - corr_mu)
    
    for (k in 1:J){
      if (k == 1){
        idx_time = which(time == T_min)
      }
      else if (k == J){
        idx_time = which(time == T_max)
      }
      else {
        idx_time = which(time %in% (k-1):k)
      }
      beta_mu_old[k,,] = beta_mu_old[k,,] + corr_mu[k]
      beta_mu_star_old[k,] = beta_mu_star_old[k,] + corr_mu[k]
      B_beta_mu_dat[idx_time,] = B_beta_mu_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_mu_u_dat[idx_time,] = B_beta_mu_u_dat[idx_time,] - 
        as.numeric(B[idx_time,] %*% corr_mu)
    }
    
    corr_b = mean(beta_b_u_old)
    beta_b_u_old = beta_b_u_old - corr_b
    
    b_old = b_old + corr_b
    B_beta_b_dat = B_beta_b_dat + corr_b
    
    
    # (7) Update the cluster specific smoothness parameters
    RSS_mu = RSS_b = 0
    for (h1 in 1:d_j[1]){
      for (h2 in 1:d_j[2]){
        RSS_mu = RSS_mu + as.numeric(crossprod(beta_mu_old[,h1,h2], P_smooth) %*% 
                                       beta_mu_old[,h1,h2])
      }
    }
    
    sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    while(sigma2_1_mu_temp > 0.2){
      sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_mu_temp) - 
                      0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_temp) - 
                      0.5/sigma2_1_mu_temp * RSS_mu + 
                      dhalfcauhy(sqrt(sigma2_1_mu_temp), hypers$s_sigma_mu, T) -
                      (log(sigma2_1_mu_old) - 
                         0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_old) -
                         0.5/sigma2_1_mu_old * RSS_mu +
                         dhalfcauhy(sqrt(sigma2_1_mu_old), hypers$s_sigma_mu, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_mu_old = sigma2_1_mu_temp
    }
    
    
    # (9a) Update random effects for mu
    reff = sample_reff_mu(tau, D, cens, beta_mu_u_old, delta_dat, b_dat,
                          B_beta_mu_dat, mu_dat, B, P_smooth, ind, time,
                          sigma2_mu_us_old, sigma2_mu_ua_old, sd_beta_mu_u,
                          acc_beta_mu_u)
    beta_mu_u_old = reff$beta_u_old
    B_beta_mu_u_dat = reff$B_beta_u_dat
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    acc_beta_mu_u = reff$acc_beta_u
    
    
    # (9b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_mu_ua_old, sigma2_mu_us_old,
                               beta_mu_u_old, P_smooth, n_ind)
    sigma2_mu_ua_old = ls_var$sigma2_ua_old
    sigma2_mu_us_old = ls_var$sigma2_us_old
    
    
    for (i in 1:n_ind) {
      beta_b_u_prop = rnorm(1, beta_b_u_old[i], sd_beta_b_u[i])
      # browser()
      idx_i = which(ind == i)
      
      logpr_prop = - 0.5 * ( beta_b_u_prop^2 / sigma2_b_ua_old)
      logpr_old = - 0.5 * ( beta_b_u_old[i]^2 / sigma2_b_ua_old)
      
      b_dat_prop = exp(B_beta_b_dat[idx_i,] + beta_b_u_prop)
      
      
      logpost_prop = log_likelihood(tau[idx_i],
                                    mu_dat[idx_i,],
                                    b_dat_prop,
                                    delta_dat[idx_i],
                                    cens[idx_i],
                                    D[idx_i,], TRUE) +
        logpr_prop
      logpost_old = log_likelihood(tau[idx_i],
                                   mu_dat[idx_i,],
                                   b_dat[idx_i,],
                                   delta_dat[idx_i],
                                   cens[idx_i],
                                   D[idx_i,], TRUE) +
        logpr_old
      
      
      alpha_acc = logpost_prop - logpost_old;
      l_u = log(runif(1))
      if (l_u < alpha_acc){
        beta_b_u_old[i] = beta_b_u_prop
        b_dat[idx_i,] = b_dat_prop
        acc_beta_b_u[i] = acc_beta_b_u[i] + 1
        B_beta_b_u_dat[idx_i,] = beta_b_u_old[i]
      }
    }
    
    sigma2_b_ua_prop = exp(rnorm(1, log(sigma2_b_ua_old), sd = 0.2))
    lu = log(runif(1))
    
    log_prop = - 0.5 * n_ind * log(sigma2_b_ua_prop) -
      0.5 * sum(beta_b_u_old^2) / sigma2_b_ua_prop -
      log(1 + sigma2_b_ua_prop^2)
    log_old = - 0.5 * n_ind * log(sigma2_b_ua_old) -
      0.5 * sum(beta_b_u_old^2) / sigma2_b_ua_old -
      log(1 + sigma2_b_ua_old^2)
    
    alpha = min(c(0, log_prop + log(sigma2_b_ua_prop) -
                    log_old - log(sigma2_b_ua_old)))
    if (lu < alpha){
      sigma2_b_ua_old = sigma2_b_ua_prop
    }
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # This is the correction for the random effects integration: we need to 
      # compute the variance of the random effects as detailed in the Supplementary
      # Materials
      cov_reff_mu = solve(diag(J) / sigma2_mu_ua_old + P_smooth / sigma2_mu_us_old)
      corr_term_gr_mu = rep(0, nrow(Bgrid))
      corr_term_mu = rep(0, T_max)
      for (k in 1:J){
        corr_term_gr_mu = corr_term_gr_mu + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_mu[,k]))
        corr_term_mu = corr_term_mu + rowSums(B_basis(1:T_max, knots)[,k] * 
                                                t(t(B_basis(1:T_max, knots)) * cov_reff_mu[,k]))
      }
      
      
      
      for (h1 in 1:d_j[1]){
        for (h2 in 1:d_j[2]){
          post_mean_mu[,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + 0.5 * corr_term_gr_mu)
          post_mean_b[,h1,h2,it] = exp(b_old[h1,h2] + 0.5 * sigma2_b_ua_old)
          
          for (i in 1:n_ind){
            post_ind_mu[,i,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + Bgrid %*% beta_mu_u_old[i,])
            post_ind_b[,i,h1,h2,it] = exp(b_old[h1,h2] + beta_b_u_old[i])
          }
        }
      }
      post_ind_delta[,,it] = delta_old
      
      
      # Sample from the predictive distributions of response category and 
      # response times (so we can use predictive checks to evaluate goodness of fit)
      for (t in 1:T_max){
        for (x_temp in 1:d_j[1]){
          mu_temp = as.numeric(exp(B_basis(t, knots) %*% beta_mu_old[,x_temp,] + 
                                     0.5 * corr_term_mu[t]))
          b_temp = exp(b_old[x_temp,] + 0.5 * sigma2_b_ua_old)
          delta_temp = mean(delta_old[x_temp,])
          pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                               b_temp^2)
          pred_ans[t,x_temp,it] = which.min(pred_temp)
          pred_time[t,x_temp,it] = min(pred_temp)
          
          for (i in 1:n_ind){
            mu_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_mu_old[,x_temp,] + 
                                                              beta_mu_u_old[i,])))
            b_temp = exp(b_old[x_temp,] + beta_b_u_old[i])
            delta_temp = delta_old[x_temp,i]
            pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                                 b_temp^2)
            pred_ans_ind[i,t,x_temp,it] = which.min(pred_temp)
            pred_time_ind[i,t,x_temp,it] = min(pred_temp)
          }
        }
      }
      
      
      # Save MCMC objects
      post_mean_delta[it,] = rowMeans(delta_old)
      sigma2_mu_us[it] = sigma2_mu_us_old
      sigma2_mu_ua[it] = sigma2_mu_ua_old
      sigma2_1_mu[it] = sigma2_1_mu_old
      sigma2_b_ua[it] = sigma2_b_ua_old
      # z[,,it] = z_temp		
      loglik_chain[it,] = log_likelihood_ind(tau, mu_dat, b_dat, delta_dat, cens, D)
      
      it = it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list(
    # 'Z' = z, 
    'post_mean_delta' = post_mean_delta, 
    'post_mean_mu' = post_mean_mu,
    'post_mean_b' = post_mean_b,
    'post_ind_delta' = post_ind_delta,
    'post_ind_mu' = post_ind_mu,
    'post_ind_b' = post_ind_b,
    'sigma2_mu_us' = sigma2_mu_us, 
    'sigma2_mu_ua' = sigma2_mu_ua,
    'sigma2_b_ua' = sigma2_b_ua,
    'sigma2_1_mu' = sigma2_1_mu, 
    'pred_ans' = pred_ans, 
    'pred_time' = pred_time,
    'pred_ans_ind' = pred_ans_ind, 
    'pred_time_ind' = pred_time_ind, 
    'loglik' = loglik_chain
  ))
}




LDDMM_fixandconst_bound = function(data, hypers, cluster = TRUE, Niter = 5000, burnin = 2000, thin = 5){
  
  # Choose the number of knots (default is between the beginning and the end of 
  # the study, at every block)
  T_min = min(data$block)
  T_max = max(data$block)
  knots = T_min:T_max
  K = length(knots)
  P_smooth = P_smooth1(K + 1)
  
  # Choose a fine grid in order to plot the curves resulting from the spline basis
  # expansion
  xgrid = seq(T_min, T_max, by = 1)
  
  # Extract quantities of interest
  tau = data$r_time
  ind = data$subject
  time = data$block
  cens = data$cens
  D = cbind(data$s, data$d)
  
  n_ind = length(unique(ind))
  ind = as.numeric(plyr::mapvalues(factor(ind),
                                   from = levels(factor(ind)),
                                   to = 1:n_ind))
  B = B_basis(data$block, knots)
  Bgrid = B_basis(xgrid, knots)
  
  samp_size = (Niter - burnin)/thin # sample size
  p = ncol(D) # number of covariates
  d_j = rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] = length(unique(D[!is.na(D[,j]),j]))
  }
  J = ncol(B) # number of locations
  n = nrow(B) # total number of observations
  n_ind = length(unique(ind)) # number of individuals
  
  
  # Rescale the time steps in \{1, ..., T_max\}
  T_max = max(time)
  T_min = min(time)
  time = time - T_min + 1
  T_max = max(time)
  T_min = min(time)
  idx_xy = t(apply(expand.grid(y = 1:d_j[2], x = 1:d_j[1]), 1, rev))
  colnames(idx_xy) = NULL
  beta_b_u_old = array(0, dim = n_ind)
  
  
  # Set HB structures
  r_HB = 1
  Z_max = if_else(cluster, min(prod(d_j), 6), prod(d_j))
  dim_HB = sum((Z_max - 1)^(0:r_HB) * choose(d_j[2], 0:r_HB))
  
  # Set MCMC objects
  z = array(NA, dim = c(prod(d_j), J, samp_size))
  post_mean_delta = array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_mean_b = array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_ind_delta = array(NA, dim = c(d_j[1], n_ind, samp_size))
  post_ind_mu = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  post_ind_b = array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  sigma2_mu_us = array(NA, dim = samp_size)
  sigma2_mu_ua = array(NA, dim = samp_size)
  sigma2_b_ua = array(NA, dim = samp_size)
  sigma2_1_mu = array(NA, dim = samp_size)
  sigma2_1_b = array(NA, dim = samp_size)
  pred_time = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans = array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind = array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  loglik_chain = array(NA, dim = c(samp_size, n))
  
  # Set initial values
  delta_old = array(NA, dim = c(d_j[1], n_ind))
  beta_mu_old = array(NA, dim = c(J, d_j))
  delta_dat = array(NA, dim = n)
  b_old = array(NA, dim = d_j[2])
  B_beta_b_dat = array(0, dim = c(n, d_j[1]))
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
      if (length(idx_temp) > 0){
        delta_old[s_temp,i_temp] = min(tau[idx_temp])/2
        delta_dat[idx_temp] = delta_old[s_temp,i_temp]
      }
      else{
        delta_old[s_temp,i_temp] = 1E-3
      }
    }
    for(j in 1:d_j[2]){
      beta_mu_old[,s_temp,j] = rep(0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
      b_old[j] = 0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp]))
      B_beta_b_dat[which((D[,1] == s_temp)),j] = b_old[j]
    }
  }
  # b_old = 1.5*log(mean(tau)) - log(sd(tau))
  # B_beta_b_dat = array(b_old, dim = c(n, d_j[1]))
  
  
  low_bound_mu = min(beta_mu_old) - 1.5
  upp_bound_mu = max(beta_mu_old) + 1
  
  beta_mu_star_prop = array(NA, dim = c(J, Z_max))
  beta_mu_u_old = array(0, dim = c(n_ind, J))
  sigma2_1_mu_old = 0.005
  sigma2_mu_us_old = 0.005
  sigma2_mu_ua_old = 0.005
  sigma2_b_ua_old = 0.005
  
  
  # Message passing structures
  beta_mess = array(NA, dim = c(J, dim_HB))
  z_old = list()
  
  for (j in 1:d_j[1]){
    z_old[[j]] = array(NA, dim = c(d_j[2], J))
    for (jj in 1:d_j[2]){
      if (j == jj){
        z_old[[j]][jj,] = if_else(cluster, j, jj + (j - 1) * d_j[1])
      }
      else{
        z_old[[j]][jj,] = if_else(cluster, sample((d_j[1] + 1):Z_max, 1), jj + (j - 1) * d_j[1])
      }
    }
  }
  z_temp = do.call(rbind, z_old)
  
  
  beta_mu_star_old = array(NA, dim = c(J, Z_max))
  for (i in 1:Z_max){
    beta_mu_star_old[,i] = beta_mu_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                       idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]]
  }
  
  rand_mat = array(rnorm(prod(d_j)), dim = d_j)
  idx_succ = which(rand_mat == diag(rand_mat))
  
  v_old = array(NA, dim = c(d_j[2], J))
  z_prop = list()
  
  # Transition dynamics objects
  alpha_S_old = alpha_F_old = 1
  Q_S_old = array(NA, dim = c(Z_max, Z_max))
  Q_F_old = array(NA, dim = c(Z_max, Z_max))
  pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                        table_int(z_temp[idx_succ,1], Z_max))
  pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                        table_int(z_temp[-idx_succ,1], Z_max))
  tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
  tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
  for (h in 1:Z_max){
    Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + tr_count_S[h,])
    Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + tr_count_F[h,])
  }
  
  
  # Auxiliary variables
  prob_mat = array(NA, dim = c(dim_HB, dim_HB))
  prob_vec = array(NA, dim = dim_HB)
  beta_b_u_chain = array(NA, dim = c(samp_size, n_ind))
  
  # MH proposal parameters
  sd_MH_delta = array(0.3, dim = c(d_j[1], n_ind))
  sd_MH_beta_mu = 0.4
  sd_MH_beta_b = array(0.05, dim = d_j[2])
  acc_b = array(0, dim = d_j[2])
  sd_beta_mu_u =  array(0.4, dim = c(n_ind, J))
  sd_beta_b_u =  array(0.4, dim = n_ind)
  acc_delta = array(0, dim = c(d_j[1], n_ind))
  acc_beta_mu_u = array(0, dim = c(n_ind, J))
  acc_beta_b_u = array(0, dim = n_ind)
  n_batch = 0
  
  # Auxiliary variables
  B_beta_mu_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_dat = array(b_old, dim = c(n, d_j[1]))
  B_beta_mu_u_dat = array(0, dim = c(n, d_j[1]))
  B_beta_b_u_dat = array(0, dim = c(n, d_j[1]))
  mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
  b_dat = exp(B_beta_b_dat + B_beta_b_u_dat)
  
  
  # Gibbs Sampler
  it = 1
  pb = txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (0) Adaptively tune the MH variance for the proposals of \delta_{s, i}, 
    # beta_u_mu, beta_u_b
    if (iter %% 20 == 0){
      n_batch = n_batch + 1
      delta_n = min(0.01, n_batch^(-0.5))
      
      for (i in 1:n_ind){
        if (acc_beta_b_u[i]/iter > 0.44){
          sd_beta_b_u[i] = exp(log(sd_beta_b_u[i]) + delta_n)
        }
        else{
          sd_beta_b_u[i] = exp(log(sd_beta_b_u[i]) - delta_n)
        }
        for (x_temp in 1:d_j[1]){
          if (acc_delta[x_temp,i]/iter > 0.44){
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) + delta_n)
          }
          else{
            sd_MH_delta[x_temp,i] = exp(log(sd_MH_delta[x_temp,i]) - delta_n)
          }
        }
        for (k in 1:J){
          if (acc_beta_mu_u[i,k]/iter > 0.44){
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) + delta_n)
          }
          else{
            sd_beta_mu_u[i,k] = exp(log(sd_beta_mu_u[i,k]) - delta_n)
          }
        }
      }
      for (d2 in 1:d_j[2]){
        if (acc_b[d2]/iter > 0.44){
          sd_MH_beta_b[d2] = exp(log(sd_MH_beta_b[d2]) + delta_n)
        }
        else{
          sd_MH_beta_b[d2] = exp(log(sd_MH_beta_b[d2]) - delta_n)
        }
        
      }
      
    }
    
    
    # (1) Update of the delta parameter: \delta_{s,i}: MH with log normal 
    #     proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp = which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp = tau[idx_temp]
        cens_temp = cens[idx_temp]
        D_temp = D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        while (delta_prop > min(tau[idx_temp])){
          delta_prop = exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        }
        loglik_prop = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                     b_dat[idx_temp,], 
                                     rep(delta_prop, length(idx_temp)), 
                                     cens_temp, D_temp, TRUE)
        loglik_old = log_likelihood(tau_temp, mu_dat[idx_temp,], 
                                    b_dat[idx_temp,], 
                                    rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                    cens_temp, D_temp, TRUE)
        
        alpha_acc = min(0, loglik_prop + log(delta_prop) -
                          loglik_old - log(delta_old[s_temp,i_temp]))
        l_u = log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] = delta_prop
          delta_dat[idx_temp] = delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] = acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    # (2) Update of mu parameters: \mu_{x,y}^{(i)}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time = T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time = T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time = (k-1):k
      }
      
      for (h in 1:Z_max){ # loop over possible latent values
        # tuples (combinations of covariates) that are clustered together via 
        # the latent h
        idx_cov = matrix(idx_xy[which(z_temp[,k] == h),], length(which(z_temp[,k] == h)), p)
        
        X_1k = unique(idx_cov[,1]) # all possible values of x in this cluster
        X_2k = unique(idx_cov[,2]) # all possible values of y in this cluster
        
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{j,k}: posterior update
          # Pick data with covariate levels of x clustered in group h and 
          # at the correct locations
          idx_i = which( (D[,1] %in% X_1k) & (time %in% idx_time) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[which((D[,1] %in% X_1k) & (time %in% idx_time)),]
          
          
          # Normal proposal distribution centered on the current value
          if (k == 1){
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_mu_star_old[k,h], sd_MH_beta_mu)
          }
          # Normal proposal from the prior
          else if (k == J){
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          # Normal proposal from the prior
          else {
            beta_pre1 = beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_mu_star_prop[,h] = beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] = rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
          
          # Modify the proposed values in the corresponding positions
          for (hh in 1:nrow(idx_cov)){
            B_beta_mu_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] = 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_mu_star_prop[,h]
          }
          
          # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
          mu_prop_dat = exp(B_beta_mu_prop_dat + B_beta_mu_u_dat[idx_i,])
          
          
          if (k == 1){
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          else if (k == J){
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE)
          }
          else {
            beta_post1 = beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop = log_likelihood(tau_temp, mu_prop_dat, 
                                          b_dat[idx_i,], delta_dat[idx_i], 
                                          cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            logpost_old = log_likelihood(tau_temp, mu_dat[idx_i,], 
                                         b_dat[idx_i,], delta_dat[idx_i], 
                                         cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          
          alpha_acc = min(0, logpost_prop - logpost_old)
          l_u = log(runif(1))
          
          if (l_u < alpha_acc){
            
            beta_mu_star_old[k,h] = beta_mu_star_prop[k,h]
            B_beta_mu_dat[idx_i,] = B_beta_mu_prop_dat
            mu_dat[idx_i,] = mu_prop_dat
          }
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_mu_star_old[k,h] = runif(1, low_bound_mu, upp_bound_mu)
        }
      }
    }
    
    # 2(a) Update the constant boundary parameters
    tau_temp = tau
    cens_temp = cens
    D_temp = D
    
    for (d2 in 1:d_j[2]){
      b_prop = rnorm(1, b_old[d2], sd_MH_beta_b[d2])
      
      # Modify the proposed values in the corresponding positions
      B_beta_b_prop_dat = B_beta_b_dat
      B_beta_b_prop_dat[,d2] = b_prop
      
      # This is the proposed value for b
      b_prop_dat = exp(B_beta_b_prop_dat)
      
      logpost_prop = log_likelihood(tau_temp, mu_dat, 
                                    b_prop_dat, delta_dat, 
                                    cens_temp, D_temp, TRUE)
      logpost_old = log_likelihood(tau_temp, mu_dat, 
                                   b_dat, delta_dat, 
                                   cens_temp, D_temp, TRUE)
      
      alpha_acc = min(0, logpost_prop - logpost_old)
      l_u = log(runif(1))
      
      if (l_u < alpha_acc){
        b_old[d2] = b_prop
        # b_dat = b_prop_dat
        B_beta_b_dat = B_beta_b_prop_dat
        b_dat = b_prop_dat
        acc_b[d2] = acc_b[d2] + 1
      }
    }
    
    
    
    # (3) Update the cluster assignments
    for (x_temp in 1:d_j[1]){ # loop over possible latent values
      if (cluster) {
        beta_mess = array(-Inf, dim = c(J, dim_HB))
        beta_mess[J,] = 1/dim_HB
        
        v_old[,J] = H_ball_unif(z_old[[x_temp]][,J], S = Z_max, r = r_HB)
        z_prop[[J]] = H_ball(v_old[,J], S = Z_max, r = r_HB)
        for (k in (J - 1):1){
          idx_i = which( (time == k) & (D[,1] == x_temp) )
          tau_temp = tau[idx_i]
          cens_temp = cens[idx_i]
          D_temp = D[idx_i,]
          
          # (i) Sample the auxiliary variables
          v_temp = H_ball(z_old[[x_temp]][,k], S = Z_max, r = r_HB)
          
          probs = rep(-Inf, dim_HB)
          for (h in 1:dim_HB){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
            
            B_beta_mu_prop_dat = 0.5 * (beta_mu_star_old[k,v_temp[,h]] +
                                          beta_mu_star_old[k+1,z_old[[x_temp]][,k+1]])
            
            mu_dat_prop = exp(t(B_beta_mu_prop_dat + t(B_beta_mu_u_dat[idx_i,])))
            
            probs[h] = g_HB(log_likelihood(tau_temp, mu_dat_prop, b_dat[idx_i,],
                                           delta_dat[idx_i], cens_temp, D_temp, TRUE))
          }
          probs = as.numeric(normalise_log(probs))
          
          v_old[,k] = v_temp[,sample(1:dim_HB, 1, prob = probs)]
          z_prop[[k]] = H_ball(v_old[,k], S = Z_max, r = r_HB)
          
          
          # (ii) Pass messages backwards only in the restricted state space given
          #      by the slice
          z_kp1_temp = which(beta_mess[k+1,] > 0)
          prob_mat = array(-Inf, dim = c(dim_HB, dim_HB))
          
          for (h1 in z_kp1_temp){
            for (h2 in 1:dim_HB){
              
              B_beta_mu_prop_dat = B_beta_mu_dat[idx_i,]
              B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k,z_prop[[k]][,h2]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k,v_old[,k]] +
                                              beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
              
              mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_i,])))
              mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_i,])))
              
              prob_mat[h2,h1] = log(beta_mess[k+1,h1]) -
                0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k,z_prop[[k]][,h2]] -
                                               beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])^2) +
                log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_i,],
                               delta_dat[idx_i], cens_temp, D_temp, TRUE) +
                g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_i,],
                                    delta_dat[idx_i], cens_temp, D_temp, TRUE)) +
                sum(log(Q_F_old[cbind(z_prop[[k]][-x_temp,h2],z_prop[[k+1]][-x_temp,h1])])) +
                log(Q_S_old[z_prop[[k]][x_temp,h2],z_prop[[k+1]][x_temp,h1]])
            }
          }
          if ( sum(is.infinite(sum_rows_log(prob_mat))) == dim_HB){
            beta_mess[k,] = 1/dim_HB
          }
          else{
            beta_mess[k,] = as.numeric(sum_rows_log(prob_mat))
            beta_mess[k,] = as.numeric(normalise_log(beta_mess[k,]))
          }
        }
        
        
        # (iii) Sample states forward (only on allowed states)
        idx_fail = (1:d_j[2])[-x_temp]
        # Sample z_1
        prob_vec = log(beta_mess[1,]) + log(pi_S_0[z_prop[[1]][x_temp,]]) +
          colSums(matrix(log(pi_F_0[z_prop[[1]][-x_temp,]]), d_j[2] - 1, dim_HB))
        prob_vec = as.numeric(normalise_log(prob_vec))
        
        idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
        z_old[[x_temp]][,1] = z_prop[[1]][,idx_samp]
        
        # Sample z_k
        for (k in 2:J){
          idx_km1 = which( (time == k - 1) & (D[,1] == x_temp) )
          tau_temp = tau[idx_km1]
          cens_temp = cens[idx_km1]
          D_temp = D[idx_km1,]
          
          prob_vec = log(beta_mess[k,]) + 
            log(Q_S_old[cbind(z_old[[x_temp]][x_temp,k-1], z_prop[[k]][x_temp,])])
          for (kkk in idx_fail){
            prob_vec = prob_vec + log(Q_F_old[cbind(z_old[[x_temp]][kkk,k-1], z_prop[[k]][kkk,])])
          }
          
          for (z_k_temp in which(is.finite(prob_vec))){
            B_beta_mu_prop_dat = B_beta_mu_dat[idx_km1,]
            
            B_beta_mu_prop_dat_1 = 0.5 * (beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            B_beta_mu_prop_dat_2 = 0.5 * (beta_mu_star_old[k-1,v_old[,k-1]] +
                                            beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
            
            mu_dat_prop_1 = exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_km1,])))
            mu_dat_prop_2 = exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_km1,])))
            
            prob_vec[z_k_temp] = prob_vec[z_k_temp] -
              0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                             beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])^2) +
              log_likelihood(tau_temp, mu_dat_prop_1, b_dat[idx_km1,],
                             delta_dat[idx_km1], cens_temp, D_temp, TRUE) +
              g_HB(log_likelihood(tau_temp, mu_dat_prop_2, b_dat[idx_km1,],
                                  delta_dat[idx_km1], cens_temp, D_temp, TRUE))
          }
          prob_vec = as.numeric(normalise_log(prob_vec))
          
          idx_samp = sample(1:dim_HB, 1, FALSE, prob_vec)
          z_old[[x_temp]][,k] = z_prop[[k]][,idx_samp]
        }
      }
      
      # (4) Assign the cluster specific curves f_{\mu}
      for (y_temp in 1:d_j[2]){
        beta_mu_old[,x_temp,y_temp] = beta_mu_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
      }
      B_beta_mu_dat[(D[,1] == x_temp),] = B[D[,1] == x_temp,] %*% beta_mu_old[,x_temp,]
    }
    z_temp = do.call(rbind, z_old)
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    
    
    # (5) Update the transition probabilities
    pi_S_0 = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + table_int(z_temp[idx_succ,1], Z_max))
    pi_F_0 = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + table_int(z_temp[-idx_succ,1], Z_max))
    tr_count_S = count_assign(z_temp[idx_succ,], Z_max)
    tr_count_F = count_assign(z_temp[-idx_succ,], Z_max)
    for (h in 1:Z_max){
      Q_S_old[h,] = rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                                 tr_count_S[h,])
      Q_F_old[h,] = rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                                 tr_count_F[h,])
    }
    
    
    alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    while ((alpha_S_prop < 0.01) | (alpha_S_prop > 10)){
      alpha_S_prop = exp(rnorm(1, log(alpha_S_old), 0.1))
    }
    alpha_acc = dgamma(alpha_S_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_S_prop) -
      Z_max^2 * lgamma(alpha_S_prop/Z_max) +
      (alpha_S_prop/Z_max - 1) * sum(log(Q_S_old)) - (
        dgamma(alpha_S_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_S_old) -
          Z_max^2 * lgamma(alpha_S_old/Z_max) +
          (alpha_S_old/Z_max - 1) * sum(log(Q_S_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_S_old = alpha_S_prop
      }
    }
    
    alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    while ((alpha_F_prop < 0.01) | (alpha_F_prop > 10)){
      alpha_F_prop = exp(rnorm(1, log(alpha_F_old), 0.1))
    }
    alpha_acc = dgamma(alpha_F_prop, 1, 1, log = T) +
      Z_max * lgamma(alpha_F_prop) -
      Z_max^2 * lgamma(alpha_F_prop/Z_max) +
      (alpha_F_prop/Z_max - 1) * sum(log(Q_F_old)) - (
        dgamma(alpha_F_old, 1, 1, log = T) +
          Z_max * lgamma(alpha_F_old) -
          Z_max^2 * lgamma(alpha_F_old/Z_max) +
          (alpha_F_old/Z_max - 1) * sum(log(Q_F_old)))
    
    l_u = log(runif(1))
    if (!is.na(alpha_acc)){
      if (l_u < alpha_acc){
        alpha_F_old = alpha_F_prop
      }
    }
    
    
    # (6) Correction term for random effects.
    corr_mu = colMeans(beta_mu_u_old)
    beta_mu_u_old = t(t(beta_mu_u_old) - corr_mu)
    
    for (k in 1:J){
      if (k == 1){
        idx_time = which(time == T_min)
      }
      else if (k == J){
        idx_time = which(time == T_max)
      }
      else {
        idx_time = which(time %in% (k-1):k)
      }
      beta_mu_old[k,,] = beta_mu_old[k,,] + corr_mu[k]
      beta_mu_star_old[k,] = beta_mu_star_old[k,] + corr_mu[k]
      B_beta_mu_dat[idx_time,] = B_beta_mu_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_mu_u_dat[idx_time,] = B_beta_mu_u_dat[idx_time,] - 
        as.numeric(B[idx_time,] %*% corr_mu)
    }
    
    
    
    corr_b = mean(beta_b_u_old)
    beta_b_u_old = beta_b_u_old - corr_b
    
    b_old = b_old + corr_b
    B_beta_b_dat = B_beta_b_dat + corr_b
    
    
    # (7) Update the cluster specific smoothness parameters
    RSS_mu = RSS_b = 0
    for (h1 in 1:d_j[1]){
      for (h2 in 1:d_j[2]){
        RSS_mu = RSS_mu + as.numeric(crossprod(beta_mu_old[,h1,h2], P_smooth) %*% 
                                       beta_mu_old[,h1,h2])
      }
    }
    
    sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    while(sigma2_1_mu_temp > 0.2){
      sigma2_1_mu_temp = exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    }
    l_alpha = min(c(0, log(sigma2_1_mu_temp) - 
                      0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_temp) - 
                      0.5/sigma2_1_mu_temp * RSS_mu + 
                      dhalfcauhy(sqrt(sigma2_1_mu_temp), hypers$s_sigma_mu, T) -
                      (log(sigma2_1_mu_old) - 
                         0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_old) -
                         0.5/sigma2_1_mu_old * RSS_mu +
                         dhalfcauhy(sqrt(sigma2_1_mu_old), hypers$s_sigma_mu, T))))
    l_u = log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_mu_old = sigma2_1_mu_temp
    }
    
    
    # (9a) Update random effects for mu
    reff = sample_reff_mu(tau, D, cens, beta_mu_u_old, delta_dat, b_dat,
                          B_beta_mu_dat, mu_dat, B, P_smooth, ind, time,
                          sigma2_mu_us_old, sigma2_mu_ua_old, sd_beta_mu_u,
                          acc_beta_mu_u)
    beta_mu_u_old = reff$beta_u_old
    B_beta_mu_u_dat = reff$B_beta_u_dat
    mu_dat = exp(B_beta_mu_dat + B_beta_mu_u_dat)
    acc_beta_mu_u = reff$acc_beta_u
    
    
    # (9b) Update random effects variances: MH with log normal proposal
    ls_var = sample_smooth_var(sigma2_mu_ua_old, sigma2_mu_us_old,
                               beta_mu_u_old, P_smooth, n_ind)
    sigma2_mu_ua_old = ls_var$sigma2_ua_old
    sigma2_mu_us_old = ls_var$sigma2_us_old
    
    
    for (i in 1:n_ind) {
      beta_b_u_prop = rnorm(1, beta_b_u_old[i], sd_beta_b_u[i])
      # browser()
      idx_i = which(ind == i)
      
      logpr_prop = - 0.5 * ( beta_b_u_prop^2 / sigma2_b_ua_old)
      logpr_old = - 0.5 * ( beta_b_u_old[i]^2 / sigma2_b_ua_old)
      
      b_dat_prop = exp(B_beta_b_dat[idx_i,] + beta_b_u_prop)
      
      
      logpost_prop = log_likelihood(tau[idx_i],
                                    mu_dat[idx_i,],
                                    b_dat_prop,
                                    delta_dat[idx_i],
                                    cens[idx_i],
                                    D[idx_i,], TRUE) +
        logpr_prop
      logpost_old = log_likelihood(tau[idx_i],
                                   mu_dat[idx_i,],
                                   b_dat[idx_i,],
                                   delta_dat[idx_i],
                                   cens[idx_i],
                                   D[idx_i,], TRUE) +
        logpr_old
      
      
      alpha_acc = logpost_prop - logpost_old;
      l_u = log(runif(1))
      if (l_u < alpha_acc){
        beta_b_u_old[i] = beta_b_u_prop
        b_dat[idx_i,] = b_dat_prop
        acc_beta_b_u[i] = acc_beta_b_u[i] + 1
        B_beta_b_u_dat[idx_i,] = beta_b_u_old[i]
      }
    }
    
    sigma2_b_ua_prop = exp(rnorm(1, log(sigma2_b_ua_old), sd = 0.2))
    lu = log(runif(1))
    
    log_prop = - 0.5 * n_ind * log(sigma2_b_ua_prop) -
      0.5 * sum(beta_b_u_old^2) / sigma2_b_ua_prop -
      log(1 + sigma2_b_ua_prop^2)
    log_old = - 0.5 * n_ind * log(sigma2_b_ua_old) -
      0.5 * sum(beta_b_u_old^2) / sigma2_b_ua_old -
      log(1 + sigma2_b_ua_old^2)
    
    alpha = min(c(0, log_prop + log(sigma2_b_ua_prop) -
                    log_old - log(sigma2_b_ua_old)))
    if (lu < alpha){
      sigma2_b_ua_old = sigma2_b_ua_prop
    }
    
    
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # This is the correction for the random effects integration: we need to 
      # compute the variance of the random effects as detailed in the Supplementary
      # Materials
      cov_reff_mu = solve(diag(J) / sigma2_mu_ua_old + P_smooth / sigma2_mu_us_old)
      corr_term_gr_mu = rep(0, nrow(Bgrid))
      corr_term_mu = rep(0, T_max)
      for (k in 1:J){
        corr_term_gr_mu = corr_term_gr_mu + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_mu[,k]))
        corr_term_mu = corr_term_mu + rowSums(B_basis(1:T_max, knots)[,k] * 
                                                t(t(B_basis(1:T_max, knots)) * cov_reff_mu[,k]))
      }
      
      
      
      for (h1 in 1:d_j[1]){
        for (h2 in 1:d_j[2]){
          post_mean_mu[,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + 0.5 * corr_term_gr_mu)
          post_mean_b[,h1,h2,it] = exp(b_old[h2] + 0.5 * sigma2_b_ua_old)
          
          for (i in 1:n_ind){
            post_ind_mu[,i,h1,h2,it] = exp(Bgrid %*% beta_mu_old[,h1,h2] + Bgrid %*% beta_mu_u_old[i,])
            post_ind_b[,i,h1,h2,it] = exp(b_old[h2] + beta_b_u_old[i])
          }
        }
      }
      post_ind_delta[,,it] = delta_old
      
      
      # Sample from the predictive distributions of response category and 
      # response times (so we can use predictive checks to evaluate goodness of fit)
      for (t in 1:T_max){
        for (x_temp in 1:d_j[1]){
          mu_temp = as.numeric(exp(B_basis(t, knots) %*% beta_mu_old[,x_temp,] + 
                                     0.5 * corr_term_mu[t]))
          b_temp = exp(b_old + 0.5 * sigma2_b_ua_old)
          delta_temp = mean(delta_old[x_temp,])
          pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                               b_temp^2)
          pred_ans[t,x_temp,it] = which.min(pred_temp)
          pred_time[t,x_temp,it] = min(pred_temp)
          
          for (i in 1:n_ind){
            mu_temp = as.numeric(exp(B_basis(t, knots) %*% (beta_mu_old[,x_temp,] + 
                                                              beta_mu_u_old[i,])))
            b_temp = exp(b_old + beta_b_u_old[i])
            delta_temp = delta_old[x_temp,i]
            pred_temp = delta_temp + LaplacesDemon::rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                                 b_temp^2)
            pred_ans_ind[i,t,x_temp,it] = which.min(pred_temp)
            pred_time_ind[i,t,x_temp,it] = min(pred_temp)
          }
        }
      }
      
      
      # Save MCMC objects
      post_mean_delta[it,] = rowMeans(delta_old)
      sigma2_mu_us[it] = sigma2_mu_us_old
      sigma2_mu_ua[it] = sigma2_mu_ua_old
      sigma2_1_mu[it] = sigma2_1_mu_old
      sigma2_b_ua[it] = sigma2_b_ua_old
      # z[,,it] = z_temp		
      loglik_chain[it,] = log_likelihood_ind(tau, mu_dat, b_dat, delta_dat, cens, D)
      
      it = it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list(
    'post_mean_delta' = post_mean_delta, 
    'post_mean_mu' = post_mean_mu,
    'post_mean_b' = post_mean_b,
    'post_ind_delta' = post_ind_delta,
    'post_ind_mu' = post_ind_mu,
    'post_ind_b' = post_ind_b,
    'sigma2_mu_us' = sigma2_mu_us, 
    'sigma2_mu_ua' = sigma2_mu_ua,
    'sigma2_b_ua' = sigma2_b_ua,
    'sigma2_1_mu' = sigma2_1_mu, 
    'pred_ans' = pred_ans, 
    'pred_time' = pred_time,
    'pred_ans_ind' = pred_ans_ind, 
    'pred_time_ind' = pred_time_ind, 
    'loglik' = loglik_chain
  ))
}
