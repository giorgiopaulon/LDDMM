#' Drift Diffusion Model Fit
#'
#' Fits the longitudinal inverse Gaussian drift-diffusion mixed model of Paulon
#' et al. (2021) by Markov chain Monte Carlo.
#'
#' On trial \eqn{\ell} of block \eqn{t}, subject \eqn{i} hears input stimulus
#' \eqn{s} and answers \eqn{d}. One accumulator per possible decision races to
#' its boundary; the first to arrive determines the decision and the response
#' time. The accumulator of decision \eqn{d} has drift
#' \eqn{\mu^{(i)}_{d,s}(t) = \exp\{f_{\mu,d,s}(t) + u^{(i)}_{\mu}(t)\}}, boundary
#' \eqn{b^{(i)}_{d,s}(t) = \exp\{f_{b,d,s}(t) + u^{(i)}_{b}(t)\}} and starts after
#' an offset \eqn{\delta^{(i)}_s}. The fixed effects \eqn{f} are quadratic
#' B-splines over blocks whose coefficients are locally clustered across
#' \eqn{(d, s)} by hidden Markov models; the random effects \eqn{u} are
#' smooth subject-specific spline curves.
#'
#' Priors are calibrated so that, for the most stable performance, most values
#' of `r_time` should lie between 0 and 10; response times in seconds work well.
#' The offset \eqn{\delta^{(i)}_s} cannot exceed the fastest response of
#' subject \eqn{i} to stimulus \eqn{s}, so implausibly fast responses (for
#' example anticipations under 100 ms) force it towards zero and are best
#' removed before fitting.
#'
#' @param data data frame with one row per trial and the columns
#' * `subject`: participant labels
#' * `block`: training blocks (longitudinal units), consecutive integers
#' * `s`: input stimuli, coded 1, ..., K
#' * `d`: decisions, coded 1, ..., K like `s` (may be `NA` for censored trials)
#' * `r_time`: response times (positive), e.g. in seconds
#' * `cens`: censoring indicators (1 censored, 0 observed); for a censored
#'   trial `r_time` is the time at which observation stopped
#' @param hypers list of hyperparameters (missing entries take default values):
#' * `s_sigma_mu`, `s_sigma_b`: scales of the half-Cauchy priors on the
#'   standard deviations \eqn{\sigma_{1}} that control the smoothness of the
#'   drift and boundary curves (default 0.1)
#' * `sigma2_1_max`, `sigma2_us_max`: upper truncation of the priors on the
#'   smoothness variances of the fixed-effect curves and of the random-effect
#'   curves (default 0.2)
#' * `sigma2_ua_max`: upper truncation of the half-Cauchy prior on the variance
#'   of the random effects (default `Inf`, no truncation)
#' * `pseudo_sd`: spread of the pseudo-prior of unused clusters (default 0.1).
#'   Coefficients of clusters that no curve uses are drawn close to those of
#'   the clusters in use, so that curves can split off into a new cluster.
#'   This only affects how well the sampler mixes, not the posterior.
#' * `mu_beta0_mu`, `sigma2_beta0_mu`, `mu_beta0_b`, `sigma2_beta0_b`: mean and
#'   variance of the normal prior on spline coefficients that start a curve or
#'   belong to an unused cluster (article eq. 9). By default they are set from
#'   moment estimates in every (block, stimulus, decision) cell.
#' @param cluster should the fixed-effect coefficients be locally clustered
#'   across input-response combinations? With `FALSE` every combination has its
#'   own curve.
#' @param boundaries boundary model:
#' * `"flexible"`: \eqn{b_{d,s}(t)} varies over time and is clustered jointly
#'   with the drifts
#' * `"fixed"`: \eqn{b_{d}(t)} varies over time but not across stimuli
#' * `"constant"`: \eqn{b_{d,s}} is constant over time
#' * `"fixed-constant"`: \eqn{b_{d}} is constant over time and across stimuli
#'
#' The random effects of the boundary are spline curves for the first two
#' models and scalars for the last two.
#' @param random_effects `"correct-incorrect"` gives every subject one
#'   random-effect curve for correct responses (\eqn{d = s}) and another for
#'   incorrect responses (article §3.2.2); `"shared"` uses a single curve for both.
#' @param z_max maximum number of clusters at each knot when `cluster = TRUE`.
#'   Defaults to `min(K^2, 2K)`, i.e. 8 for 4 categories as in the article.
#' @param informed use locally informed proposals in the Hamming ball sampler
#'   of the cluster labels (supplement §S.4)?
#' @param Niter total number of iterations
#' @param burnin burnin of the chain; proposal scales are tuned only during burn-in
#' @param thin thinning factor
#' @param save_z store the cluster labels of every saved draw (element `Z`)?
#' @param save_loglik store the log-likelihood of every observation for every
#'   saved draw (element `loglik`, needed by [compute_WAIC()])?
#' @param verbose show a progress bar?
#' @param debug_checks verify after every update that the sampler's cached
#'   quantities match their definitions (slow; for testing).
#' @return List with the following MCMC posterior samples:
#' * `post_mean_delta`: population offsets (draw x stimulus)
#' * `post_mean_mu`, `post_mean_b`: population drift and boundary curves
#'   (block x stimulus x decision x draw), integrating out the random effects
#' * `post_ind_delta`: individual offsets (stimulus x subject x draw)
#' * `post_ind_mu`, `post_ind_b`: individual drift and boundary curves
#'   (block x subject x stimulus x decision x draw)
#' * `sigma2_mu_us`, `sigma2_mu_ua`, `sigma2_b_us`, `sigma2_b_ua`: smoothness
#'   and variance parameters of the random effects
#' * `sigma2_1_mu`, `sigma2_1_b`: smoothness parameters of the fixed effects
#'   (`sigma2_b_us` and `sigma2_1_b` are `NA` for constant boundaries)
#' * `pred_ans`, `pred_time`: posterior predictive decisions and response
#'   times at the population level (block x stimulus x draw)
#' * `pred_ans_ind`, `pred_time_ind`: the same for each subject
#'   (subject x block x stimulus x draw)
#' * `loglik`: log-likelihood of each observation (draw x observation)
#' * `Z`: cluster labels (stimulus x decision x knot x draw), if `save_z = TRUE`
#' * `acceptance`: acceptance rates of the Metropolis-Hastings steps after burn-in
#' * `settings`: model options and hyperparameters used
#' @references Paulon, G., Llanos, F., Chandrasekaran, B., Sarkar, A. (2021).
#'   Bayesian semiparametric longitudinal drift-diffusion mixed models for tone
#'   learning in adults. Journal of the American Statistical Association 116,
#'   1114-1127. \doi{10.1080/01621459.2020.1801448}
#' @export
LDDMM = function(data, hypers = NULL, cluster = TRUE,
                 boundaries = c('flexible', 'constant', 'fixed', 'fixed-constant'),
                 random_effects = c('correct-incorrect', 'shared'),
                 z_max = NULL, informed = TRUE,
                 Niter = 5000, burnin = 2000, thin = 5,
                 save_z = FALSE, save_loglik = TRUE, verbose = TRUE,
                 debug_checks = FALSE){
  boundaries = match.arg(boundaries)
  random_effects = match.arg(random_effects)
  if (burnin >= Niter || thin < 1) stop("Need burnin < Niter and thin >= 1.", call. = FALSE)

  data = validate_data(data)
  design = build_design(data, boundaries, random_effects, cluster, z_max, informed)
  st = init_state(design, hypers)
  run_sampler(st, Niter, burnin, thin, save_z, save_loglik, verbose, debug_checks)
}
