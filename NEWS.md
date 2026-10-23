# LDDMM 1.0.0

This major release introduces an improved sampler.

## Model and defaults

* Every participant has one random effect curve for correct responses and one for incorrect responses (`random_effects = "correct-incorrect"`, the new default). A single curve for both is available with `random_effects = "shared"`.
* The maximum number of clusters at each block defaults to `min(K^2, 2K)` (8 for 4 categories) and can be set with `z_max`.
* Spline coefficients that start a curve have the normal prior of article eq. (9), with mean and variance set from the data (`hypers$mu_beta0_*`, `hypers$sigma2_beta0_*`).
* The smoothness of the fixed effects has a half-Cauchy prior on its standard deviation (scale `hypers$s_sigma_*`). The upper limits on the variance parameters are documented hyperparameters (`hypers$sigma2_1_max`, `hypers$sigma2_us_max`, `hypers$sigma2_ua_max`).

## Improved sampler

* A single implementation serves all boundary models, and the cluster-label and random-effect updates run in C++. Fits are about 12-21 times faster on the tone learning data.
* The cluster labels are updated with an exact locally informed Hamming ball sampler, and the cluster-specific spline coefficients with proposals built from the prior of article eq. (9).
* New joint moves improve mixing: a translation between fixed and random effects, a rescaling of the random effects together with their variances, and shifts of the offsets together with the drift and boundary levels.
* Coefficients of unused clusters are drawn close to those of the clusters in use (`hypers$pseudo_sd`), so that curves move between clusters more easily.
* Random effects start from each participant's moment estimates, and proposal scales are tuned during burn-in only, in batches of 50 iterations.
* Response-time distributions are evaluated on the log scale, which keeps them accurate for slow responses.

## Other changes

* New `plot_predictive()` for posterior predictive checks of the decisions and the response times.
* Two worked examples are installed with the package (`system.file("examples", package = "lddmm")`): `toy_example.R` analyses the bundled data, and `toy_simulation.R` simulates data with known parameters and compares the posterior with the truth. The README and the vignette walk through fitting, posterior predictive checks and convergence diagnostics.
* New arguments `save_z` (store cluster labels), `save_loglik`, `informed`, `verbose` and `debug_checks`; the fit also reports acceptance rates and the settings used.
* `LDDMM()` checks its input and reports problems with the data.
* `log_likelihood_ind()` returns a plain vector, and `compute_WAIC()` is computed on the log scale.
* `plot_accuracy()` works without attaching dplyr.
* The package no longer depends on gtools, LaplacesDemon, plyr, RColorBrewer, RcppProgress and rgen.
* Unit tests, and validation scripts in `tests/validation/` (a Geweke joint distribution test, invariance tests of individual updates, and a simulate-and-recover study).
* The documentation states that response times are in seconds.

# LDDMM 0.4.2

## Major changes

* Fixed a bug in the version `boundaries = "fixed"`.
* Added the function input `cluster` to determine whether to use the clustering model or the unrestricted model.

# LDDMM 0.4.1

## Major changes

* Added the option for constant and fixed boundaries over time (`boundaries = "fixed-constant"`).

# LDDMM 0.4.0

## Major changes

* Now the option for constant boundaries over time (`boundaries = "constant"`) varies across input tone and decision. The same is true for the random effects. 
* Now the option for fixed boundaries (`boundaries = "fixed"`) varies across time and is fixed across input tone s. The same is true for the random effects.

# LDDMM 0.3.0

## Major changes

* Now the option for constant boundaries over time (`boundaries = "constant"`) includes same boundaries per each input tone, and varies across possible decisions
* Added a function to calculate the Watanabe-Akaike information criterion to compare model performance

# LDDMM 0.2.1


## Minor changes

* Updated documentation


# LDDMM 0.2.0

## Major changes

* Now the boundary parameter in the model can either be unrestricted (`boundaries = "flexible"`), constant over time (`boundaries = "constant"`) or fixed at the same value across predictors (`boundaries = "fixed"`)

## Minor changes

* Added the DOI for the article reference
* Clarified that response times should be the log transformed of time (in milliseconds) for numerical stability

## Bug fixes


# LDDMM 0.1.0

* First version
