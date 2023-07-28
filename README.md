
# LDDMM

An R package for Longitudinal Drift-Diffusion Mixed Models (LDDMM),
v0.4.

**Authors**: [Giorgio Paulon](https://giorgiopaulon.github.io/), [Abhra
Sarkar](https://abhrastat.github.io/)

### Overview

Codes accompanying “Bayesian Semiparametric Longitudinal Drift-Diffusion
Mixed Models for Tone Learning in Adults” by Paulon, Llanos,
Chandrasekaran, Sarkar.

This package implements a novel generic framework for longitudinal
functional mixed models that allows automated assessment of an
associated predictor’s local time-varying influence. We build on this to
develop a novel inverse-Gaussian drift-diffusion mixed model for
multi-alternative decision-making processes in longitudinal settings.
Our proposed model and associated computational machinery make use of
B-spline mixtures, hidden Markov models (HMM) and factorial hidden
Markov models (fHMM), locally informed Hamming ball samplers etc. to
address statistical challenges.

The main function is `LDDMM`; please see the following vignette for
details, as well as the main article:

Paulon, G., Llanos, F., Chandrasekaran, B., Sarkar, A. (2021). [Bayesian
semiparametric longitudinal drift-diffusion mixed models for tone
learning in
adults](https://doi.org/10.1080/01621459.2020.1801448).
Journal of the American Statistical Association **116**, 1114-1127

The data included in this package was analyzed in:

Roark, C. L., Paulon, G., Sarkar, A., Chandrasekaran, B. (2021).
[Comparing perceptual category learning across modalities in the same
individuals](https://link.springer.com/article/10.3758/s13423-021-01878-0).
Psychonomic Bulletin & Review **28**, 898-909

and is available [here](https://osf.io/msnq2/).

### Installation

To install the package in R, first install the `devtools` package, and
then use the commands

``` r
library(devtools)
install_github('giorgiopaulon/lddmm')
```

If you are using a Windows machine, you might have to also install and
configure `Rtools` using the following
[instructions](https://CRAN.R-project.org/bin/windows/Rtools/).

### Usage

The following is a minimal example of a simple model fit. 
For numerical stability, the unit of measurement should be such that the numerical values of most response times should lie in $[0, 10]$.


``` r
# Load libraries
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(reshape2)
library(latex2exp)
library(lddmm)

theme_set(theme_bw(base_size = 14))
cols <- brewer.pal(9, "Set1")

# Load the data
data('data')

# Descriptive plots
plot_accuracy(data)
plot_RT(data)

# Run the model
hypers <- NULL
hypers$s_sigma_mu <- hypers$s_sigma_b <- 0.1

# Change the number of iterations when running the model
# Here the number is small so that the code can run in less than 1 minute
Niter <- 25
burnin <- 15
thin <- 1
samp_size <- (Niter - burnin) / thin

set.seed(123)
fit <- LDDMM(data = data, 
             hypers = hypers, 
             Niter = Niter, 
             burnin = burnin, 
             thin = thin)

# Plot the results
plot_post_pars(data, fit, par = 'drift')
plot_post_pars(data, fit, par = 'boundary')

# Compute the WAIC to compare models
compute_WAIC(fit)
```


To extract relevant posterior draws or posterior summaries instead of simply plotting them, one can use the functions `extract_post_mean` or `extract_post_draws`. 
Auxiliary functions that assume constant boundary parameters over time ($b_{d,s}^{(i)}(t) = b_{d,s} + u_{d,s}^{(i)}$ using the article notation), fix the boundaries to the same level across input predictors ($b_{d,s}^{(i)}(t) = b_{d}(t) + u^{(i)}_{d}(t)$ using the article notation), or both ($b_{d,s}^{(i)}(t) = b_{d} + u_{d}^{(i)}$ using the article notation) can be called with the options `boundaries = "constant"`, `boundaries = "fixed"`, and `boundaries = "fixed-constant"` respectively.



### Questions or bugs

For bug reporting purposes, e-mail Giorgio Paulon
(<giorgio.paulon@utexas.edu>).

### Citation

Please cite the following publication if you use this package in your
research: Paulon, G., Llanos, F., Chandrasekaran, B., Sarkar, A. (2021).
[Bayesian semiparametric longitudinal drift-diffusion mixed models for
tone learning in
adults](https://doi.org/10.1080/01621459.2020.1801448).
Journal of the American Statistical Association **116**, 1114-1127
