# LDDMM

An R package for Longitudinal Drift-Diffusion Mixed Models (LDDMM), v0.1.

### Overview

Codes accompanying "Bayesian Semiparametric Longitudinal Drift-Diffusion Mixed Models for Tone Learning in Adults" by Paulon, Llanos, Chandrasekaran, Sarkar.

This package implements a novel generic framework for longitudinal functional mixed models that allows automated assessment of an associated predictor's local time-varying influence. We build on this to develop a novel inverse-Gaussian drift-diffusion mixed model for multi-alternative decision-making processes in longitudinal settings. Our proposed model and associated computational machinery make use of B-spline mixtures, hidden Markov models (HMM) and factorial hidden Markov models (fHMM), locally informed Hamming ball samplers etc. to address statistical challenges.

The main function is `LDDMM`; please see the corresponding help files and vignette for details, as well as the main article:

Paulon, G., Llanos, F., Chandrasekaran, B., Sarkar, A. (2021). [Bayesian semiparametric longitudinal drift-diffusion mixed models for tone learning in adults](https://www.tandfonline.com/doi/abs/10.1080/01621459.2020.1801448?journalCode=uasa20). Journal of the American Statistical Association **116**, 1114-1127

### Installation

To install the package in R, first install the `devtools` package, and then use the commands
`````````
library(devtools)
install_github('giorgiopaulon/lddmm')
`````````

### Contact

For bug reporting purposes, e-mail Giorgio Paulon (giorgio.paulon@utexas.edu).
