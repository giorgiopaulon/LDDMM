
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
