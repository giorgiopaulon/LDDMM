

#' Example dataset
#'
#' A toy dataset in the correct format for the LDDMM function call. This dataset
#' has two possible response categories.
#'
#' * `subject`: vector of size n containing the participant labels
#' * `block`: vector of size n containing the training blocks (longitudinal units)
#' * `s`: vector of size n containing the stimuli
#' * `d`: vector of size n containing the decisions
#' * `r_time`: vector of size n containing the response times (log transformed)
#' * `cens`: vector of size n containing the censoring indicators (1 censored, 0 non censored)
#'
#' @format A data frame with 24,254 rows and 6 columns
#'
"data"

