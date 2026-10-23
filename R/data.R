#' Example dataset
#'
#' Data in the format required by [LDDMM()]: 24,254 trials of 99 participants
#' over 5 training blocks, with two possible input stimuli and decisions. The
#' data were analyzed in Roark, Paulon, Sarkar and Chandrasekaran (2021),
#' Psychonomic Bulletin & Review 28, 898-909.
#'
#' * `subject`: participant labels
#' * `block`: training blocks (longitudinal units)
#' * `s`: input stimuli
#' * `d`: decisions
#' * `r_time`: response times, in seconds
#' * `cens`: censoring indicators (1 censored, 0 non censored); all 0 here
#'
#' @format A data frame with 24,254 rows and 6 columns.
"data"
