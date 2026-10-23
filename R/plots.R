#' Descriptive plots
#'
#' Proportion of each decision for every input stimulus and block, by subject
#' (thin lines) and pooled over subjects (thick lines). Censored trials are
#' excluded.
#'
#' @param data data frame in the format described in [LDDMM()].
#' @return A `ggplot` object.
#' @export
plot_accuracy = function(data){
  block = d = s = cens = subject = n = freq = NULL

  proportions = function(...){
    data |>
      dplyr::filter(cens == 0) |>
      dplyr::count(..., block, s, d) |>
      dplyr::group_by(..., block, s) |>
      dplyr::mutate(freq = n / sum(n), s = factor(s), d = factor(d)) |>
      dplyr::ungroup()
  }
  by_subject = proportions(subject) |> dplyr::mutate(subject = factor(subject))
  pooled = proportions()

  ggplot(by_subject) +
    geom_point(aes(x = block, y = freq, group = subject, col = d), alpha = 0.2) +
    geom_line(aes(x = block, y = freq, group = interaction(subject, d), col = d), alpha = 0.2) +
    geom_line(aes(x = block, y = freq, col = d), linewidth = 1.5, data = pooled) +
    facet_wrap( ~ s, nrow = 2) +
    labs(x = 'block', y = 'probability') +
    scale_color_brewer(name = "", palette = 'Set1') +
    scale_x_continuous(breaks = sort(unique(data$block))) +
    theme(legend.position = "top")
}


#' Descriptive plots
#'
#' Mean response time of each decision for every input stimulus and block,
#' pooled over subjects. Censored trials are excluded.
#'
#' @inheritParams plot_accuracy
#' @return A `ggplot` object.
#' @export
plot_RT = function(data){
  mean_r_time = block = d = s = cens = r_time = NULL

  data |>
    dplyr::filter(cens == 0) |>
    dplyr::group_by(block, s, d) |>
    dplyr::summarise(mean_r_time = mean(r_time), .groups = "drop") |>
    ggplot() +
    geom_line(aes(x = block, y = mean_r_time, col = factor(d)), linewidth = 1.5) +
    facet_wrap( ~ s) +
    labs(x = 'block', y = 'response time') +
    scale_color_brewer(name = "", palette = "Set1") +
    scale_x_continuous(breaks = sort(unique(data$block))) +
    theme(legend.position = "top")
}
