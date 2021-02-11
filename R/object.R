#' nimue simulation plot
#'
#' @param x An squire_simulation object
#' @param replicates Plot replicates
#' @param summarise Logical, add summary line
#' @param ci logical add confidence interval ribbon
#' @param q Quantiles for upper and lower of interval ribbon
#' @param var_select Vector of variable names to plot (default is all)
#' @param summary_f Function to summarise each compartment
#'   passed to the \code{fun} argument of \code{\link[ggplot2]{stat_summary}}
#' @param x_var X variable to use for plotting (default is \code{"t"},
#'   but can be set to, \code{"date"}, if \code{date_0} provided), which will
#'   cause the date to be plotted rather than time.
#' @param particle_fit If the squire_simulation provided is the result of
#'   running the particle filter, do we want to just plot the fit. Default =
#'   FALSE
#' @param date_0 Date of time 0 (e.g. "2020-03-01"), if specified a date column
#'   will be added
#' @param ... additional arguments affecting the plot produced.
#'
#' @export
#'
plot.nimue_simulation <- function(x,
                                  var_select = NULL,
                                  replicates = FALSE,
                                  summarise = TRUE,
                                  ci = TRUE,
                                  q = c(0.025, 0.975),
                                  summary_f = mean,
                                  x_var = "t",
                                  date_0 = NULL,
                                  particle_fit = FALSE,
                                  ...) {

  # work out what compartments are being plotted
  compartments = c("S", "E",
                   "IMild", "ICase", "IICU", "IHospital",
                   "IRec", "R", "D")
  summaries = c("N",
                "hospitalisations",
                "hospital_demand","hospital_occupancy",
                "ICU_demand", "ICU_occupancy",
                "vaccines", "unvaccinated", "vaccinated", "priorvaccinated",
                "infections", "deaths")

  comps <- var_select[var_select %in% compartments]
  summs <- var_select[var_select %in% summaries]

  # get the compartments requried
  pd <- do.call(rbind, lapply(seq_len(dim(x$output)[3]), function(i) {
    format(x, compartments = comps, summaries = summs, replicate = i)
    })) %>%
    dplyr::rename(y = .data$value)

  # replacing time with date if date_0 is provided
  if(!is.null(date_0)){
    assert_date(date_0)
    pd$date <- as.Date(pd$t + as.Date(date_0),
                        format = "%Y-%m-%d")
  }

  # make the x label be the x axis requested
  pd <- pd %>% dplyr::mutate(x = .data[[x_var]])


  # remove any NA rows (due to different start dates)
  if(sum(is.na(pd$t) | is.na(pd$y))>0) {
    pd <- pd[-which(is.na(pd$t) | is.na(pd$y)),]
  }

  # Format summary data
  pds <- pd %>%
    dplyr::group_by(.data$x, .data$compartment) %>%
    dplyr::summarise(ymin = stats::quantile(.data$y, q[1]),
                     ymax = stats::quantile(.data$y, q[2]),
                     y = summary_f(.data$y))


  # Plot
  p <- ggplot2::ggplot()

  # Add lines for individual draws
  if(replicates){
    p <- p + ggplot2::geom_line(data = pd,
                                ggplot2::aes(x = .data$x,
                                             y = .data$y,
                                             col = .data$compartment,
                                             group = interaction(.data$compartment, .data$replicate)),
                                alpha = max(0.2, 1 / x$parameters$replicates))
  }

  if(summarise){
    if(x$parameters$replicates < 10){
      warning("Summary statistic estimated from <10 replicates")
    }
    p <- p + ggplot2::geom_line(data = pds,
                                ggplot2::aes(x = .data$x, y = .data$y,
                                             col = .data$compartment))
  }

  if(ci){
    if(x$parameters$replicates < 10){
      warning("Confidence bounds estimated from <10 replicates")
    }
    p <- p + ggplot2::geom_ribbon(data = pds,
                                  ggplot2::aes(x = .data$x,
                                               ymin = .data$ymin,
                                               ymax = .data$ymax,
                                               fill = .data$compartment),
                                  alpha = 0.25, col = NA)
  }

  # Add remaining formatting
  p <- p +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::scale_fill_discrete(guide = FALSE) +
    ggplot2::xlab("Time") +
    ggplot2::ylab("N") +
    ggplot2::theme_bw() +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 2))

  return(p)
}
