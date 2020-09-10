#' Format vaccine model output
#'
#' Take raw odin vaccine model output and formats in long format with the option to select
#' variables and summarise over age groups. Output variables are ordered as in argument ordering.
#'
#' @param x squire_simulation object
#' @param compartments Vector of compartment names, e.g. \code{c("S", "R")}, or sub-compartment names, e.g. \code{c("S", "E1", "E2")}
#' @param summaries Vector of summary names, which may be:
#' \itemize{
#'       \item{"deaths"}{ Deaths per day }
#'       \item{"infections"}{ Infections per day. New infections (note this is currently a slightly different definitionto the main Squire mode)}
#'       \item{"hospitilisations"}{ Hospitalisations per day (Note this takes into account hospital capacity)}
#'       \item{"hospital_occupancy"}{ Occupied Hospital Beds }
#'       \item{"ICU_occupancy"}{ Occupied ICU Beds }
#'       \item{"hospital_demand}{ Required Hospital Beds }
#'       \item{"ICU_demand}{ Required ICU Beds }
#'       \item{"vaccinated"}{ Vaccines administered per day}
#'       }
#' @param reduce_age Collapse age-dimension, calculating the total in the
#'   compartment.
#' @param date_0 Date of time 0 (e.g. "2020-03-01"), if specified a date column will be added
#'
#' @return Formatted long data.frame
#' @export
format <- function(x,
                   compartments = c("S", "E",
                                    "IMild", "ICase", "IICU", "IHospital",
                                    "IRec", "R", "D"),
                   summaries = c("N",
                                 "hospitalisations",
                                 "hospital_demand","hospital_occupancy",
                                 "ICU_demand", "ICU_occupancy",
                                 "vaccines", "unvaccinated", "vaccinated", "priorvaccinated",
                                 "infections", "deaths"),
                   reduce_age = TRUE,
                   date_0 = NULL){

  # Arg checks
  assert_custom_class(x, "nimue_simulation")
  assert_logical(reduce_age)

  # Standardise output dimensions
  if(length(dim(x$output)) == 4){
    x$output <- abind::adrop(x$output, drop = c(FALSE, FALSE, FALSE, TRUE))
  }

  # Get columns indices of variables
  index <- odin_index(x$model)
  if(!all(compartments %in% names(index))){
    stop("Some compartments specified not output by model")
  }

  # Extract time
  time <- x$output[,index$t,1]

  # N replicates
  replicates = dim(x$output)[3]
  # Format over each replicate
  output <- list()
  for(i in 1:replicates){
    output[[i]] <- format_internal(x = x, compartments = compartments, summaries = summaries,
                                   reduce_age = reduce_age, index = index,
                                   time = time, replicate = i)
  }
  output <- dplyr::bind_rows(output)
  # Set levels (order) of output variables
  output$compartment <- factor(output$compartment, levels = c(compartments, summaries))

  # Add date
  if(!is.null(date_0)){
    assert_date(date_0)
    output$date <- as.Date(output$t + as.Date(date_0),
                           format = "%Y-%m-%d")
  }

  # Add age-groups if present
  if("age_index" %in% names(output)){
    ag <- c(paste0(seq(0, 75, 5), "-", seq(5, 80, 5)), "80+")
    output$age_group = factor(ag[output$age_index], levels = ag)
    output <- output  %>%
      dplyr::select(-.data$age_index)
  }

  return(output)
}

#' Internals of Format vaccine model output as data.frame
#' @inheritParams format
#' @param index odin ouput index
#' @param time time vector
#' @param replicate outpu replicate number
format_internal <- function(x, compartments, summaries, reduce_age, index, time,
                            replicate){

  # Convert cumulative outputs
  i_convert <- unlist(index[grepl("_cumu", names(index))])
  x$output[, i_convert, replicate] <- apply(x$output[,i_convert, replicate], 2, function(x){
    x - dplyr::lag(x)
  })
  names(index)[grepl("_cumu", names(index))] <- sapply(strsplit(names(index)[grepl("_cumu", names(index))], "_"), `[`, 1)

  # Select variables and summary outputs
  get <- c(compartments, summaries)
  get <- get[get %in% names(index)]
  i_select <- index[get]
  # Select outputs, collapsing over vaccine dimension where required
  o <- lapply(i_select, function(x, y){
    if(is.matrix(x)){
      apply(x, 1, function(a, b){
        rowSums(b[,a,replicate])
      }, b = y)
    } else {
      y[,x,replicate]
    }
  }, y = x$output)

  # Collapse age
  if(reduce_age){
    o <- lapply(o, collapse_age)
  } else {
    o <- lapply(o, add_age)
  }

  # Add names of variables
  for(i in 1:length(o)){
    o[[i]] <- data.frame(o[[i]]) %>%
      dplyr::mutate(compartment = names(o)[i])
  }

  # Add time and replicate columns
  o <- dplyr::bind_rows(o) %>%
    dplyr::mutate(t = rep(time, dplyr::n() / length(time)),
                  replicate = replicate)

  return(o)
}

#' Keep age groups
#'
#' @param x age-disaggregated odin output matrix
#'
#' @return age-disaggregated output matrix
add_age <- function(x){
  m <- matrix(c(rep(1:ncol(x), each = (nrow(x))), as.vector(x)), ncol = 2)
  colnames(m) <- c("age_index", "value")
  return(m)
}

#' Collapse age groups
#'
#' @param x age-disaggregated odin output matrix
#'
#' @return age-aggregated output matrix
collapse_age <- function(x){
  m <- matrix(rowSums(x), ncol = 1)
  colnames(m) <- "value"
  return(m)
}

## Index locations of outputs in odin model
#' @noRd
odin_index <- function(model) {
  n_out <- environment(model$initialize)$private$n_out %||% 0
  n_state <- length(model$initial(0))
  model$transform_variables(seq_len(1L + n_state + n_out))
}

#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

