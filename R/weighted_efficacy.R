#' Assign N doses to age groups based on weightings by number of people eligible to be vaccinated
#'
#' @param doses Total available doses
#' @param target_pop Number of people eligible to be vaccinated
assign_doses <- function(doses, target_pop){
  if(sum(target_pop) < doses){
    return(target_pop)
  } else{
    group_weights <- target_pop / sum(target_pop)
    assigned <- floor(doses * group_weights)
    if(sum(assigned) != doses){
      assigned <- assigned + (rank(-group_weights, ties.method = "last") <= (doses %% length(group_weights)))
    }
  }
  return(assigned)
}

#' Estimate coverage by age group
#'
#' @param dose_times  Dose monitoring data
#' @param dose_number  Dose number: 1 or 2
coverage <- function(dose_times, dose_number){
  sapply(dose_times, function(x){
    mean(!is.na(x[,dose_number]))
  })
}

#' Target number of individual in each age group to vaccinate
#'
#' @param dose_number Dose number: 1 or 2
#' @param dose_times Dose monitoring data
#' @param prioritisation Current row of the prioritisation matrix
#' @param t Current time
#' @param dose_period Days between dose 1 and dose 2
#' @param d2_prioritise Boolean vector indicating which age groups to prioritise for dose 2
target_pop <- function(dose_number, dose_times, prioritisation, t, dose_period, d2_prioritise){
  # Current coverage of specified dose number
  current_coverage <- coverage(dose_times, dose_number)
  # Remaining population left to cover with current dose number to reach target coverage in prioritisation step
  n_to_cover <- ceiling(pmax(0, (prioritisation - current_coverage)) * sapply(dose_times, nrow))
  if(dose_number == 2){
    # Population eligible for current dose number - must be a certain period after dose 1
    eligable <- sapply(eligable_for_second(dose_times, t, dose_period), sum)
    n_to_cover <- pmin(n_to_cover, eligable) * d2_prioritise
  }
  return(n_to_cover)
}

#' Identify those eligible for a second dose
#'
#' Those eligible are those who have had a first dose more than a given period prior and
#' who have not yet received a second dose.
#'
#' @param dose_times  Dose monitoring data
#' @param t Current time
#' @param dose_period  Days between dose 1 and dose 2
eligable_for_second <- function(dose_times, t, dose_period){
  lapply(dose_times, function(x){
    had_first_beyond_threshold <- t - x[, 1] >= dose_period
    had_first_beyond_threshold[is.na(had_first_beyond_threshold)] <- FALSE
    not_had_2nd <- is.na(x[, 2])
    had_first_beyond_threshold & not_had_2nd
  })
}

#' Administer first doses
#'
#' @param ad Number of doses to administer in each age group
#' @param dose_times Dose monitoring data
#' @param t Current time
administer_first_dose <- function(ad, dose_times, t){
  for(i in seq_along(ad)){
    if(ad[i] > 0){
      to_jab <- which(is.na(dose_times[[i]][,1]))[1:ad[i]]
      dose_times[[i]][to_jab,1] <- t
    }
  }
  return(dose_times)
}

#' Administer second doses
#'
#' @param ad Number of doses to administer in each age group
#' @param dose_times Dose monitoring data
#' @param dose_period Days between dose 1 and dose 2
#' @param t Current time
administer_second_dose <- function(ad, dose_times, dose_period, t){
  eligable <- eligable_for_second(dose_times, t, dose_period)
  for(i in seq_along(ad)){
    if(ad[i] > 0){
      to_jab <- which(eligable[[i]])[1:ad[i]]
      dose_times[[i]][to_jab,2] <- t
    }
  }
  return(dose_times)
}


#' Extract dose numbers administers
#'
#' @param dose_times Dose monitoring data
#' @inheritParams weighted_efficacy
extract_dose_number <- function(dose_times, maxt){
  output_1dose <- matrix(0, nrow = maxt, ncol = length(dose_times))
  colnames(output_1dose) <- 1:length(dose_times)
  output_2dose <- matrix(0, nrow = maxt, ncol = length(dose_times))
  colnames(output_2dose) <- 1:length(dose_times)
  for(t in 1:maxt){
    for(a in 1:length(dose_times)){
      output_1dose[t, a] <- sum(dose_times[[a]][,1] == t, na.rm = TRUE)
      output_2dose[t, a] <- sum(dose_times[[a]][,2] == t, na.rm = TRUE)
    }
  }
  output_1dose <- as.data.frame(output_1dose) %>%
    dplyr::mutate(t = 1:maxt) %>%
    tidyr::pivot_longer(-.data$t, names_to = "age_group", values_to = "n_1dose") %>%
    dplyr::mutate(age_group = factor(.data$age_group, levels = 1:length(dose_times)))
  output_2dose <- as.data.frame(output_2dose) %>%
    dplyr::mutate(t = 1:maxt) %>%
    tidyr::pivot_longer(-.data$t, names_to = "age_group", values_to = "n_2dose") %>%
    dplyr::mutate(age_group = factor(.data$age_group, levels = 1:length(dose_times)))
  output <- dplyr::left_join(output_1dose, output_2dose, by = c("t", "age_group"))
  return(output)
}


#' Extract weighted efficacies
#'
#' @param dose_number  Dose number: 1 or 2
#' @inheritParams weighted_efficacy
add_weighted_efficacy <- function(dose_number, infection_efficacy, disease_efficacy, v1v2){
  dose_number  %>%
    dplyr::group_by(.data$age_group) %>%
    dplyr::mutate(c2 = cumsum(.data$n_2dose),
                  c1 = cumsum(.data$n_1dose) - .data$c2) %>%
    # Number vaccine protected - lag between administration of dose 1 and protection
    dplyr::mutate(protected1 = dplyr::lag(.data$c1, v1v2, default = 0),
                  protected2 = .data$c2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(weighted_infection_efficacy = purrr::map2_dbl(.data$protected1, .data$protected2, function(x, y, efficacy){
      ifelse(x + y > 0, stats::weighted.mean(efficacy, c(x, y)), efficacy[1])
    }, efficacy = infection_efficacy),
    weighted_disease_efficacy = purrr::map2_dbl(.data$protected1, .data$protected2, function(x, y, efficacy){
      ifelse(x + y > 0, stats::weighted.mean(efficacy, c(x, y)), efficacy[1])
    }, efficacy = disease_efficacy)) %>%
    dplyr::select(-.data$protected1, -.data$protected2, -.data$c1, -.data$c2)
}

#' Estimate weighted efficacies
#'
#' Estimate weighted age-disaggregated efficacies and numbers of doses over time
#' to facilitate an approximation of 2 doses with varying 1st and 2nd dose efficacy.
#' Current assumptions are that the duration in v1 and v2 (vaccinated but not protected)
#' is fixed and that the number individuals that cannot be vaccinated (due to being
#' in the hospitalised flow) does not substantivly impact vaccine roll out rates. It is
#' recommended to run with reasonably large \code{N} to ensure that there are enough
#' individuals assigned to each age group to get stable results.
#'
#' @param iso3c Country iso3c code
#' @param N Number to simulate. Country population is rescaled to this number
#' @param maxt Max time
#' @param doses_per_day Vector of available doses per day
#' @param dose_period Days between dose 1 and dose 2
#' @param v1v2 Duration in v1v2 (vaccinated but not protected)
#' @param prioritisation_matrix Prioritisation matrix
#' @param d2_prioritise Boolean vector indicating which age groups are prioritised to receive 2nd dose
#' @param infection_efficacy Vector of length 2 of infection efficacy for first and second dose
#' @param disease_efficacy Vector of length 2 of disease efficacy for first and second dose
#'
#' @return Age disaggregated dose numbers and efficacies over time.
#' @export
weighted_efficacy <- function(iso3c,
                              N,
                              maxt,
                              doses_per_day,
                              dose_period,
                              v1v2,
                              prioritisation_matrix,
                              d2_prioritise,
                              infection_efficacy,
                              disease_efficacy){

  stopifnot(v1v2 < dose_period)
  # Get population
  pop <- squire::get_population(iso3c = iso3c)
  # Rescale population and define age group for each individual
  ages <- rep(1:nrow(pop), round(N * (pop$n / sum(pop$n))))
  # Record of when dose 1 and 2 were given for each individual
  dose_times <- data.frame(dose1_t = rep(NA, length(ages)), dose2_t = rep(NA, length(ages)))
  dose_times <- split(dose_times, ages)
  # Start with phase 1
  phase <- 1
  # Start at first row of prioritisation matrix
  p_step <- 1
  for(t in 1:maxt){
    # Available doses to distribute
    current_doses <- doses_per_day[t]

    # Phase 1
    ## For each row of the prioritisation matrix satisfy:
    ### All first dose coverage >= prioritisation matrix target
    ### All second dose coverage for age groups prioritised for second doses >= prioritisation matrix target
    if(phase == 1){
      # Individuals left to target in each age group
      tp <- target_pop(dose_number = 1, dose_times = dose_times, prioritisation = prioritisation_matrix[p_step, ], t = t, dose_period = dose_period, d2_prioritise = d2_prioritise)
      # Assign available doses
      ad <- assign_doses(current_doses, tp)
      # Administer available doses
      dose_times <- administer_first_dose(ad, dose_times, t)
      # Recalculate remaining doses
      current_doses <- current_doses - sum(ad)

      if(current_doses > 0){
        # Individuals left to target in each age group
        tp <- target_pop(dose_number = 2, dose_times = dose_times, prioritisation = prioritisation_matrix[p_step, ], t = t, dose_period = dose_period, d2_prioritise = d2_prioritise)
        # Assign available doses
        ad <- assign_doses(current_doses, tp)
        # Administer available doses
        dose_times <- administer_second_dose(ad, dose_times, dose_period, t)
        # Recalculate remaining doses
        current_doses <- current_doses - sum(ad)
      }

      # Move to next row of the prioritisation matrix?
      dose_1_targets_met <- all(coverage(dose_times, 1) >= prioritisation_matrix[p_step, ])
      dose_2_prioritised_targets_met <-  all(coverage(dose_times, 2)[d2_prioritise] >= prioritisation_matrix[p_step, d2_prioritise])
      if(dose_1_targets_met & dose_2_prioritised_targets_met){
        p_step <- p_step + 1
        # If we have achieved all goals of phase 1 we set up for phase 2:
        if(p_step > nrow(prioritisation_matrix)){
          phase <- 2
          p_step <- 1
          d2_prioritise <- rep(TRUE, nrow(pop))
        }
      }
    }
    # Phase 2
    ## For each row of the prioritisation matrix satisfy:
    ### All second dose coverage >= prioritisation matrix target
    if(phase == 2) {
      # Individuals left to target in each age group
      tp <- target_pop(dose_number = 2, dose_times = dose_times, prioritisation = prioritisation_matrix[p_step, ], t = t, dose_period = dose_period, d2_prioritise = d2_prioritise)
      # Assign available doses
      ad <- assign_doses(current_doses, tp)
      # Administer available doses
      dose_times <- administer_second_dose(ad, dose_times, dose_period, t)
      # Recalculate remaining doses
      current_doses <- current_doses - sum(ad)

      # Move to next row of the prioritisation matrix?
      dose_2_targets_met <- all(coverage(dose_times, 2) >= prioritisation_matrix[p_step, ])
      if(dose_2_targets_met){
        p_step <- p_step + 1
        if(p_step > nrow(prioritisation_matrix)) break
      }
    }
  }
  # Extract outputs of interest
  output <- extract_dose_number(dose_times, maxt) %>%
    add_weighted_efficacy(infection_efficacy, disease_efficacy, v1v2)
  return(output)
}
