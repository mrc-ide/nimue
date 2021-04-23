#' Return the default probabilities for modelling defined in \code{squire}
#' For more info see \href{squire parameters vignette}{https://mrc-ide.github.io/squire/articles/parameters.html}
#' @return list of default probabilities
default_probs <- function() {
  c(squire::default_probs(), list(rel_infectiousness = rep(1, 17)))
}

probs <- default_probs()

#' Return the default hospital durations for modelling defined in \code{squire}
#' For more info see \href{squire parameters vignette}{https://mrc-ide.github.io/squire/articles/parameters.html}
#' @return list of default durations
default_durations <- function() {
  squire::default_durations()
}

durs <- default_durations()

#' Return the default vaccine parameters for modelling
#' @return list of default vaccine parameters
default_vaccine_pars <- function() {
  list(dur_R = Inf,
       dur_V = 365,
       vaccine_efficacy_infection = rep(0.95, 17),
       tt_vaccine_efficacy_infection = 0,
       vaccine_efficacy_disease = rep(0.95, 17),
       tt_vaccine_efficacy_disease = 0,
       max_vaccine = 1000,
       tt_vaccine = 0,
       dur_vaccine_delay = 14,
       vaccine_coverage_mat = matrix(0.8, ncol = 17, nrow = 1))
}

vaccine_pars <- default_vaccine_pars()

#' Vaccine parameters
#'
#' @details All durations are in days.
#'
#' @inheritParams run
parameters <- function(

  # Demography
  country = NULL,
  population = NULL,
  tt_contact_matrix = 0,
  contact_matrix_set = NULL,

  # Transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  # Initial state, duration, reps
  time_period = 365,
  seeding_cases,
  seeding_age_order = NULL,
  init = NULL,

  # Parameters
  # Probabilities
  prob_hosp = probs$prob_hosp,
  prob_severe = probs$prob_severe,
  prob_non_severe_death_treatment = probs$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = probs$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = probs$prob_severe_death_treatment,
  prob_severe_death_no_treatment = probs$prob_severe_death_no_treatment,
  p_dist = probs$p_dist,

  rel_infectiousness = probs$rel_infectiousness,

  # Durations
  dur_E,
  dur_IMild,
  dur_ICase,

  dur_get_ox_survive,
  dur_get_ox_die,
  dur_not_get_ox_survive,
  dur_not_get_ox_die,

  dur_get_mv_survive,
  dur_get_mv_die,
  dur_not_get_mv_survive,
  dur_not_get_mv_die,

  dur_rec,
  dur_R,

  # Vaccine
  dur_V,
  vaccine_efficacy_infection,
  tt_vaccine_efficacy_infection,
  vaccine_efficacy_disease,
  tt_vaccine_efficacy_disease,
  max_vaccine,
  tt_vaccine,
  dur_vaccine_delay,
  vaccine_coverage_mat,

  # Health system capacity
  hosp_bed_capacity,
  ICU_bed_capacity,
  tt_hosp_beds,
  tt_ICU_beds


) {

  # Handle country population args
  cpm <- squire:::parse_country_population_mixing_matrix(country = country,
                                                         population = population,
                                                         contact_matrix_set = contact_matrix_set)
  country <- cpm$country
  population <- cpm$population
  contact_matrix_set <- cpm$contact_matrix_set

  # Standardise contact matrix set
  if(is.matrix(contact_matrix_set)){
    contact_matrix_set <- list(contact_matrix_set)
  }

  # populate contact matrix set if not provided
  if (length(contact_matrix_set) == 1) {
    baseline <- contact_matrix_set[[1]]
    contact_matrix_set <- vector("list", length(tt_contact_matrix))
    for(i in seq_along(tt_contact_matrix)) {
      contact_matrix_set[[i]] <- baseline
    }
  }


  # populate hospital and ICU bed capacity if not provided
  if (is.null(hosp_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      hosp_beds <- beds$hosp_beds
      hosp_bed_capacity <- rep(round(hosp_beds * sum(population)/1000), length(tt_hosp_beds))
    } else {
      hosp_bed_capacity <- round(5 * sum(population)/1000)
    }
  }
  if (is.null(ICU_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      ICU_beds <- beds$ICU_beds
      ICU_bed_capacity <- rep(round(ICU_beds * sum(population)/1000), length(tt_ICU_beds))
    } else {
      ICU_bed_capacity <- round(3 * hosp_bed_capacity/100)
    }
  }

  # Initial state and matrix formatting
  # ----------------------------------------------------------------------------

  # Initialise initial conditions
  mod_init <- init(population, seeding_cases, seeding_age_order, init)

  # Convert contact matrices to input matrices
  matrices_set <- squire:::matrix_set_explicit(contact_matrix_set, population)

  # If a vector is put in for matrix targeting
  if(is.vector(vaccine_coverage_mat)){
    vaccine_coverage_mat <- matrix(vaccine_coverage_mat, ncol = 17)
  }

  # Input checks
  # ----------------------------------------------------------------------------
  mc <- squire:::matrix_check(population[-1], contact_matrix_set)
  stopifnot(length(R0) == length(tt_R0))
  stopifnot(length(contact_matrix_set) == length(tt_contact_matrix))
  stopifnot(length(hosp_bed_capacity) == length(tt_hosp_beds))
  stopifnot(length(ICU_bed_capacity) == length(tt_ICU_beds))
  stopifnot(length(max_vaccine) == length(tt_vaccine))
  stopifnot(ncol(vaccine_coverage_mat) == 17)

  assert_pos(dur_E)
  assert_pos(dur_IMild)
  assert_pos(dur_ICase)
  assert_pos(dur_get_ox_survive)
  assert_pos(dur_get_ox_die)
  assert_pos(dur_not_get_ox_survive)
  assert_pos(dur_not_get_ox_die)
  assert_pos(dur_get_mv_survive)
  assert_pos(dur_get_mv_die)
  assert_pos(dur_not_get_mv_survive)
  assert_pos(dur_not_get_mv_die)
  assert_pos(dur_R)
  assert_pos(dur_V)
  assert_pos(time_period)
  assert_pos(hosp_bed_capacity)
  assert_pos(ICU_bed_capacity)
  assert_pos(max_vaccine)
  assert_pos(dur_vaccine_delay)

  assert_length(prob_hosp, length(population))
  assert_length(prob_severe, length(population))
  assert_length(prob_non_severe_death_treatment, length(population))
  assert_length(prob_non_severe_death_no_treatment, length(population))
  assert_length(prob_severe_death_treatment, length(population))
  assert_length(prob_severe_death_no_treatment, length(population))
  assert_length(rel_infectiousness, length(population))
  assert_length(p_dist, length(population))

  assert_numeric(prob_hosp, length(population))
  assert_numeric(prob_severe, length(population))
  assert_numeric(prob_non_severe_death_treatment, length(population))
  assert_numeric(prob_non_severe_death_no_treatment, length(population))
  assert_numeric(prob_severe_death_treatment, length(population))
  assert_numeric(prob_severe_death_no_treatment, length(population))
  assert_numeric(rel_infectiousness, length(population))
  assert_numeric(p_dist, length(population))

  assert_leq(prob_hosp, 1)
  assert_leq(prob_severe, 1)
  assert_leq(prob_non_severe_death_treatment, 1)
  assert_leq(prob_non_severe_death_no_treatment, 1)
  assert_leq(prob_severe_death_treatment, 1)
  assert_leq(prob_severe_death_no_treatment, 1)
  assert_leq(rel_infectiousness, 1)
  assert_leq(p_dist, 1)

  assert_greq(prob_hosp, 0)
  assert_greq(prob_severe, 0)
  assert_greq(prob_non_severe_death_treatment, 0)
  assert_greq(prob_non_severe_death_no_treatment, 0)
  assert_greq(prob_severe_death_treatment, 0)
  assert_greq(prob_severe_death_no_treatment, 0)
  assert_greq(rel_infectiousness, 0)
  assert_greq(p_dist, 0)


  # Convert and Generate Parameters As Required
  # ----------------------------------------------------------------------------

  # durations
  gamma_E = 2 * 1/dur_E
  gamma_IMild = 1/dur_IMild
  gamma_ICase = 2 * 1/dur_ICase
  gamma_get_ox_survive = 2 * 1/dur_get_ox_survive
  gamma_get_ox_die = 2 * 1/dur_get_ox_die
  gamma_not_get_ox_survive = 2 * 1/dur_not_get_ox_survive
  gamma_not_get_ox_die = 2 * 1/dur_not_get_ox_die
  gamma_get_mv_survive = 2 * 1/dur_get_mv_survive
  gamma_get_mv_die = 2 * 1/dur_get_mv_die
  gamma_not_get_mv_survive = 2 * 1/dur_not_get_mv_survive
  gamma_not_get_mv_die = 2 * 1/dur_not_get_mv_die
  gamma_rec = 2 * 1/dur_rec
  gamma_R <- 2 * 1/dur_R
  gamma_V <- 2 * 1/dur_V
  gamma_vaccine_delay <- 2 * 1 / dur_vaccine_delay

  if (is.null(beta_set)) {
    baseline_matrix <- squire:::process_contact_matrix_scaled_age(contact_matrix_set[[1]], population)
    beta_set <- beta_est_infectiousness(dur_IMild = dur_IMild,
                                          dur_ICase = dur_ICase,
                                          prob_hosp = prob_hosp,
                                          mixing_matrix = baseline_matrix,
                                          rel_infectiousness = rel_infectiousness,
                                          R0 = R0)
  }

  # normalise to sum to 1
  p_dist <- matrix(rep(p_dist, 6), nrow = 17, ncol = 6)
  p_dist <- p_dist/mean(p_dist)

  # Format vaccine-specific parameters
  gamma_vaccine <- c(0, gamma_vaccine_delay, gamma_vaccine_delay, gamma_V, gamma_V, 0)

  # Vaccine efficacies are now time changing (if specified),
  # so we need to convert these to be interpolated by odin
  # These functions also check that efficacies are correct length
  # both in terms of age groups and in terms of required timepoints

  # First the vaccine efficacy infection
  vaccine_efficacy_infection_odin_array <- format_ve_i_for_odin(
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection
    )

  # Second the vaccine efficacy disease affecting prob_hosp
  prob_hosp_odin_array <- format_ve_d_for_odin(
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
    prob_hosp = prob_hosp
  )

  # Collate Parameters Into List
  pars <- c(mod_init,
            list(N_age = length(population),
                 gamma_E = gamma_E,
                 gamma_IMild = gamma_IMild,
                 gamma_ICase = gamma_ICase,
                 gamma_get_ox_survive = gamma_get_ox_survive,
                 gamma_get_ox_die = gamma_get_ox_die,
                 gamma_not_get_ox_survive = gamma_not_get_ox_survive,
                 gamma_not_get_ox_die = gamma_not_get_ox_die,
                 gamma_get_mv_survive = gamma_get_mv_survive,
                 gamma_get_mv_die = gamma_get_mv_die,
                 gamma_not_get_mv_survive = gamma_not_get_mv_survive,
                 gamma_not_get_mv_die = gamma_not_get_mv_die,
                 gamma_rec = gamma_rec,
                 gamma_R = gamma_R,
                 prob_hosp = prob_hosp_odin_array,
                 prob_severe = prob_severe,
                 prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                 prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                 prob_severe_death_treatment = prob_severe_death_treatment,
                 prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                 rel_infectiousness = rel_infectiousness,
                 p_dist = p_dist,
                 hosp_beds = hosp_bed_capacity,
                 ICU_beds = ICU_bed_capacity,
                 tt_hosp_beds = tt_hosp_beds,
                 tt_ICU_beds = tt_ICU_beds,
                 tt_matrix = tt_contact_matrix,
                 mix_mat_set = matrices_set,
                 tt_beta = tt_R0,
                 beta_set = beta_set,
                 population = population,
                 contact_matrix_set = contact_matrix_set,
                 max_vaccine = max_vaccine,
                 tt_vaccine = tt_vaccine,
                 vaccine_efficacy_infection = vaccine_efficacy_infection_odin_array,
                 tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                 tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                 vaccine_coverage_mat = vaccine_coverage_mat,
                 N_vaccine = 6,
                 N_prioritisation_steps = nrow(vaccine_coverage_mat),
                 gamma_vaccine = gamma_vaccine))

  class(pars) <- c("vaccine_parameters", "nimue_parameters")

  return(pars)
}

#' Estimate beta parameter for explicit model
#'
#' @param dur_IMild Duration of mild infectiousness (days)
#' @param dur_ICase Delay between symptom onset and requiring hospitalisation (days)
#' @param prob_hosp Probability of hospitilisation by ages
#' @param rel_infectiousness Relative infectiousness of age categories relative
#'   to maximum infectiousness age category
#' @param mixing_matrix Mixing matrix
#' @param R0 Basic reproduction number
#'
#' @return Beta parameter
#' @export
#'
# #' @examples
beta_est_infectiousness <- function(dur_IMild,
                                    dur_ICase,
                                    prob_hosp,
                                    rel_infectiousness,
                                    mixing_matrix,
                                    R0) {

  # assertions
  assert_single_pos(dur_ICase, zero_allowed = FALSE)
  assert_single_pos(dur_IMild, zero_allowed = FALSE)
  assert_numeric(prob_hosp)
  assert_numeric(rel_infectiousness)
  assert_same_length(prob_hosp, rel_infectiousness)
  assert_numeric(mixing_matrix)
  assert_square_matrix(mixing_matrix)
  assert_same_length(mixing_matrix[,1], prob_hosp)
  assert_pos(R0, zero_allowed = FALSE)

  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }

  if(sum(is.na(rel_infectiousness)) > 0) {
    stop("rel_infectiousness must not contain NAs")
  }

  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }

  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  adjusted_eigen <- Re(eigen(mixing_matrix*relative_R0_by_age*rel_infectiousness)$values[1])
  R0 / adjusted_eigen

}

#' @noRd
format_ve_i_for_odin <- function(vaccine_efficacy_infection,
                                 tt_vaccine_efficacy_infection) {

# If just provided as a vector then we put into a list ready for formatting
if(!is.list(vaccine_efficacy_infection)){
  vaccine_efficacy_infection <- list(vaccine_efficacy_infection)
}

# check that the correct length agreement between tt_vaccine_efficacy_infection
assert_length(vaccine_efficacy_infection, length(tt_vaccine_efficacy_infection))

# now check that each vaccine efficacy is correct length (1 or 17)
vaccine_efficacy_infection <- lapply(vaccine_efficacy_infection, function(ve_i) {

  if(length(ve_i) == 1){
    ve_i <- rep(ve_i, 17)
  }

  if(length(ve_i) != 17){
    stop("Parameter vaccine_efficacy_infection must be length 1 or length 17")
  }

  return(ve_i)

})

# and now format so each list is the vaccine_efficacy_infection at each time
# point for the 6 vaccine classes
ve_i_list <- lapply(vaccine_efficacy_infection, function(ve_i) {

  ve_i = 1 - ve_i
  return(matrix(c(rep(1, 17 * 3),
                  ve_i, ve_i,
                  rep(1, 17)), nrow = 17, ncol = 6))

})

# and use this list to create an array that is in right format for odin
vaccine_efficacy_infection_odin_array <- aperm(
  array(unlist(ve_i_list), dim = c(dim(ve_i_list[[1]]), length(ve_i_list))),
  c(3, 1, 2)
)

return(vaccine_efficacy_infection_odin_array)

}


#' @noRd
format_ve_d_for_odin <- function(vaccine_efficacy_disease,
                                 tt_vaccine_efficacy_disease,
                                 prob_hosp) {


  # If just provided as a vector then we put into a list ready for formatting
  if(!is.list(vaccine_efficacy_disease)){
    vaccine_efficacy_disease <- list(vaccine_efficacy_disease)
  }

  # check that the correct length agreement between tt_vaccine_efficacy_disease
  assert_length(vaccine_efficacy_disease, length(tt_vaccine_efficacy_disease))

  # now check that each vaccine efficacy is correct length (1 or 17)
  vaccine_efficacy_disease <- lapply(vaccine_efficacy_disease, function(ve_d) {

    if(length(ve_d) == 1){
      ve_d <- rep(ve_d, 17)
    }

    if(length(ve_d) != 17){
      stop("Parameter vaccine_efficacy_disease must be length 1 or length 17")
    }

    return(ve_d)

  })

  # and now format so each list is the prob_hosp at each time
  # point for the 6 vaccine classes
  prob_hosp_list <- lapply(vaccine_efficacy_disease, function(ve_d) {

    prob_hosp_vaccine = (1 - ve_d) * prob_hosp

    # age X vaccine efficacy parameters
    prob_hosp <- matrix(c(prob_hosp, prob_hosp, prob_hosp,
                          prob_hosp_vaccine, prob_hosp_vaccine,
                          prob_hosp), nrow = 17, ncol = 6)

    return(prob_hosp)

  })

  # and use this list to create an array that is in right format for odin
  prob_hosp_odin_array <- aperm(
    array(unlist(prob_hosp_list), dim = c(dim(prob_hosp_list[[1]]), length(prob_hosp_list))),
    c(3, 1, 2)
  )

  return(prob_hosp_odin_array)

}
