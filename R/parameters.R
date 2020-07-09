#' Return the default probabilities for modelling
#' @return list of default probabilities
default_probs <- function() {
  prob_hosp <- c(
    0.000744192, 0.000634166,0.001171109, 0.002394593, 0.005346437 ,
    0.010289885, 0.016234604, 0.023349169, 0.028944623, 0.038607042 ,
    0.057734879, 0.072422135, 0.101602458, 0.116979814, 0.146099064,
    0.176634654 ,0.180000000)
  list(
    prob_hosp = prob_hosp,
    prob_severe = c(
      0.05022296,	0.05022296,	0.05022296,	0.05022296,	0.05022296,
      0.05022296,	0.05022296,	0.053214942, 0.05974426,	0.074602879,
      0.103612417, 0.149427991, 0.223777304,	0.306985918,
      0.385779555, 0.461217861, 0.709444444),
    prob_non_severe_death_treatment = c(
      0.0125702,	0.0125702,	0.0125702,	0.0125702,
      0.0125702,	0.0125702,	0.0125702,	0.013361147,
      0.015104687,	0.019164124,	0.027477519,	0.041762108,
      0.068531658,	0.105302319,	0.149305732,	0.20349534,	0.5804312),
    prob_non_severe_death_no_treatment = rep(0.6, length(prob_hosp)),
    prob_severe_death_treatment = rep(0.5, length(prob_hosp)),
    prob_severe_death_no_treatment = rep(0.95, length(prob_hosp)),
    p_dist = rep(1, length(prob_hosp))
  )
}

probs <- default_probs()


#' Return the default vaccine parameters for modelling
#' @return list of default vaccine parameters
default_vaccine_pars <- function() {
  list(dur_R = Inf,
       vaccination_target = rep(1, 17),
       dur_V = 365,
       vaccine_efficacy_infection = rep(0.95, 17),
       vaccine_efficacy_disease = rep(0.95, 17),
       max_vaccine = 1000,
       tt_vaccine = 0,
       dur_vaccine_delay = 14)
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
  dt = 0.1,
  seeding_cases,

  # Parameters
  # Probabilities
  prob_hosp = probs$prob_hosp,
  prob_severe = probs$prob_severe,
  prob_non_severe_death_treatment = probs$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = probs$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = probs$prob_severe_death_treatment,
  prob_severe_death_no_treatment = probs$prob_severe_death_no_treatment,
  p_dist = probs$p_dist,

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
  vaccination_target,
  dur_V,
  vaccine_efficacy_infection,
  vaccine_efficacy_disease,
  max_vaccine,
  tt_vaccine,
  dur_vaccine_delay,

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
  mod_init <- init(population, seeding_cases)

  # Convert contact matrices to input matrices
  matrices_set <- squire:::matrix_set_explicit(contact_matrix_set, population)

  # Input checks
  # ----------------------------------------------------------------------------
  mc <- squire:::matrix_check(population[-1], contact_matrix_set)
  stopifnot(length(R0) == length(tt_R0))
  stopifnot(length(contact_matrix_set) == length(tt_contact_matrix))
  stopifnot(length(hosp_bed_capacity) == length(tt_hosp_beds))
  stopifnot(length(ICU_bed_capacity) == length(tt_ICU_beds))
  stopifnot(length(max_vaccine) == length(tt_vaccine))
  tc <- lapply(list(tt_R0/dt, tt_contact_matrix/dt), squire:::check_time_change, time_period/dt)
  tc2 <- lapply(list(tt_hosp_beds/dt, tt_ICU_beds/dt), squire:::check_time_change, time_period/dt)
  stopifnot(all(vaccination_target %in% 0:1))

  assert_pos(dt)
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
  assert_length(p_dist, length(population))

  assert_numeric(prob_hosp, length(population))
  assert_numeric(prob_severe, length(population))
  assert_numeric(prob_non_severe_death_treatment, length(population))
  assert_numeric(prob_non_severe_death_no_treatment, length(population))
  assert_numeric(prob_severe_death_treatment, length(population))
  assert_numeric(prob_severe_death_no_treatment, length(population))
  assert_numeric(p_dist, length(population))

  assert_leq(prob_hosp, 1)
  assert_leq(prob_severe, 1)
  assert_leq(prob_non_severe_death_treatment, 1)
  assert_leq(prob_non_severe_death_no_treatment, 1)
  assert_leq(prob_severe_death_treatment, 1)
  assert_leq(prob_severe_death_no_treatment, 1)
  assert_leq(p_dist, 1)

  assert_greq(prob_hosp, 0)
  assert_greq(prob_severe, 0)
  assert_greq(prob_non_severe_death_treatment, 0)
  assert_greq(prob_non_severe_death_no_treatment, 0)
  assert_greq(prob_severe_death_treatment, 0)
  assert_greq(prob_severe_death_no_treatment, 0)
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
    beta_set <- squire::beta_est_explicit(dur_IMild = dur_IMild,
                                          dur_ICase = dur_ICase,
                                          prob_hosp = prob_hosp,
                                          mixing_matrix = baseline_matrix,
                                          R0 = R0)
  }

  # normalise to sum to 1
  p_dist <- matrix(rep(p_dist, 6), nrow = 17, ncol = 6)
  p_dist <- p_dist/mean(p_dist)

  # Format vaccine-specific parameters
  vaccine_efficacy_infection = 1 - vaccine_efficacy_infection
  prob_hosp_vaccine = (1 - vaccine_efficacy_disease) * prob_hosp

  # age X vaccine efficacy parameters
  prob_hosp <- matrix(c(prob_hosp, prob_hosp, prob_hosp,
                        prob_hosp_vaccine, prob_hosp_vaccine,
                        prob_hosp), nrow = 17, ncol = 6)
  vaccine_efficacy_infection <- matrix(c(rep(1, 17 * 3),
                                         vaccine_efficacy_infection, vaccine_efficacy_infection,
                                         rep(1, 17)), nrow = 17, ncol = 6)
  gamma_vaccine <- c(0, gamma_vaccine_delay, gamma_vaccine_delay, gamma_V, gamma_V, 0)


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
                 prob_hosp = prob_hosp,
                 prob_severe = prob_severe,
                 prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                 prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                 prob_severe_death_treatment = prob_severe_death_treatment,
                 prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                 p_dist = p_dist,
                 hosp_beds = hosp_bed_capacity,
                 ICU_beds = ICU_bed_capacity,
                 tt_hosp_beds = tt_hosp_beds,
                 tt_ICU_beds = tt_ICU_beds,
                 tt_matrix = tt_contact_matrix,
                 mix_mat_set = matrices_set,
                 tt_beta = tt_R0,
                 beta_set = beta_set,
                 dt = dt,
                 population = population,
                 contact_matrix_set = contact_matrix_set,
                 vaccination_target = vaccination_target,
                 max_vaccine = max_vaccine,
                 tt_vaccine = tt_vaccine,
                 vaccine_efficacy_infection = vaccine_efficacy_infection,
                 N_vaccine = 6,
                 gamma_vaccine = gamma_vaccine))

  class(pars) <- c("vaccine_parameters", "nimue_parameters")

  return(pars)
}
