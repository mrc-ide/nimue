#' Run the vaccine model
#'
#' @param population Population vector (for each age group). Default = NULL,
#'   which will cause population to be sourced from \code{country}
#' @param country Character for country beign simulated. WIll be used to
#'   generate \code{population} and \code{contact_matrix_set} if
#'   unprovided. Either \code{country} or \code{population} and
#'   \code{contact_matrix_set} must be provided.
#' @param contact_matrix_set Contact matrices used in simulation. Default =
#'   NULL, which will generate this based on the \code{country}.
#' @param tt_contact_matrix Time change points for matrix change. Default = 0
#' @param R0 Basic Reproduction Number. Default = 3
#' @param tt_R0 Change time points for R0. Default = 0
#' @param beta_set Alternative parameterisation via beta rather than R0.
#'   Default = NULL, which causes beta to be estimated from R0
#' @param time_period Length of simulation. Default = 365
#' @param replicates  Number of replicates. Default = 10
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param seed Random seed used for simulations. Deafult = runif(1, 0, 10000)
#' @param prob_hosp probability of hospitalisation by age.
#'   Default = c(0.000744192, 0.000634166,0.001171109, 0.002394593, 0.005346437,
#'   0.010289885, 0.016234604, 0.023349169, 0.028944623, 0.038607042,
#'   0.057734879, 0.072422135, 0.101602458, 0.116979814, 0.146099064,
#'   0.176634654 ,0.180000000)
#' @param prob_severe Probability of developing severe symptoms by age.
#'   Default = c(0.05022296,	0.05022296,	0.05022296,	0.05022296,	0.05022296,
#'   0.05022296,	0.05022296,	0.053214942, 0.05974426,	0.074602879,
#'   0.103612417, 0.149427991, 0.223777304,	0.306985918,
#'   0.385779555, 0.461217861, 0.709444444)
#' @param prob_non_severe_death_treatment Probability of death from non severe
#'   treated infection.
#'   Default = c(0.0125702,	0.0125702,	0.0125702,	0.0125702,
#'   0.0125702,	0.0125702,	0.0125702,	0.013361147,
#'   0.015104687,	0.019164124,	0.027477519,	0.041762108,
#'   0.068531658,	0.105302319,	0.149305732,	0.20349534,	0.5804312)
#' @param prob_severe_death_treatment Probability of death from severe infection
#'   that is treated. Default = rep(0.5, 17)
#' @param prob_non_severe_death_no_treatment Probability of death in non severe
#'   hospital inections that aren't treated
#' @param prob_severe_death_no_treatment Probability of death from severe infection
#'   that is not treated. Default = rep(0.95, 17)
#' @param p_dist Preferentiality of age group receiving treatment relative to
#'   other age groups when demand exceeds healthcare capacity.
#' @param rel_infectiousness Relative infectiousness per age category relative
#'   to maximum infectiousness category. Default = rep(1, 17)
#' @param rel_infectiousness_vaccinated  Relative infectiousness per age
#'   category  of vaccinated individuals relative to unvaccinated individuals.
#'   Default = rep(1, 17), which is no impact of vaccination on onwards
#'   transmissions
#' @param dur_E Mean duration of incubation period (days). Default = 4.6
#' @param dur_IMild Mean duration of mild infection (days). Default = 2.1
#' @param dur_ICase Mean duration from symptom onset to hospitil admission (days).
#'   Default = 4.5
#' @param dur_get_ox_survive Mean duration of oxygen given survive. Default = 5
#' @param dur_get_ox_die Mean duration of oxygen given death. Default = 5
#' @param dur_not_get_ox_survive Mean duration without oxygen given survive.
#'   Default = 5
#' @param dur_not_get_ox_die Mean duration without  oxygen given death.
#'  Default = 5
#' @param dur_get_mv_survive Mean duration of ventilation given survive.
#'   Default = 7.3
#' @param dur_get_mv_die Mean duration of ventilation given death. Default = 6
#' @param dur_not_get_mv_survive Mean duration without ventilation given
#'   survive. Default = 7.3
#' @param dur_not_get_mv_die Mean duration without ventilation given
#'   death. Default = 1
#' @param dur_rec Duration of recovery after coming off ventilation. Default = 2
#' @param hosp_bed_capacity General bed capacity. Can be single number of vector if capacity time-varies.
#' @param ICU_bed_capacity ICU bed capacity. Can be single number of vector if capacity time-varies.
#' @param tt_hosp_beds Times at which hospital bed capacity changes (Default = 0 = doesn't change)
#' @param tt_ICU_beds Times at which ICU bed capacity changes (Default = 0 = doesn't change)
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param seeding_age_order Vector specifying the order in which seeds are allocated to ages.
#'   If NULL, seeds are distributed randomly within working ages. If specified, must be a vector
#'   of length 17 specifying the order seeds are allocated, e.g. 1:17 will allocate first seed
#'   to the youngest age group, then the second youngest and so on. Default = NULL
#' @param init Initial conditions for simulation provided. Allows overriding
#'   if initial conditions start with an already infected population etc.
#'   Default = NULL.
#' @param dur_R Mean duration of naturally acquired immunity (days)
#' @param dur_V Mean duration of vaccine-derived immunity (days)
#' @param vaccine_efficacy_infection Efficacy of vaccine against infection.
#'   This parameter must either be length 1 numeric (a single efficacy for all
#'   age groups) or length 17 numeric vector (an efficacy for each age group).
#'   An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will
#'   reduce FOI by 20 percent etc.
#'   To specify changes in vaccine efficacy over time, vaccine efficacies must
#'   be provided as a list, with each list element being the efficacy at each
#'   time point specified by \code{tt_vaccine_efficacy_infection}. These
#'   efficacies must also be length 1 numeric (a single efficacy for all age
#'   groups)  or length 17 numeric vector (an efficacy for each age group)
#' @param tt_vaccine_efficacy_infection Timing of changes in vaccine efficacy
#'   against infection. Default = 0, which assumes fixed efficacy over time.
#'   Must be the same length as the length of \code{vaccine_efficacy_infection}
#'   when provided as a list. Time changing efficacies can occur in response to
#'   changing vaccines being  given and dosing strategy changes.
#' @param vaccine_efficacy_disease Efficacy of vaccine against severe
#'   (requiring hospitilisation) disease (by age). This parameter must either be
#'   length 1 numeric (a single efficacy for all age groups) or length 17
#'   numeric vector (an efficacy for each age group). An efficacy of 1 will
#'   reduce the probability of hospitalisation by 100 percent, an efficacy of
#'   0.2 will reduce the probability of hospitalisation by 20 percent etc.
#'   To specify changes in vaccine efficacy over time, vaccine efficacies must
#'   be provided as a list, with each list element being the efficacy at each
#'   time point specified by \code{tt_vaccine_efficacy_disease}. These
#'   efficacies must also be length 1 numeric (a single efficacy for all age
#'   groups)  or length 17 numeric vector (an efficacy for each age group).
#' @param tt_vaccine_efficacy_disease Timing of changes in vaccine efficacy
#'   against disease. Default = 0, which assumes fixed efficacy over time.
#'   Must be the same length as the length of \code{vaccine_efficacy_disease}
#'   when provided as a list. Time changing efficacies can occur in response to
#'   changing vaccines being  given and dosing strategy changes.
#' @param max_vaccine The maximum number of individuals who can be vaccinated per day.
#' @param tt_vaccine Time change points for vaccine capacity (\code{max_vaccine}).
#' @param dur_vaccine_delay Mean duration of period from vaccination to vaccine protection.
#' @param vaccine_coverage_mat Vaccine coverage targets by age (columns) and priority (row)
#' @param use_dde Use the dde solver (default is \code{TRUE})
#' @param ... Additional arguments for solver
#'
#' @return Simulation output
#' @export
run <- function(

  # demography
  country = NULL,
  population = NULL,
  tt_contact_matrix = 0,
  contact_matrix_set = NULL,

  # transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  # initial state, duration, reps
  time_period = 365,
  replicates = 10,
  seed = stats::runif(1, 0, 100000000),

  # parameters
  # probabilities
  prob_hosp = probs$prob_hosp,
  prob_hosp_multiplier = probs$prob_hosp_multiplier,
  tt_prob_hosp_multiplier = probs$tt_prob_hosp_multiplier,
  prob_severe = probs$prob_severe,
  prob_non_severe_death_treatment = probs$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = probs$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = probs$prob_severe_death_treatment,
  prob_severe_death_no_treatment = probs$prob_severe_death_no_treatment,
  p_dist = probs$p_dist,

  # onward infectiousness
  rel_infectiousness = probs$rel_infectiousness,
  rel_infectiousness_vaccinated = probs$rel_infectiousness_vaccinated,

  # durations
  dur_E  = durs$dur_E,
  dur_IMild = durs$dur_IMild,
  dur_ICase = durs$dur_ICase,

  # hospital durations
  dur_get_ox_survive = durs$dur_get_ox_survive,
  dur_get_ox_die = durs$dur_get_ox_die,
  dur_not_get_ox_survive = durs$dur_not_get_ox_survive,
  dur_not_get_ox_die = durs$dur_not_get_ox_die,

  dur_get_mv_survive = durs$dur_get_mv_survive,
  dur_get_mv_die = durs$dur_get_mv_die,
  dur_not_get_mv_survive = durs$dur_not_get_mv_survive,
  dur_not_get_mv_die = durs$dur_not_get_mv_die,

  dur_rec = durs$dur_rec,

  # vaccine
  dur_R = vaccine_pars$dur_R,
  tt_dur_R = vaccine_pars$tt_dur_R,
  dur_V = vaccine_pars$dur_V,
  vaccine_efficacy_infection = vaccine_pars$vaccine_efficacy_infection,
  tt_vaccine_efficacy_infection = vaccine_pars$tt_vaccine_efficacy_infection,
  vaccine_efficacy_disease = vaccine_pars$vaccine_efficacy_disease,
  tt_vaccine_efficacy_disease = vaccine_pars$tt_vaccine_efficacy_disease,
  max_vaccine = vaccine_pars$max_vaccine,
  tt_vaccine = vaccine_pars$tt_vaccine,
  dur_vaccine_delay = vaccine_pars$dur_vaccine_delay,
  vaccine_coverage_mat = vaccine_pars$vaccine_coverage_mat,

  # health system capacity
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  tt_hosp_beds = 0,
  tt_ICU_beds = 0,

  seeding_cases = 20,
  seeding_age_order = NULL,
  init = NULL,
  use_dde = TRUE,
  ...
) {

  # Grab function arguments
  args <- as.list(environment())
  set.seed(seed)

  # create parameter list
  pars <- parameters(country = country,
                     population = population,
                     tt_contact_matrix = tt_contact_matrix ,
                     contact_matrix_set = contact_matrix_set,
                     R0 = R0,
                     tt_R0 = tt_R0 ,
                     beta_set = beta_set,
                     time_period = time_period,
                     seeding_cases = seeding_cases,
                     seeding_age_order = seeding_age_order,
                     prob_hosp = prob_hosp,
                     tt_prob_hosp_multiplier = tt_prob_hosp_multiplier,
                     prob_hosp_multiplier = prob_hosp_multiplier,
                     prob_severe = prob_severe,
                     prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                     prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                     prob_severe_death_treatment = prob_severe_death_treatment,
                     prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                     p_dist = p_dist,
                     rel_infectiousness = rel_infectiousness,
                     rel_infectiousness_vaccinated = rel_infectiousness_vaccinated,
                     dur_E = dur_E,
                     dur_IMild = dur_IMild,
                     dur_ICase = dur_ICase,
                     dur_get_ox_survive = dur_get_ox_survive,
                     dur_get_ox_die = dur_get_ox_die,
                     dur_not_get_ox_survive = dur_not_get_ox_survive,
                     dur_not_get_ox_die = dur_not_get_ox_die,
                     dur_get_mv_survive = dur_get_mv_survive,
                     dur_get_mv_die = dur_get_mv_die,
                     dur_not_get_mv_survive = dur_not_get_mv_survive,
                     dur_not_get_mv_die = dur_not_get_mv_die,
                     dur_rec = dur_rec,
                     dur_R = dur_R,
                     tt_dur_R = tt_dur_R,
                     hosp_bed_capacity = hosp_bed_capacity,
                     ICU_bed_capacity = ICU_bed_capacity,
                     tt_hosp_beds = tt_hosp_beds,
                     tt_ICU_beds = tt_ICU_beds,
                     dur_V = dur_V,
                     vaccine_efficacy_infection = vaccine_efficacy_infection,
                     tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
                     vaccine_efficacy_disease = vaccine_efficacy_disease,
                     tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
                     max_vaccine = max_vaccine,
                     tt_vaccine = tt_vaccine,
                     dur_vaccine_delay = dur_vaccine_delay,
                     vaccine_coverage_mat = vaccine_coverage_mat,
                     init = init)

  # Set model type
  replicates <- 1
  mod_gen = vaccine

  # Running the Model
  mod <- mod_gen(user = pars, unused_user_action = "ignore",
                 use_dde = use_dde)

  # Daily output by default
  t <- round(seq(from = 1, to = time_period))

  results <- mod$run(t, ...)

  # coerce to array
  results <- array(results, dim = c(dim(results), 1), dimnames = dimnames(results))

  # Summarise inputs
  parameters <- args
  parameters$population <- pars$population
  parameters$hosp_bed_capacity <- pars$hosp_beds
  parameters$ICU_bed_capacity <- pars$ICU_beds
  parameters$beta_set <- pars$beta_set
  parameters$seeding_cases <- pars$E1_0
  parameters$contact_matrix_set <- pars$contact_matrix_set

  out <- list(output = results, parameters = parameters, model = mod, odin_parameters = pars)
  out <- structure(out, class = "nimue_simulation")
  return(out)

}
