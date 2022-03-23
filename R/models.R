#' Create a nimue model class for fitting with squire tools
#'
#' @title nimue model creation.
#' @param use_dde Logical for using dde to solve. Default = TRUE
#' We will use this structure to ensure that model fitting is flexible in the
#' future as more models are added
#'
#' @details Wraps the squire pmcmc fitting infrastructure.
#'
#' @export
nimue_deterministic_model <- function(use_dde = TRUE) {

  model_class <- "nimue_model"

  compare_model <- function(model, pars_obs, data) {
    squire:::compare_output(model, pars_obs, data, type=model_class)
  }

  # wrap param func in order to remove unused arguments (dt)
  # and then add in all the default that are passed to params usually
  # from run so have to add here
  parameters_func <- function(country, population, dt,
                              contact_matrix_set, tt_contact_matrix,
                              hosp_bed_capacity, tt_hosp_beds,
                              ICU_bed_capacity, tt_ICU_beds,

                              # vaccine defaults that are just empty in parms so declare here
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

                              # durations
                              dur_E  = durs$dur_E,
                              dur_IMild = durs$dur_IMild,
                              dur_ICase = durs$dur_ICase,

                              # hospital durations
                              dur_get_ox_survive = durs$dur_get_ox_survive,
                              tt_dur_get_ox_survive = durs$tt_dur_get_ox_survive,
                              dur_get_ox_die = durs$dur_get_ox_die,
                              tt_dur_get_ox_die = durs$tt_dur_get_ox_die,
                              dur_not_get_ox_survive = durs$dur_not_get_ox_survive,
                              dur_not_get_ox_die = durs$dur_not_get_ox_die,

                              dur_get_mv_survive = durs$dur_get_mv_survive,
                              tt_dur_get_mv_survive = durs$tt_dur_get_mv_survive,
                              dur_get_mv_die = durs$dur_get_mv_die,
                              tt_dur_get_mv_die = durs$tt_dur_get_mv_die,
                              dur_not_get_mv_survive = durs$dur_not_get_mv_survive,
                              dur_not_get_mv_die = durs$dur_not_get_mv_die,

                              dur_rec = durs$dur_rec,

                              # seeding cases default
                              seeding_cases = 5,

                              ...) {

    pars <- parameters(
      country = country,
               population = population,
               contact_matrix_set = contact_matrix_set,
               tt_contact_matrix = tt_contact_matrix,
               hosp_bed_capacity = hosp_bed_capacity,
               tt_hosp_beds = tt_hosp_beds,
               ICU_bed_capacity = ICU_bed_capacity,
               tt_ICU_beds = tt_ICU_beds,
               dur_E = dur_E,
               dur_IMild = dur_IMild,
               dur_ICase = dur_ICase,
               dur_get_ox_survive = dur_get_ox_survive,
               tt_dur_get_ox_survive = tt_dur_get_ox_survive,
               dur_get_ox_die = dur_get_ox_die,
               tt_dur_get_ox_die = tt_dur_get_ox_die,
               dur_not_get_ox_survive = dur_not_get_ox_survive,
               dur_not_get_ox_die = dur_not_get_ox_die,
               dur_get_mv_survive = dur_get_mv_survive,
               tt_dur_get_mv_survive = tt_dur_get_mv_survive,
               dur_get_mv_die = dur_get_mv_die,
               tt_dur_get_mv_die = tt_dur_get_mv_die,
               dur_not_get_mv_survive = dur_not_get_mv_survive,
               dur_not_get_mv_die = dur_not_get_mv_die,
               dur_rec = dur_rec,
               dur_R = dur_R,
               tt_dur_R = tt_dur_R,
               dur_V = dur_V,
               vaccine_efficacy_infection = vaccine_efficacy_infection,
               tt_vaccine_efficacy_infection = tt_vaccine_efficacy_infection,
               vaccine_efficacy_disease = vaccine_efficacy_disease,
               tt_vaccine_efficacy_disease = tt_vaccine_efficacy_disease,
               max_vaccine = max_vaccine,
               tt_vaccine = tt_vaccine ,
               dur_vaccine_delay = dur_vaccine_delay,
               vaccine_coverage_mat = vaccine_coverage_mat,
               seeding_cases = seeding_cases,
               ...)

    # append extra pars for fitting
    pars$dt <- dt
    pars$prob_hosp_baseline <- pars$prob_hosp[1, ,1]
    pars$use_dde <- use_dde

    class(pars) <- c("vaccine_parameters", "squire_parameters")
    return(pars)

  }

  # wrap run func correctly
  run_func <- function(country, population, dt,
                       contact_matrix_set, tt_contact_matrix,
                       hosp_bed_capacity, tt_hosp_beds,
                       ICU_bed_capacity, tt_ICU_beds,
                       replicates = 1,
                       day_return = TRUE,
                       time_period = 365,
                       ...) {

    out <- nimue::run(country = country,
               contact_matrix_set = contact_matrix_set,
               tt_contact_matrix = tt_contact_matrix,
               hosp_bed_capacity = hosp_bed_capacity,
               tt_hosp_beds = tt_hosp_beds,
               ICU_bed_capacity = ICU_bed_capacity,
               tt_ICU_beds = tt_ICU_beds,
               population = population,
               replicates = 1,
               time_period = time_period,
               use_dde = use_dde,
               ...)

    return(out)

  }

  odin_model <- function(user, unused_user_action) {
    vaccine$new(user = user, use_dde = use_dde, unused_user_action = "ignore")
  }

  model <- list(odin_model = odin_model,
                generate_beta_func = beta_est_infectiousness,
                parameter_func = parameters_func,
                run_func = run_func,
                compare_model = compare_model,
                use_dde = use_dde)
  class(model) <- c(model_class, "deterministic", "squire_model")
  model

}

