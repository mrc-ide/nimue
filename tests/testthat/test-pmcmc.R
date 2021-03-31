
#------------------------------------------------
test_that("pmcmc nimue", {

  Sys.setenv("SQUIRE_PARALLEL_DEBUG" = "TRUE")
  data <- read.csv(squire:::squire_file("extdata/example.csv"),stringsAsFactors = FALSE)
  interventions <- read.csv(squire:::squire_file("extdata/example_intervention.csv"))
  int_unique <- squire:::interventions_unique(interventions)
  reporting_fraction = 1
  country <- "Algeria"
  pars_init <- list('start_date'     = as.Date("2020-02-07"),
                    'R0'             = 2.5,
                    'Meff'           = 2)
  pars_min <- list('start_date'      = as.Date("2020-02-01"),
                   'R0'              = 1e-10,
                   'Meff'            = 0.1)
  pars_max <- list('start_date'      = as.Date("2020-02-20"),
                   'R0'              = 5,
                   'Meff'            = 5)
  pars_discrete <- list('start_date' = TRUE,
                        'R0'         = FALSE,
                        'Meff'       = FALSE)
  pars_obs <- list(phi_cases = 0.1,
                   k_cases = 2,
                   phi_death = 1,
                   k_death = 2,
                   exp_noise = 1e6)

  steps_per_day <- 1
  R0_change <- int_unique$change
  date_R0_change <- as.Date(int_unique$dates_change)
  date_contact_matrix_set_change <- NULL
  n_particles <- 2

  # proposal kernel covariances
  proposal_kernel <- matrix(0.5, ncol=length(pars_init), nrow = length(pars_init))
  diag(proposal_kernel) <- 1
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)

  out <- squire::pmcmc(data = data,
                       n_mcmc = 20,
                       log_likelihood = NULL,
                       log_prior = NULL,
                       n_particles = 2,
                       steps_per_day = 1,
                       output_proposals = FALSE,
                       n_chains = 1,
                       replicates = 20,
                       burnin = 0,
                       squire_model = nimue_deterministic_model(use_dde = FALSE),
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       country = country)

  # basic argument out checks
  expect_named(
    out,
    c("output","parameters","model","odin_parameters",
      "replicate_parameters","pmcmc_results","interventions")
  )

  # and check plotting works
  # just replicates
  expect_s3_class(
    plot(out, "deaths", replicates = TRUE, ci = FALSE, summarise = FALSE),
    "gg")

  # and ci
  expect_s3_class(
    plot(out, "deaths", x_var = "date", date_0 = max(out$pmcmc_results$inputs$data$date)),
    "gg"
  )

})




#------------------------------------------------
test_that("pmcmc with vaccine pars", {

  Sys.setenv("SQUIRE_PARALLEL_DEBUG" = "TRUE")
  data <- read.csv(squire:::squire_file("extdata/example.csv"),stringsAsFactors = FALSE)
  interventions <- read.csv(squire:::squire_file("extdata/example_intervention.csv"))
  int_unique <- squire:::interventions_unique(interventions)
  reporting_fraction = 1
  country <- "Algeria"
  pars_init <- list('start_date'     = as.Date("2020-02-07"),
                    'R0'             = 2.5,
                    'Meff'           = 2)
  pars_min <- list('start_date'      = as.Date("2020-02-01"),
                   'R0'              = 1e-10,
                   'Meff'            = 0.1)
  pars_max <- list('start_date'      = as.Date("2020-02-20"),
                   'R0'              = 5,
                   'Meff'            = 5)
  pars_discrete <- list('start_date' = TRUE,
                        'R0'         = FALSE,
                        'Meff'       = FALSE)
  pars_obs <- list(phi_cases = 0.1,
                   k_cases = 2,
                   phi_death = 1,
                   k_death = 2,
                   exp_noise = 1e6)

  steps_per_day <- 1
  R0_change <- rep(1, 4)
  date_R0_change <- as.Date(int_unique$dates_change)
  date_vaccine_change <- date_R0_change - 30
  max_vaccine <- c(100000, 200000, 300000, 400000)
  baseline_max_vaccine <- 0
  date_contact_matrix_set_change <- NULL
  n_particles <- 2

  # proposal kernel covariances
  proposal_kernel <- matrix(0.5, ncol=length(pars_init), nrow = length(pars_init))
  diag(proposal_kernel) <- 1
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)

  out <- squire::pmcmc(data = data,
                       n_mcmc = 100,
                       log_likelihood = NULL,
                       log_prior = NULL,
                       n_particles = 2,
                       steps_per_day = 1,
                       output_proposals = FALSE,
                       n_chains = 1,
                       replicates = 20,
                       burnin = 0,
                       squire_model = nimue_deterministic_model(use_dde = TRUE),
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       max_vaccine = max_vaccine*10,
                       baseline_max_vaccine = baseline_max_vaccine,
                       date_vaccine_change = date_vaccine_change,
                       baseline_vaccine_efficacy_infection = rep(0.8,17),
                       date_vaccine_efficacy_infection_change = date_R0_change,
                       vaccine_efficacy_infection = list(rep(0.6,17), rep(0.7,17), rep(0.8,17), rep(0.9,17)),
                       baseline_hosp_bed_capacity = 1000,
                       country = country)

  # basic argument out checks
  expect_named(
    out,
    c("output","parameters","model","odin_parameters",
      "replicate_parameters","pmcmc_results","interventions")
  )

  # and check plotting works
  # just replicates
  expect_s3_class(
    plot(out, "deaths", replicates = TRUE, ci = FALSE, summarise = FALSE),
    "gg")

  # and ci
  expect_s3_class(
    plot(out, "deaths", x_var = "date", date_0 = max(out$pmcmc_results$inputs$data$date)),
    "gg"
  )

  expect_s3_class(
    plot(out, particle_fit = TRUE),
    "gg"
  )

  expect_true("date_vaccine_change" %in% names(out$interventions))
  expect_true("max_vaccine" %in% names(out$interventions))
  expect_true("max_vaccine" %in% names(out$parameters))
  expect_equal(
    out$pmcmc_results$inputs$model_params$max_vaccine,
    c(baseline_max_vaccine, max_vaccine*10)
  )


})


#------------------------------------------------
test_that("pmcmc vaccine can find", {

  # curve with vaccine
  max_vaccine <- c(1000000)
  r <- run("Iran", R0 = 3, max_vaccine = c(0,max_vaccine), tt_vaccine = c(0, 30), vaccine_efficacy_infection = rep(0.99,17))
  data <- format(r, summaries = "deaths", date_0 = "2020-01-01")
  data <- data[data$compartment == "deaths",]
  data$value <- as.integer(data$value)
  data <- na.omit(data)
  data <- data[1:150,c("date", "value")]
  names(data)[2] <- "deaths"
  data_vacc <- data
  date_vaccine_change <- as.Date(c("2020-01-01"))+30

  # curve due to R changes just to show similar thing could be explained without vaccine
  r2 <- run("Iran", R0 = c(2.38,0.6), max_vaccine = 0, tt_R0 = c(0,91))
  data <- format(r2, summaries = "deaths", date_0 = "2019-12-25")
  data <- data[data$compartment == "deaths",]
  data$value <- as.integer(data$value)
  data <- na.omit(data)
  data <- data[1:150,c("date", "value")]
  names(data)[2] <- "deaths"
  data_no_vacc <- data
  data_no_vacc <- data_no_vacc[data_no_vacc$deaths > 0,]

  # check that they are similarish:
  # plot(r2, "deaths", date_0 = "2019-12-21", x_var="date") +
  # geom_point(aes(date, value), format(r, summaries = "deaths",date_0 = "2020-01-01") %>%
  # filter(compartment == "deaths"), inherit.aes = FALSE)

  # First let's see if it can find the correct vaccine solution if given R0_change initials

  # what meff is needed for a given R0_change to get 0.6.
  # Answer 2.43
  # squire:::evaluate_Rt_pmcmc(c(1,0.2),2.4,as.Date(c("2020-01-01", "2020-02-02")),pars = list("Meff" = 2.43),Rt_args = list())

  pars_init = list('start_date'     = as.Date("2019-12-25"),
                   'R0'             = 2.38,
                   'Meff'           = 2.42)
  pars_min = list('start_date'      = as.Date("2019-12-10"),
                  'R0'              = 1e-10,
                  'Meff'            = -2)
  pars_max = list('start_date'      = as.Date("2020-01-07"),
                  'R0'              = 5,
                  'Meff'            = 5)
  pars_discrete = list('start_date' = TRUE,
                       'R0'         = FALSE,
                       'Meff'       = FALSE)
  pars_obs = list(phi_cases = 0.1,
                  k_cases = 2,
                  phi_death = 1,
                  k_death = 2,
                  exp_noise = 1e6)
  steps_per_day = 1
  R0_change = c(0.2)
  date_R0_change <- as.Date(c("2019-12-25")) + c(91)

  # proposal kernel covriance
  proposal_kernel <- matrix(0.5, ncol=length(pars_init), nrow = length(pars_init))
  diag(proposal_kernel) <- 1
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)

  # now fit to this data that we know can be explained by the initails in
  # the absence of vaccines, but with vaccines present has to be explained by
  # a different set of pars
  n_mcmc <- 10
  out <- squire::pmcmc(data = data_no_vacc,
                       n_mcmc = n_mcmc,
                       log_likelihood = NULL,
                       log_prior = NULL,
                       n_particles = 2,
                       steps_per_day = 1,
                       output_proposals = FALSE,
                       n_chains = 3,
                       replicates = 20,
                       seeding_cases = 20,
                       burnin = 0,
                       squire_model = nimue_deterministic_model(use_dde = TRUE),
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       R0_change = 0.2,
                       date_R0_change = date_R0_change,
                       max_vaccine = max_vaccine,
                       baseline_max_vaccine = 0,
                       date_vaccine_change = date_vaccine_change,
                       baseline_vaccine_efficacy_infection = rep(0.99,17),
                       start_adaptation = as.integer(n_mcmc/5),
                       country = "Iran")

  n_mcmc <- 10
  out2 <- squire::pmcmc(data = data_no_vacc,
                       n_mcmc = n_mcmc,
                       log_likelihood = NULL,
                       log_prior = NULL,
                       n_particles = 2,
                       steps_per_day = 1,
                       output_proposals = FALSE,
                       n_chains = 3,
                       replicates = 20,
                       burnin = 0,
                       seeding_cases = 20,
                       squire_model = nimue_deterministic_model(use_dde = TRUE),
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       R0_change = 0.2,
                       date_R0_change = date_R0_change,
                       max_vaccine = 0,
                       baseline_max_vaccine = 0,
                       date_vaccine_change = date_vaccine_change,
                       baseline_vaccine_efficacy_infection = rep(0.99,17),
                       start_adaptation = as.integer(n_mcmc/5),
                       country = "Iran")

  expect_lt(max(out2$pmcmc_results$inputs$model_params$max_vaccine),
            max(out$pmcmc_results$inputs$model_params$max_vaccine))

})
