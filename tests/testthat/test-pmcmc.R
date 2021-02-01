
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

  # proposal kernel covriance
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
                       squire_model = nimue_deterministic_model(),
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
