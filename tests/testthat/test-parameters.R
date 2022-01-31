context("Parameters")
#library(nimue)
library(squire)

test_that("default parameter list works", {
  d1 <- default_probs()
  expect_type(d1, "list")
  expect_named(d1, c("prob_hosp", "prob_severe", "prob_non_severe_death_treatment",
                     "prob_non_severe_death_no_treatment", "prob_severe_death_treatment",
                     "prob_severe_death_no_treatment",  "p_dist",
                     "rel_infectiousness", "rel_infectiousness_vaccinated",
                     'prob_hosp_multiplier', 'tt_prob_hosp_multiplier',
                     'prob_severe_multiplier', 'tt_prob_severe_multiplier'))
  v1 <- default_vaccine_pars()
  expect_type(v1, "list")
  expect_named(v1, c("dur_R", "tt_dur_R", "dur_V",
                     "vaccine_efficacy_infection", "tt_vaccine_efficacy_infection",
                     "vaccine_efficacy_disease", "tt_vaccine_efficacy_disease",
                     "max_vaccine", "tt_vaccine", "dur_vaccine_delay",
                     "vaccine_coverage_mat"))
})


test_that("beta_est_infectiousness input checks work", {
  mm <- matrix(runif(4), ncol = 2)
  mm_na <- mm
  mm_na[1] <- NA

  expect_error(beta_est_infectiousness(dur_IMild = 0, dur_ICase = 5,
                                       prob_hosp = c(0.2,0.1),
                                       rel_infectiousness = rep(1,2),
                                       mixing_matrix = mm, R0 = 2),
               "dur_IMild must be greater than zero")

  expect_error(beta_est_infectiousness(dur_IMild = 2, dur_ICase = 0,
                                       prob_hosp = c(0.2,0.1),
                                       rel_infectiousness = rep(1,2),
                                       mixing_matrix = mm, R0 = 2),
               "dur_ICase must be greater than zero")

  expect_error(beta_est_infectiousness(dur_IMild = c(2,2), dur_ICase = 5,
                                       prob_hosp = c(0.2,0.1),
                                       rel_infectiousness = rep(1,2),
                                       mixing_matrix = mm, R0 = 2),
               "dur_IMild must be of length 1")

  expect_error(beta_est_infectiousness(dur_IMild = 2, dur_ICase = 5,
                                       prob_hosp = c(0.2,NA),
                                       rel_infectiousness = rep(1,2),
                                       mixing_matrix = mm, R0 = 2),
               "prob_hosp must not contain NAs")

  expect_error(beta_est_infectiousness(dur_IMild = 2, dur_ICase = 5,
                                       prob_hosp = c(0.2,0.2),
                                       rel_infectiousness = rep(1,2),
                                       mixing_matrix = mm_na, R0 = 2),
               "mixing_matrix must not contain NAs")



})

test_that("seeding_age_order works", {

  pop <- data.frame(n = rep(1000000, 17))
  mm <- matrix(1, 16, 16)
  mm[,16] <- 2

  # Increasing age
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = sum(pop$n),
    ICU_bed_capacity = sum(pop$n),
    dur_R = Inf,
    max_vaccine = 0,
    seed = 1,
    replicates = 1,
    R0 = 2,
    seeding_cases = 20,
    seeding_age_order = 1:17,
    time_period = 10
  )

  expect_true(all(m1$model$contents()$initial_E1[1:5,1] == c(2,2,2,1,1)))

  # Decreasing age
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = sum(pop$n),
    ICU_bed_capacity = sum(pop$n),
    dur_R = Inf,
    max_vaccine = 0,
    seed = 1,
    replicates = 1,
    R0 = 2,
    seeding_cases = 20,
    seeding_age_order = rev(1:17),
    time_period = 10
  )

  expect_true(all(m2$model$contents()$initial_E1[13:17,1] == c(1,1,2,2,2)))

})

test_that("beta_est_infectiousness value return checks", {
  mm <- matrix(runif(4), ncol = 2)
  mm_na <- mm

  # same infectiousness
  beta_1 <- beta_est_infectiousness(dur_IMild = 3,
                                    dur_ICase = 5,
                                    prob_hosp = c(0.1, 0.5),
                                    rel_infectiousness = c(1, 1),
                                    mixing_matrix = mm,
                                    R0 = 2)

  # less in one age group
  beta_2 <- beta_est_infectiousness(dur_IMild = 3,
                                    dur_ICase = 5,
                                    prob_hosp = c(0.1, 0.5),
                                    rel_infectiousness = c(0.5, 1),
                                    mixing_matrix = mm,
                                    R0 = 2)
  expect_lt(beta_1, beta_2)

  # and test that values make sense in odin runs in homogenous population
  pop <- data.frame(n = rep(1000000, 17))
  mm <- matrix(1, 16, 16)
  mm[,16] <- 2 # this then creates a homogenous matrix in 17d

  set.seed(123L)
  R0 <- 1

  # Starting model
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = sum(pop$n),
    ICU_bed_capacity = sum(pop$n),
    dur_R = Inf,
    max_vaccine = 0,
    seed = 1,
    replicates = 1,
    R0 = R0,
    seeding_cases = 17,
    seeding_age_order = 1:17,
    time_period = 50,
  )

  # One with infectiousness in
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = sum(pop$n),
    ICU_bed_capacity = sum(pop$n),
    dur_R = Inf,
    rel_infectiousness = c(rep(0.5,2), rep(1,15)),
    max_vaccine = 0,
    seed = 1,
    replicates = 1,
    R0 = R0,
    seeding_cases = 17,
    seeding_age_order = 1:17,
    time_period = 50,
  )

  # double check this is the homogenous model in demography
  expect_true(unique(m1$parameters$population) == 1e6)

  # double check this is the homogenous model in mixing
  expect_true(unique(as.numeric(m1$model$contents()$m)) == 1e-6)

  # grab infections
  i1 <- format(m1, summaries = "infections")
  i1 <- i1[i1$compartment == "infections",]
  i2 <- format(m2, summaries = "infections")
  i2 <- i2[i2$compartment == "infections",]

  # check that they are comparable totals
  expect_lt(abs(sum(i1$value[-1])-sum(i2$value[-1])), 0.5)

  # check that they are staying approximately flat over 50 days
  expect_lt(diff(i1$value[c(10,20)]),0.01)
  expect_lt(diff(i1$value[c(10,30)]),0.01)
  expect_lt(diff(i1$value[c(10,40)]),0.01)

  # check that they are staying approximately flat over 50 days
  expect_lt(diff(i2$value[c(10,20)]),0.01)
  expect_lt(diff(i2$value[c(10,30)]),0.01)
  expect_lt(diff(i2$value[c(10,40)]),0.01)

})

# time varying vaccine efficacy parameters
test_that("correct lengths on time varying efficacies", {

  # check for correct length agreements
  expect_error(
    r <- nimue::run(country = "Iran", vaccine_efficacy_infection = list(rep(0.5,17),rep(0.9,17))) ,
    "vaccine_efficacy_infection must be of length 1"
  )

  expect_error(
    r <- nimue::run(country = "Iran", vaccine_efficacy_disease = list(rep(0.5,17),rep(0.9,17))) ,
    "vaccine_efficacy_disease must be of length 1"
  )

  # check to see that it is being correctly interpolated in outputs
  r <- run("Iran",
           R0 = 1.5,
           vaccine_efficacy_disease = list(rep(0,17),rep(1,17)),
           tt_vaccine_efficacy_disease = c(0, 200),
           vaccine_efficacy_infection = rep(0,17),
           max_vaccine = c(1000000),
           dur_V = Inf,
           vaccine_coverage_mat = nimue::strategy_matrix("All", 1)
  )

  # if correct we should have deaths falling sharp the day after they have come in
  output <- format(r, summaries = c("deaths", "infections"))
  deaths <- output[output$compartment == "deaths",]
  expect_gt(deaths$value[deaths$t == 201], deaths$value[deaths$t == 200])
  expect_lt(deaths$value[deaths$t == 202], deaths$value[deaths$t == 201])

  # but infections should stay increasing
  infections <- output[output$compartment == "infections",]
  expect_gt(infections$value[infections$t == 201], infections$value[infections$t == 200])
  expect_gt(infections$value[infections$t == 202], infections$value[infections$t == 201])

  # but if we flip the profiles round
  r2 <- run("Iran",
            R0 = 1.5,
            vaccine_efficacy_disease = list(rep(0,17),rep(0,17)),
            tt_vaccine_efficacy_disease = c(0, 200),
            vaccine_efficacy_infection = list(rep(0,17), rep(1,17)),
            tt_vaccine_efficacy_infection = c(0, 200),
            max_vaccine = c(1000000),
            dur_V = Inf,
            vaccine_coverage_mat = nimue::strategy_matrix("All", 1)
  )

  # with them flipped deaths will still keep increasing due to the delay till admission etc
  output2 <- format(r2, summaries = c("deaths", "infections"))
  deaths2 <- output2[output2$compartment == "deaths",]
  expect_gt(deaths2$value[deaths2$t == 201], deaths2$value[deaths2$t == 200])
  expect_gt(deaths2$value[deaths2$t == 202], deaths2$value[deaths2$t == 201])

  # but the way infections have been encoded these will fall immediately to 0
  infections2 <- output2[output2$compartment == "infections",]
  expect_gt(infections2$value[infections2$t == 200], infections2$value[infections2$t == 199])
  expect_lt(infections2$value[infections2$t == 201], infections2$value[infections2$t == 200])


})



test_that("rel_infectiousness_vaccine works", {

  r1 <- nimue::run("Iran", max_vaccine = 200000)
  infs1 <- nimue::format(r1, summaries = "infections")
  r2 <- nimue::run("Iran", max_vaccine = 200000, rel_infectiousness_vaccinated = 0.2)
  infs2 <- nimue::format(r2, summaries = "infections")

  expect_gt(
    max(infs1$value[infs1$compartment == "infections"], na.rm=TRUE),
    max(infs2$value[infs2$compartment == "infections"], na.rm=TRUE)
  )

})

