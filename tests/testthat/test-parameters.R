context("Parameters")
library(nimue)
library(squire)

test_that("default parameter list works", {
  d1 <- default_probs()
  expect_type(d1, "list")
  expect_named(d1, c("prob_hosp", "prob_severe", "prob_non_severe_death_treatment",
                     "prob_non_severe_death_no_treatment", "prob_severe_death_treatment",
                     "prob_severe_death_no_treatment",  "p_dist"))
  v1 <- default_vaccine_pars()
  expect_type(v1, "list")
  expect_named(v1, c("dur_R", "vaccination_target", "dur_V", "vaccine_efficacy_infection",
                     "vaccine_efficacy_disease", "max_vaccine", "tt_vaccine", "dur_vaccine_delay"))
})
