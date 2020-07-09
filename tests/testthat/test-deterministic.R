context("Deterministic model")
library(nimue)
library(squire)

test_that("compare deterministic vaccine model to SEEIR model", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Reference model
  m1 <- run_deterministic_SEIR_model(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    seed = 1,
    dt = 0.5,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
  )
  oi1 <- odin_index(m1$model)

  # Vaccine model, no vaccine
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = 0,
    seed = 1,
    dt = 0.5,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
  )
  oi2 <- odin_index(m2$model)

  # Compare shared compartments
  compare_compartments <- names(oi1)[names(oi1) %in% names(oi2)][-1]

  for(i in seq_along(compare_compartments)){
    # Isolare squire output
    t1 <- m1$output[,unlist(oi1[compare_compartments[i]]),]
    # Isolate nimue output and collapse vaccine compartments if needed
    if(is.matrix(oi2[[compare_compartments[i]]])){
      t2 <- apply(oi2[[compare_compartments[i]]], 1, function(x, y){
        rowSums(y[1:199,x,1])
      }, y = m2$output)
    } else {
      t2 <- m2$output[1:199,unlist(oi2[compare_compartments[i]]),]
    }
    # Clear attributes
    attributes(t1) <- NULL
    attributes(t2) <- NULL

    expect_equal(t1, t2, tol = 0.00001)
  }

  # Check all vaccine-related compartments are 0
  expect_equal(sum(m2$output[,unlist(oi2[c("vaccinated", "priorvaccinated", "vaccines")]),]), 0)

  # Check population size is constant at specified level
  expect_equal(format(m2, "N", NULL)$value,
               rep(sum(pop$n), 200))

})


test_that("Vaccine on works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model 100% efficacy against infection
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = 10000,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    time_period = 10,
    seeding_cases = 20
  )

  # Check individuals reaching V
  expect_gt(sum(format(m1, "vaccinated", NULL)$value), 0)
  expect_gt(sum(format(m1, "priorvaccinated", NULL)$value), 0)

  # Check population size is constant at specified level
  expect_equal(format(m1, "N", NULL)$value,
               rep(sum(pop$n), 100))


  # Vaccine model 100% efficacy against infection (infinite vaccine duration)
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = 10000,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    time_period = 10,
    seeding_cases = 20,
    dur_V = Inf
  )

  # Check individuals reaching V
  expect_gt(sum(format(m2, "vaccinated", NULL)$value), 0)
  # But not leaving v
  expect_equal(sum(format(m2, "priorvaccinated", NULL)$value), 0)

  # Check population size is constant at specified level
  expect_equal(format(m2, "N", NULL)$value,
               rep(sum(pop$n), 100))

  # Vaccine model 100% efficacy against infection (infinite vaccine delay)
  m3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = 10000,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    time_period = 10,
    seeding_cases = 20,
    dur_vaccine_delay  = Inf
  )
  S_index <- odin_index(m3$model)$S

  # Unvaccainted S
  expect_gt(sum(m3$output[,S_index[,1],1]), 0)
  # Vaccinated not protected S
  expect_gt(sum(m3$output[,S_index[,2:3],1]), 0)
  # Vaccinated and protected S
  expect_equal(sum(m3$output[,S_index[,4:6],1]), 0)
})


test_that("Age targeting works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model 100% efficacy against infection yougest age group
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = 10000,
    vaccination_target = c(1, rep(0, 16)),
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )

  # Check individuals in youngest age group reaching V
  age_v <- format(m1, NULL, "vaccinated",  reduce_age = FALSE)
  expect_gt(sum(dplyr::filter(age_v, age_group == "0-5")$value), 0)
  expect_equal(sum(dplyr::filter(age_v, age_group != "0-5")$value), 0)
})

test_that("Time-varying works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model time varying
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    max_vaccine = c(0, 1000, 0),
    tt_vaccine = c(0, 10, 20),
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )

  # Check individuals in youngest age group reaching V
  t_v <- format(m1, "vaccines", NULL)
  expect_equal(sum(dplyr::filter(t_v, t < 10)$value), 0)
  expect_gt(sum(dplyr::filter(t_v, t >= 10, t <20)$value), 0)
  expect_equal(sum(dplyr::filter(t_v, t >= 21)$value), 0)
})

test_that("Efficacy against infection works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # No vaccine
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )
  # Vaccine 50% efficacy against infection
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    max_vaccine = 1000000,
    vaccine_efficacy_infection = rep(0.5, 17),
    vaccine_efficacy_disease = rep(0, 17)
  )
  # Vaccine 100% efficacy against infection
  m3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    max_vaccine = 1000000,
    vaccine_efficacy_infection = rep(1, 17),
    vaccine_efficacy_disease = rep(0, 17)
  )

  i1 <- sum(format(m1, NULL, "infections")$value)
  i2 <- sum(format(m2, NULL, "infections")$value)
  i3 <- sum(format(m3, NULL, "infections")$value)

  expect_gt(i1, i2)
  expect_gt(i2, i3)
})

test_that("Efficacy against disease works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # No vaccine
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )
  # Vaccine 50% efficacy against disease
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    max_vaccine = 1000000,
    vaccine_efficacy_disease = rep(0.5, 17),
    vaccine_efficacy_infection = rep(0, 17)
  )
  # Vaccine 100% efficacy against disease
  m3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    dt = 0.1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    max_vaccine = 1000000,
    vaccine_efficacy_disease = rep(1, 17),
    vaccine_efficacy_infection = rep(0, 17)
  )

  i1 <- sum(format(m1, NULL, "deaths")$value)
  i2 <- sum(format(m2, NULL, "deaths")$value)
  i3 <- sum(format(m3, NULL, "deaths")$value)

  expect_gt(i1, i2)
  expect_gt(i2, i3)
})


