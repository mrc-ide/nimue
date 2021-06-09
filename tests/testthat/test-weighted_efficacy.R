test_that("assign doses", {
  expect_equal(assign_doses(0, c(10, 10, 10)), c(0, 0, 0))
  expect_equal(assign_doses(30, c(10, 10, 10)), c(10, 10, 10))
  expect_equal(assign_doses(15, c(10, 10, 10)), c(5, 5, 5))
  expect_equal(assign_doses(16, c(10, 10, 10)), c(5, 5, 6))
})

test_that("coverage", {
  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, 2, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  expect_equal(coverage(dose_times, 1), c(2/3, 0, 1))
  expect_equal(coverage(dose_times, 2), c(2/3, 1, 0))
})

test_that("eligible for second", {
  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, NA, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  expect_equal(eligable_for_second(dose_times, 1, 14),
               list(c(FALSE, FALSE, FALSE),
                    c(FALSE, FALSE, FALSE),
                    c(FALSE, FALSE, FALSE)))
  expect_equal(eligable_for_second(dose_times, 1, 0),
               list(c(TRUE, FALSE, FALSE),
                    c(FALSE, FALSE, FALSE),
                    c(TRUE, FALSE, FALSE)))
  expect_equal(eligable_for_second(dose_times, 200, 14),
               list(c(TRUE, FALSE, FALSE),
                    c(FALSE, FALSE, FALSE),
                    c(TRUE, TRUE, TRUE)))
})

test_that("target_pop", {
  dose_times <- list(matrix(c(1, 2, NA, NA, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, NA, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  # Dose 1
  expect_equal(target_pop(dose_number = 1, dose_times, prioritisation = rep(1, 3),
                          t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3)), c(1, 3, 0))
  # Dose 1 as a function of prioritisation matrix
  expect_equal(target_pop(dose_number = 1, dose_times, prioritisation = c(0, 1, 0),
                          t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3)), c(0, 3, 0))
  # Dose 2 - none as all d2_prioritise set to FALSE
  expect_equal(target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                          t = 1, dose_period = 14, d2_prioritise = rep(FALSE, 3)), c(0, 0, 0))
  # Dose 2 - none as too soon after dose 1
  expect_equal(target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                          t = 1, dose_period = 14, d2_prioritise = rep(TRUE, 3)), c(0, 0, 0))
  # Dose 2
  expect_equal(target_pop(dose_number = 2, dose_times, prioritisation = c(1, 1, 1),
                          t = 15, dose_period = 14, d2_prioritise = rep(TRUE, 3)), c(1, 0, 1))
})

test_that("administer first dose", {
  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, 2, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  expect_equal(administer_first_dose(c(1, 3, 0), dose_times, t = 10),
               list(matrix(c(1, 2, 10, 2, 3, NA), nrow = 3),
                    matrix(c(10, 10, 10, 2, 3, 4), nrow = 3),
                    matrix(c(1, 2, 2, NA, NA, NA), nrow = 3)))
  expect_equal(administer_first_dose(c(0, 0, 0), dose_times, t = 10),
               dose_times)
})

test_that("administer second dose", {
  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(NA, NA, NA, 2, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  expect_equal(administer_second_dose(c(0, 0, 3), dose_times, t = 100, dose_period = 14),
               list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                    matrix(c(NA, NA, NA, 2, 3, 4), nrow = 3),
                    matrix(c(1, 2, 2, 100, 100, 100), nrow = 3)))
  expect_equal(administer_second_dose(c(0, 0, 0), dose_times, t = 100, dose_period = 14),
               dose_times)
  expect_equal(administer_second_dose(c(0, 0, 3), dose_times, t = 4, dose_period = 14),
               dose_times)
})

test_that("extract dose number", {
  dose_times <- list(matrix(c(1, 2, NA, 2, 3, NA), nrow = 3),
                     matrix(c(1, 1, 1, 2, 3, 4), nrow = 3),
                     matrix(c(1, 2, 2, NA, NA, NA), nrow = 3))
  disease_efficacy <-  c(0.5, 0.5)
  infection_efficacy <- c(0.01, 1)
  dn <- extract_dose_number(dose_times, 3)
  expect_equal(dn$t, rep(1:3, each = 3))
  expect_equal(dn$age_group, factor(rep(1:3, 3), levels = 1:3))
  expect_equal(dn$n_1dose, c(1, 3, 1, 1, 0, 2, 0, 0, 0))
  expect_equal(dn$n_2dose, c(0, 0, 0, 1, 1, 0, 1, 1, 0))

  we <- add_weighted_efficacy(dn, infection_efficacy, disease_efficacy, 0, 0)
  expect_equal(we$t, rep(1:3, each = 3))
  expect_equal(we$age_group, factor(rep(1:3, 3), levels = 1:3))
  for(age in 1:3){
    dn_temp <- dplyr::filter(dn, age_group == age)
    we_temp <- dplyr::filter(we, age_group == age)
    # Number vaccinated
    vx2 = cumsum(dn_temp$n_2dose) # Vaccinated with 2 doses
    vx1 = cumsum(dn_temp$n_1dose) - vx2# Vaccinated with 1 dose only (received 1 but not yet 2)
    # Number vaccine protected - lag between administration of dose 1 and protection
    protected2 = vx2
    protected1 = pmax(0, dplyr::lag(vx1, 0, default = 0) - vx2)
    protected <- cbind(protected1, protected2)
    expect_equal(we_temp$weighted_infection_efficacy,
                 apply(protected, 1, function(x){
                   weighted.mean(infection_efficacy, x)
                 })
    )
    expect_equal(we_temp$weighted_disease_efficacy,
                 apply(protected, 1, function(x){
                   weighted.mean(disease_efficacy, x)
                 })
    )
  }

  # With a lag from first dose to protection
  we_lag <- add_weighted_efficacy(dn, infection_efficacy, disease_efficacy, 1, 0)
  expect_equal(we$t, rep(1:3, each = 3))
  expect_equal(we$age_group, factor(rep(1:3, 3), levels = 1:3))
  for(age in 1:3){
    dn_temp <- dplyr::filter(dn, age_group == age)
    we_temp <- dplyr::filter(we_lag, age_group == age)
    # Number vaccinated
    vx2 = cumsum(dn_temp$n_2dose) # Vaccinated with 2 doses
    vx1 = cumsum(dn_temp$n_1dose) - vx2# Vaccinated with 1 dose only (received 1 but not yet 2)
    # Number vaccine protected - lag between administration of dose 1 and protection
    protected2 = vx2
    protected1 = pmax(0, dplyr::lag(vx1, 1, default = 0) - vx2)
    protected <- cbind(protected1, protected2)
    expect_equal(we_temp$weighted_infection_efficacy,
                 apply(protected, 1, function(x){
                   ifelse(sum(x) == 0, infection_efficacy[1], weighted.mean(infection_efficacy, x))
                 })
    )
    expect_equal(we_temp$weighted_disease_efficacy,
                 apply(protected, 1, function(x){
                   ifelse(sum(x) == 0, disease_efficacy[1], weighted.mean(disease_efficacy, x))
                 })
    )
  }

  # With a lag from second dose to protection
  we_lag <- add_weighted_efficacy(dn, infection_efficacy, disease_efficacy, 0, 1)
  expect_equal(we$t, rep(1:3, each = 3))
  expect_equal(we$age_group, factor(rep(1:3, 3), levels = 1:3))
  for(age in 1:3){
    dn_temp <- dplyr::filter(dn, age_group == age)
    we_temp <- dplyr::filter(we_lag, age_group == age)
    # Number vaccinated
    vx2 = cumsum(dn_temp$n_2dose) # Vaccinated with 2 doses
    vx1 = cumsum(dn_temp$n_1dose) - vx2# Vaccinated with 1 dose only (received 1 but not yet 2)
    # Number vaccine protected - lag between administration of dose 1 and protection
    protected2 = dplyr::lag(vx2, 1, default = 0)
    protected1 = pmax(0, dplyr::lag(vx1, 0, default = 0) - protected2)
    protected <- cbind(protected1, protected2)
    expect_equal(we_temp$weighted_infection_efficacy,
                 apply(protected, 1, function(x){
                   ifelse(sum(x) == 0, infection_efficacy[1], weighted.mean(infection_efficacy, x))
                 })
    )
    expect_equal(we_temp$weighted_disease_efficacy,
                 apply(protected, 1, function(x){
                   ifelse(sum(x) == 0, disease_efficacy[1], weighted.mean(disease_efficacy, x))
                 })
    )
  }
})

test_that("weighted efficacy", {
  # Without 2nd dose priority
  t1 <- weighted_efficacy(iso3c = "GHA",
                          N = 1000,
                          maxt = 365,
                          doses_per_day = rep(15, 365),
                          dose_period = 12 * 7,
                          delay_dose1 = 28,
                          delay_dose2 = 7,
                          prioritisation_matrix = nimue::strategy_matrix("Elderly"),
                          d2_prioritise = rep(FALSE, 17),
                          infection_efficacy = c(0.1, 0.9),
                          disease_efficacy = c(0.2, 0.8))
  # With 2nd dose priority
  t2 <- weighted_efficacy(iso3c = "GHA",
                          N = 1000,
                          maxt = 365,
                          doses_per_day = rep(15, 365),
                          dose_period = 12 * 7,
                          delay_dose1 = 28,
                          delay_dose2 = 7,
                          prioritisation_matrix = nimue::strategy_matrix("Elderly"),
                          d2_prioritise = rep(TRUE, 17),
                          infection_efficacy = c(0.1, 0.9),
                          disease_efficacy = c(0.2, 0.8))
  expect_type(t1, "list")
  expect_named(t1, c("t", "age_group", "n_1dose", "n_2dose", "weighted_infection_efficacy", "weighted_disease_efficacy"))
  expect_type(t2, "list")
  expect_named(t2, c("t", "age_group", "n_1dose", "n_2dose", "weighted_infection_efficacy", "weighted_disease_efficacy"))
})
