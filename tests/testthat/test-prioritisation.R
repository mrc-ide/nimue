

test_that("Tracking of proportion who have recieved vaccine works", {
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyr))
  # Run with no vaccine
  r0 <- run("Angola", max_vaccine = 0, time_period = 50)
  # Run with small amount of vaccine
  r1 <- run("Angola", max_vaccine = 100, time_period = 50)
  # Run with large amount of vaccine
  r2 <- run("Angola", max_vaccine = 10e6, time_period = 50)
  # Run with age targeting
  r3 <- run("Angola", max_vaccine = 10e6, time_period = 50,
            vaccine_coverage_mat = c(1, rep(0, 15), 1))
  # Run with "All" strategy
  m1 <- strategy_matrix("All")
  r4 <- run("Angola", max_vaccine = 10e6, time_period = 50,
            vaccine_coverage_mat = m1)
  # Run with "Elderly" strategy
  m2 <- strategy_matrix("Elderly")
  r5 <- run("Angola", max_vaccine = 10e6, time_period = 50,
            vaccine_coverage_mat = m2)
  # Run with "Working Elderly Children" strategy
  m3 <- strategy_matrix("Working Elderly Children")
  r6 <- run("Angola", max_vaccine = 1e6, time_period = 50,
            vaccine_coverage_mat = m3)
  # Run with "Risk Elderly Working Children" strategy
  m4 <- strategy_matrix("Risk Elderly Working Children", risk_proportion = 0.3)
  r7 <- run("Angola", max_vaccine = 1e6, time_period = 50,
            vaccine_coverage_mat = m4)

  # Extract the proportion vaccinated (by age)
  get_pv <- function(x){
    format(x, compartments = c(), summaries = c("N", "unvaccinated"), reduce_age = FALSE) %>%
      pivot_wider(id_cols = c(t, replicate, age_group), values_from = value, names_from = compartment) %>%
      mutate(prop_vaccinated = 1 - unvaccinated / N)
  }

  p0 <- get_pv(r0)
  expect_equal(sum(p0$prop_vaccinated), 0)

  p1 <- get_pv(r1)
  p2 <- get_pv(r2)
  expect_gt(sum(p1$prop_vaccinated), 0)
  expect_gt(sum(p2$prop_vaccinated), sum(p1$prop_vaccinated))
  p3 <- get_pv(r3)
  a1 <- c("0-5", "80+")
  expect_gt(sum(filter(p3, age_group == a1[1])$prop_vaccinated), 0)
  expect_gt(sum(filter(p3, age_group == a1[2])$prop_vaccinated), 0)
  expect_equal(sum(filter(p3, !age_group %in% a1)$prop_vaccinated), 0)
  p4 <- get_pv(r4)
  expect_true(all(filter(p4, t == 2)$prop_vaccinated > 0))
  expect_equal(var(filter(p4, t == 2)$prop_vaccinated > 0), 0)
  expect_equal(var(filter(p4, t == max(t))$prop_vaccinated > 0), 0)
  p5 <- get_pv(r5)
  expect_equal(sum(filter(p5, t == 2, age_group != "80+")$prop_vaccinated), 0)
  expect_gt(filter(p5, t == 2, age_group == "80+")$prop_vaccinated, 0)
  expect_true(all(filter(p5, t == max(t))$prop_vaccinated > 0))
  p6 <- get_pv(r6)
  d1 <- which(diff(filter(p6, age_group == "35-40")$prop_vaccinated)>0)[1]
  d2 <- which(diff(filter(p6, age_group == "80+")$prop_vaccinated)>0)[1]
  d3 <- which(diff(filter(p6, age_group == "0-5")$prop_vaccinated)>0)[1]
  expect_gt(d2, d1)
  expect_gt(d3, d2)
  p7 <- get_pv(r7)
  d4 <- which(diff(filter(p7, age_group == "35-40")$prop_vaccinated)>0)[1]
  d5 <- which(diff(filter(p7, age_group == "80+")$prop_vaccinated)>0)[1]
  d6 <- which(diff(filter(p7, age_group == "0-5")$prop_vaccinated)>0)[1]
  expect_gt(d5, d4)
  expect_gt(d6, d5)
})
