context("Output")
library(nimue)
library(squire)

test_that("output works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model, no vaccine
  m <- run(
    time_period = 2,
    population = pop$n,
    contact_matrix_set = mm,
    seed = 1,
    dt = 1,
    replicates = 1,
    seeding_cases = 20
  )
  o1 <- format(m)
  expect_identical(colnames(o1), c("value", "compartment", "t", "replicate"))
  o2 <- format(m, reduce_age = FALSE)
  expect_identical(colnames(o2), c("value", "compartment", "t", "replicate", "age_group"))
  o3 <- format(m, date_0 = "2020-01-01")
  expect_identical(colnames(o3), c("value", "compartment", "t", "replicate", "date"))
  o4 <- format(m, compartments = "E", summaries = NULL)
  expect_identical(levels(o4$compartment), "E")
  o5 <- format(m, compartments = NULL, summaries = "infections")
  expect_identical(levels(o5$compartment), "infections")

  expect_error(format(m, compartments = "Not there"), "Some compartments specified not output by model")
})
