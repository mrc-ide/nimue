test_that("Pre-programmed strategies work", {
  expect_error(strategy_matrix("a"), "Strategy must be one of: All, Elderly, Working Elderly Children, Risk Elderly Working Children")
  expect_error(strategy_matrix("All", -1), "max_coverage must be between 0 and 1")
  expect_error(strategy_matrix("All", 2), "max_coverage must be between 0 and 1")
  expect_error(strategy_matrix("All", 0.5, -1), "risk_proportion must be between 0 and 1")
  expect_error(strategy_matrix("All", 0.5, 2), "risk_proportion must be between 0 and 1")

  expect_equal(ncol(strategy_matrix("All")), 17)
  expect_equal(nrow(strategy_matrix("All")), 1)

  expect_equal(ncol(strategy_matrix("Elderly")), 17)
  expect_equal(nrow(strategy_matrix("Elderly")), 17)

  expect_equal(ncol(strategy_matrix("Working Elderly Children")), 17)
  expect_equal(nrow(strategy_matrix("Working Elderly Children")), 3)

  expect_equal(ncol(strategy_matrix("Risk Elderly Working Children")), 17)
  expect_equal(nrow(strategy_matrix("Risk Elderly Working Children")), 4)
})
