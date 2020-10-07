test_that("Tracking of proportion who have recieved vaccine works", {
  # Run with no vaccine
  r0 <- run("Angola", max_vaccine = 0, time_period = 50)
  # Run with small amount of vaccine
  r1 <- run("Angola", max_vaccine = 100, time_period = 50)
  # Run with large amount of vaccine
  r2 <- run("Angola", max_vaccine = 100000, time_period = 50,
            vaccine_coverage_mat = rep(0.05, 17))
  # Run with age targeting
  r3 <- run("Angola", max_vaccine = 100, time_period = 50,
            vaccine_coverage_mat = c(1, 1, rep(0, 13), 1, 1))

  index <- odin_index(r0$model)$prop_received

  expect_equal(as.vector(colSums(r0$output[,index,1])), rep(0, 17))
  expect_true(all(as.vector(colSums(r1$output[,index,1])) > rep(0, 17)))
  expect_true(all(as.vector(colSums(r2$output[,index,1])) > as.vector(colSums(r1$output[,index,1]))))
  expect_equal(round(as.vector(apply(r2$output[,index,1], 2, max)), 2), rep(0.05, 17))
  expect_equal(as.vector(colSums(r3$output[,index[3:15],1])), rep(0, 13))
  expect_true(all(as.vector(colSums(r3$output[,index[c(1,2,16,17)],1])) > rep(0, 4)))


 # prioritisation_matrix <- matrix(0, ncol = 17, nrow = 4)
 # prioritisation_matrix[1,5] <- 0.2
 #  prioritisation_matrix[2,17] <- 0.8
 #  prioritisation_matrix[3,5] <- 0.8
 #  prioritisation_matrix[4,] <- 0.8

})
