test_that("Output type", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineSBM(A,"poisson",type,directed=TRUE,dimLabels=c("Actor"))
  expect_equal(netA$dimLabels[[1]],"Actor")
  expect_true(is_SBM(netA))
})
