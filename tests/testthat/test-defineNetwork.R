test_that("Output type", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineNetwork(A,"poisson",type,directed=TRUE,dimLabels=list("Actor","Actor"))
  expect_equal(netA$dimLabels[[1]],"Actor")
})
