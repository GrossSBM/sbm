test_that("Output type", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type,directed=TRUE,dimLabels=c("Actor"))
  expect_equal(netA$dimLabels[[1]],"Actor")
  expect_true(is_SBM(netA))
  })


test_that("Model versus data", {
  A <- matrix(rbinom(100,1,.2),10,10)
  type <- "simple"
  netA <- defineSBM(A,"bernoulli",type,directed=TRUE,dimLabels=c("Actor"))
  expect_error(defineSBM(A,"poisson",type,directed=TRUE,dimLabels=c("Actor")))
  expect_error(defineSBM(A*5,"bernoulli",type,directed=TRUE,dimLabels=c("Actor")))
  expect_error(defineSBM(A,"gaussian",type,directed=TRUE,dimLabels=c("Actor")))
  expect_error(defineSBM(A*0.4,"poisson",type,directed=TRUE,dimLabels=c("Actor")))
  expect_error(defineSBM(A*0.4,"bernoulli",type,directed=TRUE,dimLabels=c("Actor")))
})

