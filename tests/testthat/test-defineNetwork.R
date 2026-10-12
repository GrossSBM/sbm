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



test_that("Validation of the model argument", {
  A <- matrix(rbinom(100,1,.2),10,10)
  S <- A; S[lower.tri(S)] <- t(S)[lower.tri(S)]; diag(S) <- 0
  Z <- matrix(rnorm(100)*rbinom(100,1,.5),10,10)
  opts <- list(verbosity = 0, plot = FALSE, nbCores = 1)

  ## estimate functions: only the models with an implemented inference
  expect_error(estimateSimpleSBM(S, "bernouilli", estimOptions = opts), "must belong to: bernoulli, poisson, gaussian$")
  expect_error(estimateSimpleSBM(S, "ZIgaussian", estimOptions = opts), "must belong to")
  expect_error(estimateSimpleSBM(S, c("poisson", "bernoulli"), estimOptions = opts), "must belong to")
  expect_error(estimateBipartiteSBM(A, "bernouilli", estimOptions = opts), "must belong to")

  ## defineSBM: all the models, including ZIgaussian (used in the multipartite inference)
  expect_error(defineSBM(A, "bernouilli", "bipartite"), "must belong to: bernoulli, poisson, gaussian, ZIgaussian")
  netZ <- defineSBM(Z, "ZIgaussian", "bipartite", dimLabels = c("a", "b"))
  expect_true(is_SBM(netZ))
  expect_equal(netZ$modelName, "ZIgaussian")
})
