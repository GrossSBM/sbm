## Behaviors of sbm 0.4.7 that must be kept

test_that("the number of blocks of a sampler is known before the memberships are drawn", {
  mySampler <- SimpleSBM$new("bernoulli", 60, FALSE, c(.3, .3, .4), list(mean = diag(.4, 3) + .05))
  expect_equal(mySampler$nbBlocks, 3)
  expect_equal(mySampler$nbConnectParam, 6)
  mySampler$rMemberships(store = TRUE)
  expect_equal(mySampler$nbBlocks, 3)
})

test_that("a model which is not bernoulli can be fitted on a binary network, with a warning", {
  A <- sampleSimpleSBM(40, c(.5, .5), list(mean = diag(.5, 2) + .1))$networkData
  expect_warning(SimpleSBM_fit$new(A, "poisson", FALSE), "binary")
  expect_warning(SimpleSBM_fit$new(A, "gaussian", FALSE), "binary")
  expect_warning(BipartiteSBM_fit$new(A[, 1:30], "poisson"), "binary")
  ## models which are not adapted to the data are still rejected
  expect_error(SimpleSBM_fit$new(A * 2, "bernoulli", FALSE), "binary")
  expect_error(SimpleSBM_fit$new(A + 0.5, "poisson", FALSE), "integer")
})

test_that("BipartiteSBM_fit: storedModels and setModel() index the models by row number", {
  set.seed(4)
  net <- sampleBipartiteSBM(c(60, 50), list(c(.5, .5), c(.3, .3, .4)),
                            list(mean = matrix(c(.8, .1, .2, .6, .5, .05), 2)))
  fit <- estimateBipartiteSBM(net$networkData, estimOptions = list(verbosity = 0, plot = FALSE, nbCores = 1))
  models <- fit$storedModels
  expect_true(all(c("indexModel", "nbParams", "rowBlocks", "colBlocks", "nbBlocks", "ICL", "loglik") %in% names(models)))
  expect_equal(models$nbBlocks, models$rowBlocks + models$colBlocks)
  for (k in seq_len(nrow(models))) {
    fit$setModel(k)
    expect_equal(unname(fit$nbBlocks), c(models$rowBlocks[k], models$colBlocks[k]))
    expect_equal(fit$ICL, models$ICL[k])
  }
  ## as used by shinySbm: select a model by its total number of blocks
  n <- models$nbBlocks[2]
  fit$setModel(which(fit$storedModels$nbBlocks == n))
  expect_equal(sum(fit$nbBlocks), n)
})
