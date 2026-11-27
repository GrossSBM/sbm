## Option classifInit of estimateSimpleSBM() / estimateBipartiteSBM(): additional initializations

## A food web with two compartments coupled by top predators (rows eat columns). The default
## exploration stops at 3 blocks (top predators / consumers / basal species): separating the two
## channels requires to split the consumers and their resources at the same time.
sim_compartments <- function(seed = 6) {
  set.seed(seed)
  alpha <- matrix(c(.05, .45, .45, .00, .00,
                    .00, .05, .00, .60, .00,
                    .00, .00, .05, .00, .60,
                    .00, .00, .00, .00, .00,
                    .00, .00, .00, .00, .00), 5, 5, byrow = TRUE)
  z <- rep(1:5, c(6, 10, 10, 12, 12))
  P <- alpha[z, z]
  Y <- matrix(rbinom(length(P), 1, P), nrow(P))
  diag(Y) <- 0
  perm <- sample(length(z))
  list(Y = Y[perm, perm], z = z[perm])
}
same_partition <- function(z1, z2) {
  length(unique(z1)) == length(unique(z2)) && nrow(unique(cbind(z1, z2))) == length(unique(z1))
}
opts <- list(verbosity = 0, plot = FALSE, nbCores = 1)

test_that("classifInit = 'kmeans' recovers the compartments of a food web", {
  sim <- sim_compartments()
  fit0 <- estimateSimpleSBM(sim$Y, "bernoulli", directed = TRUE, estimOptions = opts)
  set.seed(1)
  fit1 <- estimateSimpleSBM(sim$Y, "bernoulli", directed = TRUE,
                            estimOptions = c(opts, list(classifInit = "kmeans")))
  expect_equal(fit0$nbBlocks, 3)                 # the default exploration misses the compartments
  expect_equal(fit1$nbBlocks, 5)
  expect_true(same_partition(fit1$memberships, sim$z))
  ## adding initializations can only improve the best model found for each number of blocks
  K <- intersect(fit0$storedModels$nbBlocks, fit1$storedModels$nbBlocks)
  icl0 <- fit0$storedModels$ICL[match(K, fit0$storedModels$nbBlocks)]
  icl1 <- fit1$storedModels$ICL[match(K, fit1$storedModels$nbBlocks)]
  expect_true(all(icl1 >= icl0 - 1e-8))
})

test_that("classifInit accepts clusterings supplied by the user", {
  sim <- sim_compartments()
  fit <- estimateSimpleSBM(sim$Y, "bernoulli", directed = TRUE,
                           estimOptions = c(opts, list(classifInit = sim$z)))
  expect_equal(fit$nbBlocks, 5)
  expect_true(same_partition(fit$memberships, sim$z))
  ## a list of clusterings, with labels of any type
  fit2 <- estimateSimpleSBM(sim$Y, "bernoulli", directed = TRUE,
                            estimOptions = c(opts, list(classifInit = list(letters[sim$z], rep(1:2, 25)))))
  expect_equal(fit2$nbBlocks, 5)
  expect_error(estimateSimpleSBM(sim$Y, "bernoulli", directed = TRUE,
                                 estimOptions = c(opts, list(classifInit = 1:3))),
               "vector of length 50")
})

test_that("classifInit works for undirected networks and other models", {
  set.seed(2)
  mySampler <- sampleSimpleSBM(60, c(.5, .5), list(mean = matrix(c(.8, .1, .1, .6), 2)), directed = FALSE)
  set.seed(3)
  fit <- estimateSimpleSBM(mySampler$networkData, "bernoulli", directed = FALSE,
                           estimOptions = c(opts, list(classifInit = "kmeans")))
  expect_equal(fit$nbBlocks, 2)
  expect_true(same_partition(fit$memberships, mySampler$memberships))
  set.seed(4)
  mySampler <- sampleSimpleSBM(60, c(.5, .5), list(mean = matrix(c(6, 1, 1, 4), 2)), model = "poisson")
  set.seed(5)
  fit <- estimateSimpleSBM(mySampler$networkData, "poisson",
                           estimOptions = c(opts, list(classifInit = "kmeans")))
  expect_equal(fit$nbBlocks, 2)
})

test_that("classifInit works for bipartite SBM", {
  set.seed(7)
  mySampler <- sampleBipartiteSBM(c(40, 50), list(c(.5, .5), c(.3, .3, .4)),
                                  list(mean = matrix(c(.9, .5, .1, .1, .5, .9), 2, 3, byrow = TRUE)))
  Y <- mySampler$networkData
  fit0 <- estimateBipartiteSBM(Y, estimOptions = opts)
  set.seed(8)
  fit1 <- estimateBipartiteSBM(Y, estimOptions = c(opts, list(classifInit = "kmeans")))
  expect_equal(fit1$nbBlocks, fit0$nbBlocks)
  expect_true(max(fit1$storedModels$ICL) >= max(fit0$storedModels$ICL) - 1e-8)
  ## clustering supplied by the user: list (row, col)
  z <- mySampler$memberships
  fit2 <- estimateBipartiteSBM(Y, estimOptions = c(opts, list(classifInit = list(z[[1]], z[[2]]))))
  expect_equal(unname(fit2$nbBlocks), c(2, 3))
  expect_error(estimateBipartiteSBM(Y, estimOptions = c(opts, list(classifInit = list(z[[1]])))),
               "list of two vectors")
})
