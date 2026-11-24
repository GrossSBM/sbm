## Interface between sbm and GREMLINS: order of the functional groups,
## initial clustering, missing values and parallel computing

simuThreeFG <- function(propNA = 0) {
  set.seed(2)
  # A has 1 block, B 2 blocks, C 3 blocks.
  # The two packages order the FGs differently: C, A, B in sbm (dimLabels),
  # C, B, A in GREMLINS (row FGs of all the networks first, then column FGs)
  zA <- rep(1, 30); zB <- rep(1:2, each = 20); zC <- rep(1:3, each = 20)
  mk <- function(z1, z2, p) {
    m <- matrix(rbinom(length(z1) * length(z2), 1, p[cbind(rep(z1, length(z2)), rep(z2, each = length(z1)))]), length(z1))
    if (propNA > 0) m[sample(length(m), propNA * length(m))] <- NA
    m
  }
  netCA <- defineSBM(mk(zC, zA, matrix(c(.8, .2, .5), 3, 1)), "bernoulli", "bipartite", dimLabels = c("C", "A"))
  netBC <- defineSBM(mk(zB, zC, matrix(c(.9, .1, .5, .1, .9, .1), 2, 3, byrow = TRUE)), "bernoulli", "bipartite", dimLabels = c("B", "C"))
  list(listSBM = list(netCA, netBC), z = list(A = zA, B = zB, C = zC))
}

test_that("multipartite fit: numbers of blocks are correctly labelled", {
  s <- simuThreeFG()
  fit <- estimateMultipartiteSBM(s$listSBM, list(verbosity = 0, nbCores = 1, nbBlocksRange = list(C = c(1, 5), A = c(1, 5), B = c(1, 5))))
  expect_equal(unname(fit$nbBlocks[c("A", "B", "C")]), c(1, 2, 3))
  best <- which.max(fit$storedModels$ICL)
  expect_equal(unlist(fit$storedModels[best, paste("nbBlocks", c("A", "B", "C"))], use.names = FALSE), c(1, 2, 3))
  # nbParams = connectivity parameters (3x1 + 2x3) + block proportions (0 + 1 + 2)
  expect_equal(fit$storedModels$nbParams[best], 3 + 6 + 3)
})

test_that("multipartite fit: classifInit (named, any order) is used", {
  s <- simuThreeFG()
  init <- s$z[c("B", "A", "C")]
  fit <- estimateMultipartiteSBM(s$listSBM, list(verbosity = 0, nbCores = 1, classifInit = init))
  expect_equal(unname(fit$nbBlocks[c("A", "B", "C")]), c(1, 2, 3))
  # fixed numbers of blocks
  fitFixed <- estimateMultipartiteSBM(s$listSBM, list(verbosity = 0, nbCores = 1, classifInit = init,
                                                      nbBlocksRange = list(C = c(3, 3), A = c(1, 1), B = c(2, 2))))
  expect_equal(unname(fitFixed$nbBlocks[c("A", "B", "C")]), c(1, 2, 3))
  expect_equal(aricode::ARI(fitFixed$memberships$C, s$z$C), 1)
})

test_that("multipartite fit with missing values", {
  s <- simuThreeFG(propNA = 0.2)
  fit <- estimateMultipartiteSBM(s$listSBM, list(verbosity = 0, nbCores = 1))
  expect_equal(unname(fit$nbBlocks[c("A", "B", "C")]), c(1, 2, 3))
  pred <- fit$predict()
  expect_false(anyNA(unlist(pred)))
})

test_that("multipartite fit in parallel restores the future plan of the user", {
  skip_on_cran()
  s <- simuThreeFG()
  future::plan(future::sequential)
  fit <- estimateMultipartiteSBM(s$listSBM, list(verbosity = 0, nbCores = 2, initBM = FALSE))
  expect_true(inherits(future::plan(), "sequential"))
  expect_equal(unname(fit$nbBlocks[c("A", "B", "C")]), c(1, 2, 3))
})

test_that("plots of a multipartite fit with a simple network and single blocks", {
  set.seed(1)
  netAA <- defineSBM(matrix(rbinom(900, 1, 0.3), 30, 30), "bernoulli", "simple", directed = TRUE, dimLabels = "A")
  netAB <- defineSBM(matrix(rpois(600, 2), 30, 20), "poisson", "bipartite", dimLabels = c("A", "B"))
  fit <- estimateMultipartiteSBM(list(netAA, netAB), list(verbosity = 0, nbCores = 1, initBM = FALSE))
  expect_equal(unname(fit$nbBlocks), c(1, 1))
  for (type in c("data", "expected", "meso")) expect_s3_class(plot(fit, type = type), "ggplot")
})
