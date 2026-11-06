set.seed(20261106)

## E[Y_ij] = sum_{q,l} Z1[i,q] Z2[j,l] logistic(logit(mean[q,l]) + covarEffect[i,j]), entry by entry
brute_force_prediction <- function(Z1, Z2, mean, covarEffect) {
  P <- matrix(0, nrow(covarEffect), ncol(covarEffect))
  for (i in seq_len(nrow(P))) for (j in seq_len(ncol(P)))
    P[i, j] <- sum(outer(Z1[i, ], Z2[j, ]) * .logistic(.logit(mean) + covarEffect[i, j]))
  P
}

soft_memberships <- function(n, Q) {
  Z <- matrix(rgamma(n * Q, shape = 0.5), n, Q)
  Z / rowSums(Z)
}

test_that("SimpleSBM: bernoulli prediction with covariates is the expectation under soft memberships", {
  nbNodes <- 20
  connectParam <- list(mean = matrix(c(.7, .2, .1, .2, .5, .3, .1, .3, .6), 3, 3))
  covarList <- list(matrix(rnorm(nbNodes^2), nbNodes, nbNodes))
  mySBM <- SimpleSBM$new("bernoulli", nbNodes, FALSE, c(.3, .3, .4), connectParam,
                         covarParam = 1.5, covarList = covarList)
  Z <- soft_memberships(nbNodes, 3)
  mySBM$probMemberships <- Z

  expected <- brute_force_prediction(Z, Z, connectParam$mean, 1.5 * covarList[[1]])
  expect_equal(unname(mySBM$predict()), expected, tolerance = 1e-12)

  ## the plug-in formula logistic(Z logit(mean) Z' + covarEffect) is not the expectation
  plugin <- .logistic(Z %*% .logit(connectParam$mean) %*% t(Z) + 1.5 * covarList[[1]])
  expect_gt(max(abs(plugin - expected)), 1e-2)

  ## hard memberships: both formulas coincide
  mySBM$probMemberships <- as_indicator(as_clustering(Z))
  Zh <- mySBM$probMemberships
  expect_equal(unname(mySBM$predict()),
               .logistic(Zh %*% .logit(connectParam$mean) %*% t(Zh) + 1.5 * covarList[[1]]),
               tolerance = 1e-12)
})

test_that("SimpleSBM: predict() uses the covariates given in argument", {
  nbNodes <- 15
  connectParam <- list(mean = matrix(c(.6, .2, .2, .5), 2, 2))
  covarList <- list(matrix(rnorm(nbNodes^2), nbNodes, nbNodes))
  newCovarList <- list(matrix(rnorm(nbNodes^2), nbNodes, nbNodes))
  mySBM <- SimpleSBM$new("bernoulli", nbNodes, FALSE, c(.5, .5), connectParam,
                         covarParam = -0.8, covarList = covarList)
  Z <- soft_memberships(nbNodes, 2)
  mySBM$probMemberships <- Z

  expect_equal(unname(mySBM$predict(newCovarList)),
               brute_force_prediction(Z, Z, connectParam$mean, -0.8 * newCovarList[[1]]),
               tolerance = 1e-12)
})

test_that("BipartiteSBM: bernoulli prediction with covariates is the expectation under soft memberships", {
  nbNodes <- c(18, 12)
  connectParam <- list(mean = matrix(c(.7, .2, .1, .3, .6, .4), 2, 3))
  covarList <- list(matrix(rnorm(prod(nbNodes)), nbNodes[1], nbNodes[2]))
  mySBM <- BipartiteSBM$new("bernoulli", nbNodes, list(c(.5, .5), c(.3, .3, .4)), connectParam,
                            covarParam = 1.2, covarList = covarList)
  Z1 <- soft_memberships(nbNodes[1], 2)
  Z2 <- soft_memberships(nbNodes[2], 3)
  mySBM$probMemberships <- list(Z1, Z2)

  expected <- brute_force_prediction(Z1, Z2, connectParam$mean, 1.2 * covarList[[1]])
  expect_equal(unname(mySBM$predict()), expected, tolerance = 1e-12)

  plugin <- .logistic(Z1 %*% .logit(connectParam$mean) %*% t(Z2) + 1.2 * covarList[[1]])
  expect_gt(max(abs(plugin - expected)), 1e-2)
})
