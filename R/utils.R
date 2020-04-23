extractParamBM <- function(BMobject, Q) {

  model <- BMobject$model_name
  membership_name <-  BMobject$membership_name

  res <- list()

  if (model %in% c('bernoulli','bernoulli_multiplex','bernoulli_covariates')) {
    res$connectParam <- BMobject$model_parameters[Q][[1]]$pi
  }

  if ( model %in% c('poisson','poisson_covariates')) {
    res$connectParam  <- BMobject$model_parameters[Q][[1]]$lambda
  }
  if (model %in% c('poisson_covariates','bernoulli_covariates')) {

    res$covariatesParam <-  BMobject$model_parameters[Q][[1]]$beta
  }

  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    res$probMemberships <-  BMobject$memberships[[Q]]$Z
    res$memberships <- apply(res$probMemberships, 1, which.max)
    n <- nrow(res$probMemberships)
    res$blockProp <-  colSums(res$probMemberships)/n
    res$nbBlocks <- length(res$blockProp)
  }



  if (membership_name == 'LBM'){
    res$probMemberships[[1]] <-  BMobject$memberships[[Q]]$Z1
    res$probMemberships[[2]] <-  BMobject$memberships[[Q]]$Z2

    res$memberships[[1]] <- apply(res$probMemberships[[1]], 1, which.max)
    res$memberships[[2]] <- apply(res$probMemberships[[2]], 1, which.max)
    nRow <- nrow(res$probMemberships[[1]])
    nCol <- nrow(res$probMemberships[[1]])
    res$blockProp <-   list()
    res$blockProp[[1]] <- colSums(res$probMemberships[[1]])/nRow
    res$blockProp[[2]] <-  colSums(res$probMemberships[[2]])/nCol
    res$nbBlocks <- c(length(res$blockProp[[1]]),length(res$blockProp[[2]]))
    names(res$nbBlocks) <- c('row','col')
  }

  ########## ordering
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    o <- switch(model,
                poisson   = order(res$connectParam %*% matrix(res$blockProp,ncol = 1),decreasing = TRUE),
                bernoulli = order(res$connectParam %*% matrix(res$blockProp,ncol = 1),decreasing = TRUE),
                1:res$nbBlocks
    )
    res$blockProp <- res$blockProp[o]
    res$connectParam <- res$connectParam[o,o]
    res$probMemberships <- res$probMemberships[, o]
    res$memberships <- apply(res$probMemberships, 1, which.max)
  }

  if (membership_name == 'LBM'){
    oRow <- switch(model,
                   poisson   = order(res$connectParam %*% matrix(res$blockProp[[2]], ncol = 1),decreasing = TRUE),
                   bernoulli = order(res$connectParam %*% matrix(res$blockProp[[2]], ncol = 1),decreasing = TRUE),
                   1:res$nbBlocks[1]
    )
    oCol <- switch(model,
                   poisson   = order(c(matrix(res$blockProp[[1]],nrow = 1) %*% res$connectParam),decreasing = TRUE),
                   bernoulli = order(c(matrix(res$blockProp[[1]],nrow = 1) %*% res$connectParam),decreasing = TRUE),
                   1:res$nbBlocks[2]
    )

    res$blockProp[[1]] <- res$blockProp[[1]][oRow]
    res$blockProp[[2]] <- res$blockProp[[2]][oCol]
    res$connectParam <- res$connectParam[oRow,oCol]
    res$probMemberships[[1]] <- res$probMemberships[[1]][ , oRow, drop = FALSE]
    res$probMemberships[[2]] <- res$probMemberships[[2]][ , oCol, drop = FALSE]
    res$memberships[[1]]  <- apply(res$probMemberships[[1]] , 1, which.max)
    res$memberships[[2]]  <- apply(res$probMemberships[[2]] , 1, which.max)

  }

  res
}


as_indicator <- function(clustering) {
  K <- length(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

as_clustering <- function(indicator) {
  cl <- apply(indicator, 1, which.max)
  cl
}

## Some utils function for math
.xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))
.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}
