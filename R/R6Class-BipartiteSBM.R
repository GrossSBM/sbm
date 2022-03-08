#' R6 class for Bipartite SBM
#'
#' @export
BipartiteSBM <-
  R6::R6Class(
    classname = "BipartiteSBM",
    inherit = SBM,
    public = list(
      #' @description constructor for SBM
      #' @param model character describing the type of model
      #' @param nbNodes number of nodes in each dimension of the network
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The dimensions of mu must match \code{blockProp} lengths
      #' @param dimLabels optional labels of each dimension (in row, in column)
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      initialize = function(model, nbNodes, blockProp, connectParam, dimLabels=c(row="row", col="col"), covarParam=numeric(length(covarList)), covarList=list()) {

        ## SANITY CHECKS (on parameters)
        stopifnot(length(dimLabels) == 2)
        stopifnot(length(blockProp) ==  2, is.list(blockProp),
                  length(blockProp[[1]]) ==  nrow(connectParam$mean), # dimensions match between vector of
                  length(blockProp[[2]]) ==  ncol(connectParam$mean)) # block proportion and connectParam$mean
        stopifnot(all(blockProp[[1]] > 0), all(blockProp[[1]] < 1))   # positive proportions
        stopifnot(all(blockProp[[2]] > 0), all(blockProp[[2]] < 1))
        names(blockProp) <- names(dimLabels)

        ## Check that connectivity parameters and model are consistent
        switch(model,
          "bernoulli"  = stopifnot(all(connectParam$mean >= 0), all(connectParam$mean <= 1)),
          "poisson"    = stopifnot(all(connectParam$mean >= 0)),
          "gaussian"   = stopifnot(length(connectParam$var) == 1 | length(connectParam$var) == length(connectParam$mean),
                                   connectParam$var > 0),
          "ZIgaussian" = stopifnot(length(connectParam$var) == 1 | length(connectParam$var) == length(connectParam$mean),
                                   connectParam$var > 0, all(connectParam$p0 >= 0), all(connectParam$p0 <= 1))
        )

        super$initialize(model, NA, nbNodes, dimLabels, blockProp, connectParam, covarParam, covarList)
      },
      #' @description a method to sample new block memberships for the current SBM
      #' @param store should the sampled blocks be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled blocks
      rMemberships = function(store = FALSE) {
        Z <- list(
          row = t(rmultinom(private$dim[1], size = 1, prob = private$pi[[1]])),
          col = t(rmultinom(private$dim[2], size = 1, prob = private$pi[[2]]))
          )
        if(!is.null(private$dimlab)){names(Z) <- private$dimlab}
        if (store) private$Z <- Z
        Z
      },
      #' @description a method to sample a network data (edges) for the current SBM
      #' @param store should the sampled edges be stored (and overwrite the existing data)? Default to FALSE
      #' @return the sampled network
      rEdges = function(store = FALSE) {
        Y <- private$sampling_func[[1]](private$dim[1]*private$dim[2], list(mean = self$expectation, var = private$theta$var)) %>%
          matrix(private$dim[1], private$dim[2])
        if (store) private$Y <- Y
        Y
      },
      #' @description prediction under the current parameters
      #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated.
      #' @param theta_p0 double for thresholding...
      predict = function(covarList = self$covarList, theta_p0 = 0) {
        stopifnot(!is.null(private$Z[[1]]),
                  !is.null(private$Z[[2]]),
                  !is.null(private$theta$mean))
        stopifnot(is.list(covarList),  self$nbCovariates == length(covarList))

        mu <- ((1-theta_p0)>0.5 ) * private$theta$mean
        if (length(covarList) > 0) {
          stopifnot(all(sapply(covarList, nrow) == self$nbNodes[1]),
                    all(sapply(covarList, ncol) == self$nbNodes[2]))
          if (self$modelName == "bernoulli") {
            res <- private$invlink[[1L]](private$Z[[1]] %*% private$link[[1L]]( mu ) %*% t(private$Z[[2]]) + self$covarEffect)
          } else {
            res <- private$invlink[[1L]](private$link[[1L]](private$Z[[1]] %*% mu %*% t(private$Z[[2]])) + self$covarEffect)
          }

        } else {
          res <- private$Z[[1]] %*% mu %*% t(private$Z[[2]])
        }
        rownames(res)<- rownames(private$Y)
        colnames(res)<- colnames(private$Y)
        res

      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Bipartite Stochastic Block Model") {super$show(type)},
      #' @description basic matrix plot method for BipartiteSBM object or mesoscopic plot
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
      #' @param plotOptions list with the parameters for the plot. See help of the corresponding S3 method for details.
      #' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g}, the \code{layout} and the \code{plotOptions} for the \code{'meso'}
      #' @import ggplot2
      plot = function(type = c('data','expected','meso'), ordered = TRUE, plotOptions = list()) {

        if(is.null(self$memberships)) {ordered <- FALSE; type <- 'data'}
        if (ordered & !is.null(self$memberships))
          clustering <- setNames(self$memberships, c('row', 'col'))
        else
          clustering <- NULL

        switch(match.arg(type),
          "meso" =
            plotMeso(
              thetaMean  = private$theta$mean,
              pi         = private$pi,
              model      = private$model,
              directed   = private$directed_,
              bipartite  = TRUE,
              nbNodes    = self$nbNodes,
              nodeLabels = as.list(private$dimlab),
              plotOptions),
          "data" =
            plotMatrix(self$networkData,
                       private$dimlab,
                       clustering,
                       plotOptions),
          "expected" =
            plotMatrix(self$expectation,
                       private$dimlab,
                       clustering,
                       plotOptions)
        )
      }
    ),
    active = list(
### field with write access
      #' @field dimLabels vector of two characters giving the label of each connected dimension (row, col)
      dimLabels    = function(value) {
        if (missing(value))
          return(private$dimlab)
        else {
          stopifnot(is.atomic(value), is.character(value), length(value) == 2)
          if (is.null(names(value))){names(value)  = c('row', 'col')}
          private$dimlab <- value
        }
      },
      #' @field blockProp list of two vectors of block proportions (aka prior probabilities of each block)
      blockProp   = function(value) {
        if (missing(value))
          return(private$pi)
        else {
          stopifnot(is.list(value), length(value) == length(private$dimlab))
          walk(value, ~stopifnot(is.numeric(.x), all(.x > 0), all(.x < 1)))
          private$pi <- setNames(value, private$dimlab)
        }
      },
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam   = function(value) {
        if (missing(value))
          return(private$theta)
        else {
          stopifnot(is.list(value))
          ## Check that connectivity parameters and model are consistent
          switch(private$model,
            "bernoulli"  = stopifnot(all(value$mean >= 0), all(value$mean <= 1)),
            "poisson"    = stopifnot(all(value$mean >= 0)),
            "gaussian"   = stopifnot(all(value$var > 0)),
            "ZIgaussian" = stopifnot(all(value$p0 >= 0), all(value$p0 <= 1))
          )
          private$theta <- value
        }
      },
      #' @field probMemberships  matrix of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {
        if (missing(value))
          return(private$Z)
        else {
          stopifnot(is.list(value), length(value) == length(private$dimlab))
          walk2(value, private$dim, ~stopifnot(nrow(.x) == .y))
          private$Z <- value
        }
      },
### field with access only
      #' @field nbBlocks vector of size 2: number of blocks (rows, columns)
      nbBlocks = function(value) {if(!is.null(private$Z)) setNames(map_int(private$pi, length), private$dimlab)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {private$dim[1] * private$dim[2]},
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {sum(map_int(private$theta, ~length(.x)))},
      #' @field memberships list of size 2: vector of memberships in row, in column.
      memberships = function(value) {
        if (!is.null(private$Z)) return(setNames(map(private$Z, as_clustering), private$dimlab))
      },
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {map(private$Z, ~as_indicator(as_clustering(.x)))}
    )
  )

