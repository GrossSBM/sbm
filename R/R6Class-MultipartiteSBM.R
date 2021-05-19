#' R6 Class definition of a Multipartite SBM
#'
#' @export
MultipartiteSBM <-
  R6::R6Class(
    classname = "MultipartiteSBM",
    inherit = SBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      arch = NULL # matrix describing the organization of the multipartite network
    ),
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param model character describing the type of model
      #' @param architecture a 2-column matrix describing interactions between the networks
      #' @param directed vector of logical: are the network directed or not?
      #' @param nbNodes number of nodes in each dimension/part of the network
      #' @param dimLabels labels of each par of the network
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam parameters of connectivity (vector of list of vectors)
      initialize = function(model = character(0), architecture = matrix(NA, 0, 2), directed = logical(0),
                            nbNodes = numeric(0), dimLabels = character(0), blockProp=list(), connectParam=list()) {

        ## SANITY CHECK
        stopifnot(is.matrix(architecture), ncol(architecture) == 2)
        stopifnot(is.logical(directed), nrow(architecture) == length(directed))

        ## Check that connectivity parameters and model are consistent
        walk2(model, connectParam,
          ~switch(.x,
            "bernoulli"  = stopifnot(all(.y$mean >= 0), all(.y$mean <= 1)),
            "poisson"    = stopifnot(all(.y$mean >= 0)),
            #"gaussian"   = stopifnot(length(.y$var) == 1, .y$var > 0),
            "gaussian"   = stopifnot(.y$var > 0),
            "ZIgaussian" = stopifnot(all(.y$p0 >= 0), all(.y$p0 <= 1))
          )
        )

        ## MODEL & PARAMETERS
        super$initialize(model, directed, nbNodes, dimLabels, blockProp, connectParam)
        private$arch <- architecture
      },
      #' @description print method
      #' @param type character to tune the displayed name
      show = function(type = "Multipartite Stochastic Block Model"){
        cat(type, "\n")
        cat(length(self$dimLabels), "parts/functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
        cat("=====================================================================\n")
        cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
        cat("distributions on each network: ", self$modelName ,"\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $architecture \n")
        cat("  $modelName, $blockProp, $connectParam, $memberships, $networkData\n")
        cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
        cat("* R6 and S3 methods \n")
        cat("  plot, print, coef, predict, fitted, $setModel, $reorder \n")
      },
      #' @description print method
      print = function() self$show(),
      #' @description plot Multipartite Network
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
      #' @param plotOptions list of plot options for the mesoscopic view or matrix view
      plot = function(type = c('data','expected','meso'), ordered = TRUE, plotOptions = list()){

        if (ordered) clustering <- self$memberships else clustering <- NULL
        switch(match.arg(type),
          "meso" =
            plotMesoMultipartite(
              private$arch, self$connectParam, private$pi, private$model,
              private$directed_, private$dim, private$dimlab, plotOptions
            ),
          "data" =
            plotMultipartiteMatrix(
              map(private$Y,"networkData"),
              private$arch, private$dim, private$dimlab, NULL,
              private$model, clustering, plotOptions
            ),
          "expected" =
            plotMultipartiteMatrix(
              self$predict(),
              private$arch, private$dim, private$dimlab, NULL,
              private$model, clustering, plotOptions
            )
        )
      }
    ),
    active = list(
      #' @field dimLabels vector of characters giving the label of each connected dimension
      dimLabels    = function(value) {
        if (missing(value))
          return(private$dimlab)
        else {
          stopifnot(is.atomic(value), is.character(value))
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
          walk2(private$model, value,
            ~switch(.x,
              "bernoulli"  = stopifnot(all(.y$mean >= 0), all(.y$mean <= 1)),
              "poisson"    = stopifnot(all(.y$mean >= 0)),
              "gaussian"   = stopifnot(length(.y$var) == 1, .y$var > 0),
              "ZIgaussian" = stopifnot(all(.y$p0 >= 0), all(.y$p0 <= 1))
            )
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
      #' @field nbBlocks : vector with the number of blocks in each FG
      nbBlocks = function(value) {if(!is.null(private$Z)) setNames(map_int(private$pi, length), private$dimlab)},
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {sum(map_int(private$theta, ~map_int(.x, length)))},
      #' @field architecture organization of the multipartite network
      architecture = function(value) {private$arch},
      #' @field nbNetworks number of networks in the multipartite network
      nbNetworks = function(value) {length(private$directed_)},
      #' @field memberships list of size 2: vector of memberships in all parts of the network
      memberships = function(value) {if (!is.null(private$Z)) map(private$Z, as_clustering)},
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {map(private$Z, ~as_indicator(as_clustering(.x)))}
    )
  )
