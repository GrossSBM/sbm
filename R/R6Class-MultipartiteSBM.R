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
            "gaussian"   = stopifnot(length(.y$var) == 1, .y$var > 0),
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
        cat(self$nbLabels, "functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
        cat("=====================================================================\n")
        cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
        cat("distributions on each network: ", self$modelName ,"\n")
        cat("=====================================================================\n")
        cat("* Useful fields \n")
        cat("  $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $architecture \n")
        cat("  $modelName, $blockProp, $connectParam, $memberships, $networkData\n")
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
              private$arch, private$dim, private$dimlab,
              private$model, clustering, plotOptions
            ),
          "expected" =
            plotMultipartiteMatrix(
              self$predict(),
              private$arch, private$dim, private$dimlab,
              private$model, clustering, plotOptions
            )
        )
      }
    ),
    active = list(
      #' @field architecture organization of the multipartite network
      architecture = function(value) {private$arch},
      #' @field nbNetworks number of networks in the multipartite network
      nbNetworks = function(value) {length(private$directed_)},
      #' @field nbLabels number of functional groups involved in the multipartite
      nbLabels  = function(value){length(private$dimlab)}
    )
  )

#' Check  if an object is MultipartiteSBM
#'
#' Auxiliary function to check the given class of an object
#' @param  Robject an R6 object inheriting from class MultipartiteSBM
#' @return TRUE or FALSE
#' @export
is_MultipartiteSBM <- function(Robject) {inherits(Robject,"MultipartiteSBM")}

