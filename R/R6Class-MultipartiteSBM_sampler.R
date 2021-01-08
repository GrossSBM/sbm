#' R6 Class definition of a Multipartite SBM sampler
#'
#' This class is designed to give a representation and sample from a Multipartite SBM
#'
#' @import R6 GREMLINS
#' @export
MultipartiteSBM_sampler <-
  R6::R6Class(
    classname = "MultipartiteSBM_sampler",
    inherit = MultipartiteSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
       Z = NULL # a list of the sampled indicator of blocks
    ),
    public = list(
      #' @description print/show method
      #' @param type character to tune the displayed name
      show = function(type = "Sampler for a Stochastic Block Model") {
        super$show(type)
        cat("  $indMemberships, $expectation\n")
        cat("* S3 methods \n")
        cat("  plot, print, coef\n")
      }
    ),
    active = list(
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {private$Z},
      #' @field expectation expected values of connection under the current model
      expectation = function() {
        # TODO: ensure that element from the lsit are sampelr object (either Bipartite or Simple SBM)
        map(private$netList, expectation)
      }
    )
  )
