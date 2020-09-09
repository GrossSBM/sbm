#' R6 virtual class for SimpleSBM representation
#'
#' @import R6
BipartiteSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(classname = "BipartiteSBM",
              inherit = SBM,
              active=list(
                memberships = function(value) {if (missing(value)) return(private$cluster) else private$cluster <- value}
              ))
# TODO initialize
