#' R6 virtual class for SBM representation (mother class of Simple and Bipartite SBM fit and sampler)
#'
#' @import R6
SimpleSBM <- # this class inherit from SBM and allow to use multipartite as a list of SimpleSBM and BipartiteSBM
  R6::R6Class(classname = "SimpleSBM",
              inherit = SBM)
