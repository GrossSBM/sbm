#' R6 Class definition of a Multipartite SBM
#'
#' R6 virtual class for Multipartite SBM representation
#'
#' @export
MultipartiteSBM <-
  R6::R6Class(classname = "MultipartiteSBM",
     # fields for internal use (referring to the mathematical notation)
     private = list(
       nbNet = NULL,
       listNet = NULL
                       ),
     public = list(
       #' @description constructor for Multiparite SBM
       #' @param listSBM list of SBM object with
       initialize = function(listSBM) {
          private$listNet = listSBM
          private$nbNet = length(listSBM)
       }
     )
     )
