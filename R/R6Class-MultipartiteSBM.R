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
       listNet = NULL,
       E = NULL,
       nbFG = NULL,
       dimFG = NULL,
       namesFG = NULL,
       allZ = NULL
       ),
     public = list(
       #' @description constructor for Multipartite SBM
       #' @param listSBM list of SimpleSBM or BipartiteSBM
       initialize = function(listSBM, memberships = NULL) {
         private$listNet <- listSBM
         private$nbNet <- length(listSBM)
         private$namesFG <- unique(unlist(lapply(listSBM, function(net){net$dimLabels})))
         private$nbFG <- length(private$namesFG)
         E_FG <- lapply(listSBM,function(net){return(c(net$dimLabels$row,net$dimLabels$col))})
         E_FG <- do.call(rbind,E_FG)
         E <- matrix(sapply(E_FG,function(a){which(private$namesFG == a)}),private$nbNet,2,byrow = FALSE)
         private$E <-  E
         private$dimFG <- sapply(1:private$nbFG ,function(k){
           u <- which(E[,1] == k); v = 1;
           if (length(u) == 0){u <- which(E[,2] == k); v = 2}
           u <- u[1]
           dim(listSBM[[u]]$netMatrix)[v]})
         private$allZ <- memberships
         },
       plot = function(normMat){
         listNetMatrix = lapply(private$listNet,function(s){s$netMatrix})
         plotMultipartite(listNetMatrix, private$E, private$dimFG, private$namesFG,normalizing = normMat, clustering = private$allZ)
       }
     ),
     active = list(
         #' @field nbNetworks : number of networks in the multipartite network
         nbNetworks    = function(value) {private$nbNet},
         #' @field listSBM : list of SimpleSBMs or BipartiteSBMs
         listSBM    = function(value) {private$listNet},
         #' @field archiMultipartite : organisation of the multipartite network
         archiMultipartite     = function(value) {private$E},
         #' @field dimLabels  : labels of the functional groups
         dimLabels   = function(value){private$namesFG},
         #' @field nbLabels  : number of Functional groups involved in the multipartite
         nbLabels   = function(value){private$nbFG},
         #' @field nbNodes  : number of Nodes in each FG,
         nbNodes  = function(value){private$dimFG}
     )
     )
