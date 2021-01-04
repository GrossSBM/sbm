#' R6 Class definition of a Multiplex SBM
#'
#' R6 virtual class for Multiplex SBM representation
#'
#' @export
MultiplexSBM <-
  R6::R6Class(classname = "MultiplexSBM",
              # fields for internal use (referring to the mathematical notation)
              private = list(
                nbNet = NULL,
                listNet = NULL,
                E = NULL,
                nbFG = NULL,
                dimFG = NULL,
                namesFG = NULL,
                allZ = NULL,
                directed_  = NULL,
                distrib = NULL,
                pi  = NULL, # list of vectors of parameters for block prior probabilities
                theta = NULL,
                tau = NULL
              ),
              public = list(
                #' @description constructor for Multiplex SBM
                #' @param listSBM list of SimpleSBM or BipartiteSBM
                #' @param memberships list of memberships for each node in each function group.Default value is NULL
                initialize = function(listSBM, memberships = NULL) {
                  private$nbNet <- length(listSBM)
                  private$listNet <- listSBM
                  private$namesFG <- unique(unlist(lapply(listSBM, function(net){net$dimLabels})))
                  private$nbFG <- length(private$namesFG)
                  E_FG <- lapply(listSBM,function(net){return(c(net$dimLabels$row,net$dimLabels$col))})
                  E_FG <- do.call(rbind,E_FG)
                  E <- matrix(sapply(E_FG,function(a){which(private$namesFG == a)}),private$nbNet,2,byrow = FALSE)
                  private$E <-  E
                  private$dimFG <- sapply(1:private$nbFG ,function(k){
                    u <- which(E[,1] == k); v = 1;
                    if (length(u) == 0) {u <- which(E[,2] == k); v = 2}
                    u <- u[1]
                    dim(listSBM[[u]]$netMatrix)[v]})
                  private$allZ <- memberships
                  private$distrib <- sapply(listSBM, function(net) {net$modelName})
                  private$directed_ <- sapply(listSBM, function(net) {if(is.null(net$directed)){return(NA)}else{return(net$directed)}})
                },
                #' @description print method
                #' @param type character to tune the displayed name
                show = function(type = "Multiplex Stochastic Block Model"){
                  cat(type, "\n")
                  cat(self$nbLabels, "functional groups (", self$dimLabels, "), ", self$nbNetworks, "networks\n")
                  cat("=====================================================================\n")
                  cat("nbNodes per FG = (", self$nbNodes, ") --  nbBlocks per FG = (",self$nbBlocks, ")\n")
                  cat("distributions on each network =(", self$modelName ,")\n")
                  cat("=====================================================================\n")
                  cat("* Useful fields \n")
                  cat(" $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $archiMultiplex \n")
                  cat(" $modelName, $blockProp, $connectParam, $memberships, \n")
                  cat("* Useful functions \n")
                  cat("$plot, $optimize \n")
                },
                #' @description print method
                print = function() self$show(),
                #' @description plot Multiplex Network
                #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
                #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
                #' @param plotOptions list of plot options for the mesoscopic view or matrix view
                plot = function(type=c('data','expected','meso'), ordered = TRUE, plotOptions = list()){

                  if (length(type)>1){type='data'}
                  if (type %in% c('data','expected')){
                    listNetMatrix = switch(type,
                                           'data'= lapply(private$listNet,function(s){s$netMatrix}),
                                           'expected' = self$predict()
                    )
                    if (ordered) { clust = private$allZ }else{ clust = NULL}
                    distrib <- private$distrib
                    if(type == 'expected'){distrib[1]='notbernoulli'}
                    outP <- plotMultiplexMatrix(listNetMatrix,
                                                   private$E,
                                                   private$dimFG,
                                                   private$namesFG,
                                                   distrib  = distrib,
                                                   clustering = clust,
                                                   plotOptions = plotOptions)
                  }
                  if (type  == 'meso'){

                    outP <- plotMesoMultiplex(private$E,private$theta, private$pi,private$distrib,private$directed_,private$dimFG,private$namesFG ,plotOptions)
                  }
                  outP
                }
              ),
              active = list(
                #' @field nbNetworks : number of networks in the Multiplex network
                nbNetworks    = function(value) {private$nbNet},
                #' @field listSBM : list of SimpleSBMs or BipartiteSBMs
                listSBM    = function(value) {private$listNet},
                #' @field archiMultiplex : organization of the Multiplex network
                archiMultiplex     = function(value) {private$E},
                #' @field dimLabels  : labels of the functional groups
                dimLabels   = function(value){private$namesFG},
                #' @field nbLabels  : number of Functional groups involved in the Multiplex
                nbLabels   = function(value){private$nbFG},
                #' @field nbNodes  : number of Nodes in each FG,
                nbNodes  = function(value){u <- private$dimFG; names(u) = private$namesFG; return(u)},
                #' @field expectation expected values of connection under the currently adjusted model
                expectation = function() {self$predict()},
                #' @field modelName vector of characters, the family of model for the distribution of the edges in each network
                modelName    = function(value) {private$distrib},
                #' @field directed : vector of boolean
                directed  = function(value){private$directed_}

              )
  )
