#' R6 Class definition of a Multipartite SBM fit
#'
#' This class is designed to give a representation and adjust a Multiparite SBM fitted with GREMLIN.
#'
#' @import R6 GREMLIN
#' @export
MultipartiteSBM_fit <-
  R6::R6Class(
    classname = "MultipartiteSBM_fit",
    inherit = MultipartiteSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      GREMLINobject       = NULL,
      import_from_GREMLIN = function(index = 1) {
        GREMLINfit <- private$GREMLINobject$fittedModel[[index]]
        listtau <- GREMLINfit$paramEstim$tau
        listpi <-  GREMLINfit$paramEstim$list_pi


        lapply(private$listNet, function(net) {
          if (substr(class(net)[1], 1, 6) == "Simple") {
            Lab <- net$dimLabels[[1]]
            net$varProb <- listtau[[Lab]]
            net$blockProp <- listpi[[Lab]]
          }
          else {
            rowLab <- net$dimLabels[[1]]
            colLab <- net$dimLabels[[2]]
            net$varProb <-
              list(listtau[[rowLab]], listtau[[colLab]])
            net$blockProp <-
              list(listpi[[rowLab]], listpi[[colLab]])
          }
        })
        #print(private$nbNet)
        for (i in 1:private$nbNet) {
          private$listNet[[i]]$connectParam <-
            GREMLINfit$paramEstim$list_theta[[i]]
        }
        #print(private$GREMLINfit$paramEstim$list_theta[[i]]);print(private$listNet[[i]]$connectParam)}
      }
    ),
    public = list(
      #' @description constructor for Multiparite SBM
      #' @param listSBM list of SBM object with
      initialize = function(listSBM) {
        super$initialize(listSBM)
      },
      #' @description estimation via GREMLIN
      #' @param currentOptions options for MultipartiteBM
      optimize = function(currentOptions) {


        # ----- formatting data for using GREMLIN
        listNetG <- lapply(private$listNet, function(net) {
          if (substr(class(net)[1], 1, 6) == "Simple") {
            ifelse(net$directed, type <- "diradj", type <- "adj")
          }
          else {
            type <-  "inc"
          }
          GREMLIN::defineNetwork(net$netMatrix,
                                 type,
                                 rowFG = net$dimLabels[[1]],
                                 colFG = net$dimLabels[[2]])
        })

        vdistrib <- sapply(private$listNet, function(net) {
          net$modelName
        })

        v_Kmin  <- sapply(1:private$nbFG, function(k){currentOptions$nbBlocksRange[[k]]})
        v_Kmax  <- sapply(1:private$nbFG, function(k){currentOptions$nbBlocksRange[[k]][2]})
        v_Kinit <- sapply(1:private$nbFG, function(k){currentOptions$nbBlocksRange[[k]]})
        verbose <- (currentOptions$verbosity>0)
        nbCores <- currentOptions$nbCores
        maxiterVE <- currentOptions$maxiterVE
        maxiterVEM <- currentOptions$maxiterVEM
        namesFG <- names(currentOptions$nbBlocksRange)

        private$GREMLINobject <- GREMLIN::multipartiteBM(list_Net = listNetG, v_distrib = vdistrib, initBM = TRUE)
        #multipartiteBM(list_Net,  v_distrib = NULL ,namesFG = NULL, v_Kmin = 1 , v_Kmax = 10 , v_Kinit = NULL , save = TRUE , verbose = TRUE, nbCores = NULL, maxiterVE = NULL , maxiterVEM = NULL)


        private$import_from_GREMLIN()
        ### TODO find what to ouput???
        #import_from_GREMLIN doit remplir listFit
      },
      #' @description getBM returns a given network from a Multipartite SBM
      #' @param i index of the asked network
      getBM = function(i) {
        private$listNet[[i]]
      }
    )
    # ,
    # active=list(
    #   getBM = function(i) {private$listNet[[i]]}
    # )
  )
