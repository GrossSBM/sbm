#' R6 Class definition of a Multipartite SBM fit
#'
#' This class is designed to give a representation and adjust a Multiparite SBM fitted with GREMLIN.
#'
#' @import R6 GREMLIN
#' @export
MultipartiteSBM_fit <-
  R6::R6Class(classname = "MultipartiteSBM_fit",
              inherit = MultipartiteSBM,
              # fields for internal use (referring to the mathematical notation)
              private = list(
                listFit = NULL, #objets estimes une liste d'une classe a creeer (simpleSBM ou bipartiteSBM...)
                # listFit serait cree et rempli par optimize
                # contains outputs of GREMLIN
                GREMLINobject       = NULL,
                import_from_GREMLIN = function(index = 1) {
                  GREMLINfit <- private$GREMLINobject$fittedModel[[index]]
                  clusterings <- lapply(GREMLINfit$paramEstim$tau,as_clustering)
                  lapply(private$listNet,function(net){
                    if (substr(class(net)[1],1,6)=="Simple") {
                      net$memberships <- clusterings[[net$dimLabels[[1]]]]
                      print(clusterings[[net$dimLabels[[1]]]])
                      print(net$memberships)
                    }
                    else {
                      net$memberships <- list(clusterings[[net$dimLabels[[1]]]],clusterings[[net$dimLabels[[2]]]])
                    }
                    #net$tau <- private$GREMLINobject$fittedModel[[1]]$paramEstim

                  })
                }
              ),
              public = list(
                #' @description constructor for Multiparite SBM
                #' @param listSBM list of SBM object with
                initialize = function(listSBM) {
                  super$initialize(listSBM)
                },
                #' @description estimation via GREMLIN
                #' @param options options for MultipartiteBM
                optimize = function(currentOptions) {
                  # formatting data for using GREMLIN
                  listNetG <- lapply(private$listNet,function(net){
                      if (substr(class(net)[1],1,6)=="Simple") {
                        ifelse(net$directed,type<-"diradj",type<-"adj")
                      }
                        else {type <-  "inc"}
                      GREMLIN::defineNetwork(net$netMatrix,type,rowFG=net$dimLabels[[1]],colFG=net$dimLabels[[2]])
                  })
                 # print(length(private$listNet))
                #  print(length(listNetG))
                  vdistrib <- sapply(private$listNet,function(net){
                    net$modelName
                  })

                 # print(vdistrib)
                  print(listNetG[[1]])
                  print(listNetG[[2]])

                  private$GREMLINobject <- GREMLIN::multipartiteBM(list_Net=listNetG,v_distrib=vdistrib)

                  private$import_from_GREMLIN()
                  ### TODO find what to ouput???
                  #import_from_GREMLIN doit remplir listFit
                },
                getBM = function(i) {private$listNet[[i]]} # a supprimer ?
              )
              # ,
              # active=list(
              #   getBM = function(i) {private$listNet[[i]]}
              # )
  )

