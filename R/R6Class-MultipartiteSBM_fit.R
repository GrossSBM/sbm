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

                  Res <- multipartiteBM(list_Net=listNetG,v_distrib=vdistrib)
                  Res[[1]]


                }

              )
  )

