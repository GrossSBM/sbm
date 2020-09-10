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
                GREMLINobject       = NULL,
                import_from_GREMLIN = function(index = 1) {
                  GREMLINfit <- private$GREMLINobject$fittedModel[[index]]
                  listtau <- GREMLINfit$paramEstim$tau
                  listpi <-  GREMLINfit$paramEstim$list_pi


                  lapply(private$listNet,function(net){
                    if (substr(class(net)[1],1,6)=="Simple") {
                      Lab <- net$dimLabels[[1]]
                      net$varProb <- listtau[[Lab]]
                      net$blockProp <- listpi[[Lab]]
                    }
                    else {
                      rowLab <- net$dimLabels[[1]]
                      colLab <- net$dimLabels[[2]]
                      net$varProb <- list(listtau[[rowLab]],listtau[[colLab]])
                      net$blockProp <- list(listpi[[rowLab]],listpi[[colLab]])
                    }
                  })
                  #print(private$nbNet)
                  for (i in 1:private$nbNet) {private$listNet[[i]]$connectParam <- GREMLINfit$paramEstim$list_theta[[i]]}
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
                  # print(listNetG[[1]])
                  # print(listNetG[[2]])

                  private$GREMLINobject <- GREMLIN::multipartiteBM(list_Net=listNetG,v_distrib=vdistrib)

                  private$import_from_GREMLIN()
                  ### TODO find what to ouput???
                  #import_from_GREMLIN doit remplir listFit
                },
             getBM = function(i) {private$listNet[[i]]}
              )
              # ,
              # active=list(
              #   getBM = function(i) {private$listNet[[i]]}
              # )
  )

