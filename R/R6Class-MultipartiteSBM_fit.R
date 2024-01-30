#' R6 Class definition of a Multipartite SBM fit
#'
#' This class is designed to give a representation and adjust a Multipartite SBM fitted with GREMLIN.
#'
#' @import R6 GREMLINS
#' @export
MultipartiteSBM_fit <-
  R6::R6Class(
    classname = "MultipartiteSBM_fit",
    inherit = MultipartiteSBM,
    # fields for internal use (referring to the mathematical notation)
    private = list(
       J = NULL,
       vICL = NULL,
       GREMLINSobject = NULL,

      #------------ function to convert GREMLINS result into a sbm object result
      import_from_GREMLINS = function(index = 1) {

        fit <- private$GREMLINSobject$fittedModel[[index]]
        private$J    <- last(fit$vJ)
        private$vICL <- last(fit$ICL)

        ## extract pi, theta, tau
        list_pi    <- fit$paramEstim$list_pi[private$dimlab]
        list_theta <- fit$paramEstim$list_theta#[paste0(private$dimlab[private$arch[,1]],
                                            #private$dimlab[private$arch[ ,2]])]
        list_tau   <- fit$paramEstim$tau[private$dimlab]
        list_theta_mean <- map_if(list_theta, is.list, "mean", ~.x)

        #----------------------------------------------------------
        ## Reordering
        ordAll <- order_mbm(list_theta_mean, list_pi, private$arch)
        list_tau <- map2(list_tau, ordAll, ~.x[, .y, drop = FALSE])
        list_pi  <- map2(list_pi, ordAll, ~.x[.y])
        for (s_ in 1:self$nbNetworks){
          o_row <- ordAll[[private$arch[s_, 1]]]
          o_col <- ordAll[[private$arch[s_, 2]]]

          l_s <- list(mean  =  list_theta_mean[[s_]][o_row , o_col, drop = FALSE])
          if (private$model[s_] %in% c('gaussian','ZIgaussian')) {
            var_s <- list_theta[[s_]]$var
            if (is.matrix(var_s)){var_s  = var_s[o_row,o_col]}
            l_s$var <- var_s
          }
          if (private$model[s_] == 'ZIgaussian'){
            p0_s <- list_theta[[s_]]$p0
            if (is.matrix(p0_s)){p0_s  = p0_s[o_row,o_col]}
            l_s$p0 <- p0_s
          }

          private$Y[[s_]]$connectParam <- l_s
        }
        private$theta <- map(private$Y, "connectParam")

        lapply(private$Y, function(net) {
          if (inherits(net, "SimpleSBM_fit")) {
            Lab <- net$dimLabels
            net$probMemberships <- list_tau[[Lab]]
            net$blockProp       <- list_pi[[Lab]]
          }
          if (inherits(net, "BipartiteSBM_fit")) {
            rowLab <- net$dimLabels[[1]]
            colLab <- net$dimLabels[[2]]
            net$probMemberships <-
              list(list_tau[[rowLab]], list_tau[[colLab]])
            net$blockProp <-
              list(list_pi[[rowLab]], list_pi[[colLab]])
          }
        })

        private$Z <- list_tau
        private$pi  <- list_pi

      }
    ),
    #-----------------------------------------------
    public = list(
      #' @description constructor for Multipartite SBM
      #' @param netList list of SBM objects
      initialize = function(netList) {
        directed  <- map(netList, "directed") %>% map_lgl(~ifelse(is.null(.x), NA, .x))
        dimLabels <- map(netList, "dimLabels") %>% unlist()
        nbNodes   <- map(netList, "nbNodes")  %>% unlist()
        dup <- duplicated(dimLabels);
        if (sum(dup) > 0){
          dimLabels <- dimLabels[-which(dup)]
          nbNodes   <- nbNodes[-which(dup)]
        }

        arch      <- map_if(netList, ~inherits(.x, "SimpleSBM_fit"),
                            function(net) setNames(c(net$dimLabels, net$dimLabels), c("from", "to")),
                    .else = function(net) setNames(unname(net$dimLabels), c("from", "to"))) %>% bind_rows() %>%
           map(factor, levels = dimLabels) %>% map_df(as.numeric) %>% as.matrix()
        super$initialize(model        = map_chr(netList, "modelName"),
                         architecture = arch,
                         directed     = directed,
                         nbNodes      = nbNodes,
                         dimLabels    = dimLabels,
                         blockProp    = map(netList, "blockProp"),
                         connectParam = map(netList, "connectParam"))
        private$Y <- netList
      },
      #' @description estimation of multipartiteSBM via GREMLINS
      #' @param estimOptions options for MultipartiteBM
      #' @inherit estimateSimpleSBM details
      optimize = function(estimOptions) {

        currentOptions <- list(
          verbosity     = 1,
          nbBlocksRange = rep(list(c(1, 10)), length(private$dimlab)),
          nbCores       = 2,
          maxiterVE     = 100,
          maxiterVEM    = 100,
          initBM = TRUE
        )
        names(currentOptions$nbBlocksRange) <- private$dimlab
        ## Current options are default expect for those passed by the user
        currentOptions[names(estimOptions)] <- estimOptions

        # ----- formatting data for using GREMLINS
        listNetG <- lapply(private$Y, function(net) {
          if (inherits(net, "SimpleSBM_fit")) {
            GREMLINS::defineNetwork(
              net$networkData, ifelse(net$directed, "diradj", "adj"),
              rowFG = net$dimLabels, colFG = net$dimLabels)
          } else {
            GREMLINS::defineNetwork(
              net$networkData, "inc",
              rowFG = net$dimLabels[[1]], colFG = net$dimLabels[[2]])
          }
        })

        vdistrib <- private$model
        v_Kmin  <- map_dbl(currentOptions$nbBlocksRange, 1)
        v_Kmax  <- map_dbl(currentOptions$nbBlocksRange, 2)
        verbose <- (currentOptions$verbosity > 0)
        nbCores <- currentOptions$nbCores
        maxiterVE <- currentOptions$maxiterVE
        maxiterVEM <- currentOptions$maxiterVEM
        namesFG <- names(currentOptions$nbBlocksRange)
        initBM <- currentOptions$initBM


        if (sum(abs(v_Kmin - v_Kmax)) > 0) {
          private$GREMLINSobject <- GREMLINS::multipartiteBM(
            list_Net = listNetG,
            v_distrib = vdistrib ,
            namesFG = namesFG,
            v_Kmin = v_Kmin  ,
            v_Kmax = v_Kmax ,
            v_Kinit = NULL ,
            initBM = initBM,
            keep = TRUE ,
            verbose = verbose,
            nbCores = nbCores,
            maxiterVE =  maxiterVE ,
            maxiterVEM =  maxiterVEM)
        } else {
          private$GREMLINSobject <- GREMLINS::multipartiteBMFixedModel(
            list_Net = listNetG,
            v_distrib = vdistrib,
            namesFG = namesFG ,
            v_K = v_Kmax,
            classifInit = NULL,
            nbCores = nbCores,
            maxiterVE = maxiterVE,
            maxiterVEM = maxiterVEM,
            verbose = verbose)
        }
        private$import_from_GREMLINS()
      },
      #' @description prediction under the currently estimated model
      #' @return a list of matrices matrix of expected values for each dyad
      predict = function() {map(private$Y, predict)},
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current MultipartiteSBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$GREMLINSobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_GREMLINS(index)
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Multipartite Stochastic Block Model"){
        super$show(type)
        cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
        cat("* R6 and S3 methods \n")
        cat("  plot, print, coef, predict, fitted, $setModel, $reorder \n")
      }
  ),
  #-----------------------------------------------
  active = list(
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik = function(value) {private$J},
    #' @field ICL double: value of the integrated classification log-likelihood
    ICL    = function(value) {private$vICL},
    #' @field storedModels data.frame of all models fitted (and stored) during the optimization
    storedModels = function(value) {
      fit <- private$GREMLINSobject$fittedModel
      nbModels <- length(fit)
      Blocks <- as.data.frame(t(sapply(fit, function(m) m$paramEstim$v_K)))
      colnames(Blocks) <- paste('nbBlocks', private$dimlab)
      nbConnectParam <- sapply(fit, function(m){
        computeNbConnectParams_MBM(m$paramEstim$v_K, private$model, private$arch, private$directed_)
      })
      nbParams  <- nbConnectParam + rowSums(Blocks) - ncol(Blocks)
      indexModel <- 1:nbModels
      U <- cbind(indexModel, nbParams, Blocks)
      U$nbBlocks <-rowSums(Blocks)
      U$ICL      <- map_dbl(fit, "ICL")
      U$loglik   <- map_dbl(fit, ~last(.x$vJ))
      U
    }
  )
)
