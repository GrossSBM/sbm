#' R6 Class definition of a Multiplex SBM fit
#'
#' This class is designed to give a representation and adjust a Multiplex SBM fitted with GREMLIN.
#'
#' @import R6 GREMLINS
#' @export
MultiplexSBM_fit <-
  R6::R6Class(
    classname = "MultiplexSBM_fit",
    inherit = MultipartiteSBM_fit,
    # fields for internal use (referring to the mathematical notation)
    private = list(
      BMobject  = NULL,
      dependent = NULL,
      names_layers_  = NULL,
      import_from_BM  = function(index = which.max(private$BMobject$ICL)) {
        private$J     <- private$BMobject$PL[index]
        private$vICL  <- private$BMobject$ICL[index]
        parameters    <- private$BMobject$model_parameters[[index]]
        private$theta <- switch(private$BMobject$model_name,
                "gaussian_multivariate" = list(mean=parameters$mu,cov=parameters$Sigma),
                "bernoulli_multiplex"   = list(prob00=parameters$pi$`00`,prob01=parameters$pi$`01`,
                                               prob10=parameters$pi$`10`,prob11=parameters$pi$`11`)
        )},
      import_from_BM_Simple = function(index = which.max(private$BMobject$ICL)) { # a function updating the Class
        private$import_from_BM(index)
        listMemberships <- list(private$BMobject$memberships[[index]]$Z)
        names(listMemberships) <- self$dimLabels
        private$Z <- listMemberships
        private$pi  <- colMeans(private$Z[[1]])
      },
      import_from_BM_Bipartite  = function(index = which.max(private$BMobject$ICL)) {
        private$import_from_BM(index)
        listMemberships <- list(
           private$BMobject$memberships[[index]]$Z1,
           private$BMobject$memberships[[index]]$Z2
        )
        names(listMemberships) <- self$dimLabels
        private$Z <- listMemberships
        private$pi  <- lapply(private$Z, colMeans)
      }
    ),
    #-----------------------------------------------
    public = list(
      #' @description constructor for Multiplex SBM
      #' @param netList list of SBM object with
      #' @param dependentNet boolean indicating whether dependence is assumed between networks beyond the common dependence on the latent variables
      initialize = function(netList, dependentNet = FALSE) {

        # check whether the multipartite at hand is actually a multiplex
        lab_per_col <- map(netList, "dimLabels") %>%  reduce(rbind) %>% as.data.frame() %>% dplyr::summarize(across(everything(), ~length(unique(.x))))
        if (any(lab_per_col > 1))
          stop("list of networks provided does not correspond to a Multiplex architecture")
        super$initialize(netList)
        # CHECKING dependence structure


        if (dependentNet) {


          if (length(lab_per_col) == 1){
            if (! ( all(self$directed == TRUE) | all(self$directed == FALSE)) )
              stop("in the dependent case, all networks should be either directed or not directed")
          }

          dBern  <- isTRUE(all.equal(self$modelName, rep("bernoulli", self$nbNetworks)))
          dGauss <- isTRUE(all.equal(self$modelName, rep("gaussian" , self$nbNetworks)))
          if (!(dGauss | (dBern&self$nbNetworks == 2)))
            stop("dependency in multiplex network is only handled for Gaussian distribution or a bivariate Bernoulli distribution")
        }
        private$dependent <- dependentNet
        # namesLayers
        private$names_layers_ <- names(netList)
      },
      #' @description estimation of multipartiteSBM via GREMLINS
      #' @param estimOptions options for MultipartiteBM
      #' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
      #'  \itemize{
      #'  \item{"nbCores"}{integer for number of cores used.  Default is 2}
      #'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
      #'  \item{"nbBlocksRange"}{List of length the number of functional groups, each element supplying the minimal and maximal number of blocks to be explored. The names of the list must be the names of the functional groups.  Default value is from 1 to 10)}
      #'  \item{"initBM"}{Boolean. True if using simple and bipartite SBM as initialisations. Default value  = TRUE}
      #'  \item{"maxiterVEM"}{Number of max. number of iterations in  the VEM. Default value  = 100}
      #'  \item{"maxiterVE"}{Number of max. number of iterations in  the VE. Default value  = 100}
      #'}
      optimize = function(estimOptions) {

        if (self$dependentNetwork == FALSE) {
          super$optimize(estimOptions)
        } else {
          currentOptions <- list(
            verbosity     = 3,
            plot          = TRUE,
            explorFactor  = 1.5,
            nbBlocksRange = c(4,Inf),
            nbCores       = 2,
            fast          = TRUE
          )
          currentOptions[names(estimOptions)] <- estimOptions

          ## Transform estimOptions to a suited for blockmodels list of options
          blockmodelsOptions <- list(
            verbosity          = currentOptions$verbosity,
            plotting           = if(currentOptions$plot) character(0) else "",
            explore_min        = currentOptions$nbBlocksRange[1],
            explore_max        = currentOptions$nbBlocksRange[2],
            ncores             = currentOptions$nbCores,
            exploration_factor = currentOptions$explorFactor
          )

          if (self$modelName[1]=="bernoulli") {model_type="bernoulli_multiplex"}
          if (self$modelName[1]=="gaussian") {model_type="gaussian_multivariate"}
          ## generating arguments for blockmodels call

          # membership type
          if (length(self$dimLabels)>1) {
            membership <-  "LBM" ; type="bipartite"
          } else {
            membership <- ifelse(!self$directed[1], "SBM_sym", "SBM") ; type="simple"
          }

          # recuperer les matrices
          Ys <- .na2zero(lapply(self$networkData,function(net) net$networkData))

          args <- list(membership_type =  membership, adj = Ys)
          args <- c(args, blockmodelsOptions)
          private$BMobject <- do.call(paste0("BM_", model_type), args)
          ## performing estimation
          private$BMobject$estimate()
          ## Exporting blockmodels output to simpleSBM_fit fields
          if (type=="simple") private$import_from_BM_Simple() else private$import_from_BM_Bipartite()

          invisible(private$BMobject)

        }
      },
      #' @description plot Multiplex Network
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection). Default to 'data'.
      #' @param ordered TRUE is the matrices are plotted after reorganization with the blocks. Default value = TRUE
      #' @param plotOptions list of plot options for the matrix view
      plot = function(type = c('data','expected'), ordered = TRUE, plotOptions = list()){

        if (ordered) clustering <- self$memberships else clustering <- NULL

        switch(match.arg(type),
               "data" =
                 plotMultipartiteMatrix(
                   map(private$Y,"networkData"),
                   private$arch, private$dim, private$dimlab,private$names_layers_,
                   private$model, clustering, plotOptions
                 ),
               "expected" ={
                  if (private$dependent==F) {expectations = self$predict()}
                 else {
                   expectations=self$predict()
                 }
                  plotMultipartiteMatrix(
                   expectations,
                   private$arch, private$dim, private$dimlab,private$names_layers_,
                   private$model, clustering, plotOptions

                 )}
        )
      },
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Multiplex Stochastic Block Model"){
          cat(type, "\n")
          type <- ifelse(length(self$nbNodes==1),'Simple','Bipartite')
          cat(type, ' - ',self$nbNetworks, "layers networks", self$namesLayers ,'\n' )
          cat("nodesNames:",self$dimLabels,"\n")
          cat("=====================================================================\n")
          cat("nbNodes = (", self$nbNodes, ") --  nbBlocks  = (",self$nbBlocks, ")\n")
          cat("distributions on each network: ", self$modelName ,"\n")
          cat("Dependent networks:", self$dependentNetwork,"\n")
          cat("=====================================================================\n")
          cat("* Useful fields \n")
          cat("  $nbNetwork, $nbNodes, $nbBlocks, $dimLabels, $dependentNetwork \n")
          cat("  $modelName, $blockProp, $connectParam, $memberships, $networkData\n")
          cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
          cat("* R6 and S3 methods \n")
          cat("  plot, print, coef, predict, fitted, $setModel, $reorder \n")
        }

  ,
  #' @description prediction under the currently estimated model
  #' @return a list of matrices matrix of expected values for each dyad
  predict = function() {
    if (private$dependent) {
      warning("provided expectations are the marginal expectations in the dependent case")

      if (self$modelName[1]=="bernoulli")
      {
        lmu <- list(self$connectParam$prob10+self$connectParam$prob11,
                   self$connectParam$prob01+self$connectParam$prob11)
      }
      if (self$modelName[1]=="gaussian")
      {
        lmu <- self$connectParam$mean
      }

      if (length(self$dimLabels)==1)
      {Z <- self$probMemberships[[1]]
         pred <- lapply(lmu,function(mu){Z %*% mu %*% t(Z)})
      }
      else {
        Z1 <-  self$probMemberships[[1]]
        Z2 <-  self$probMemberships[[2]]
        pred <- lapply(lmu,function(mu){Z1 %*% mu %*% t(Z2)})
      }
    }
    else pred = map(private$Y, predict)
    pred }
  ),
  active = list(
    ### field with access only
    #' @field nbBlocks vector of size 2: number of blocks (rows, columns)
    nbBlocks    = function(value) {
      if(is.list(private$pi)){map_int(private$pi, length)}else{length(private$pi)}},
    #' @field dependentNetwork : connection parameters in each network
    dependentNetwork = function(value) {private$dependent},
    #' @field storedModels data.frame of all models fitted (and stored) during the optimization
    storedModels = function(value) {

      type <- ifelse(length(self$nbNodes)==1,'Simple','Bipartite')
      dep <- self$dependentNetwork
      if (dep){

        nbConnectParam <- c(unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters)))
        if(type == 'Simple'){
          nbBlocks <- unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z)))
          nbConnectParam <- unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters))
          U <- data.frame(
            indexModel  = 1:length(nbBlocks),
            nbParams = nbConnectParam + nbBlocks - 1 -1 * (self$modelName[1]=='gaussian'),
            nbBlocks = nbBlocks,
            ICL      = private$BMobject$ICL,
            loglik   = private$BMobject$PL
          )
          }else{
            rowBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z1))))
            colBlocks <- c(0, unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z2))))
            nbConnectParam <- c(NA, unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters)))
            U <- data.frame(
              indexModel = rowBlocks + colBlocks,
              nbParams  = nbConnectParam  + rowBlocks + colBlocks - 2  -1 * (self$modelName[1]=='gaussian'),
              rowBlocks = rowBlocks,
              colBlocks = colBlocks,
              nbBlocks  = rowBlocks + colBlocks,
              ICL       = private$BMobject$ICL,
              loglik    = private$BMobject$PL
            )
            U[!is.na(U$nbParams), , drop = FALSE]
        }
      }else{
        fit <- private$GREMLINSobject$fittedModel
        nbModels <- length(fit)
        if(type == 'Simple'){
          nbBlocks <- vapply(fit, function(m) m$paramEstim$v_K,1)
          }else{
          Blocks_ <- as.data.frame(t(sapply(fit, function(m) m$paramEstim$v_K)))
          rowBlocks  <- Blocks_[,1]
          colBlocks <- Blocks_[,2]
        }
        nbConnectParam <- sapply(fit, function(m){
          computeNbConnectParams_MBM(m$paramEstim$v_K, private$model, private$arch, private$directed_)
        })


        U <- switch(type,
                    Simple=
          data.frame(
          indexModel =  1:nbModels,
          nbParams  = nbConnectParam  +   nbBlocks-1,
          nbBlocks  = nbBlocks,
          ICL       = map_dbl(fit, "ICL"),
          loglik    = map_dbl(fit, ~last(.x$vJ))),
          Bipartite =
            data.frame(
            indexModel =  1:nbModels,
            nbParams  = nbConnectParam  +  rowBlocks +  colBlocks- 2 ,
            rowBlocks  = rowBlocks,
            colBlocks =  colBlocks,
            nbBlocks  = rowBlocks + colBlocks,
            ICL       = map_dbl(fit, "ICL"),
            loglik    = map_dbl(fit, ~last(.x$vJ)))
          )

      }

    U[order(U$ICL,decreasing = TRUE),]
    },
    #########
    #' @field namesLayers : names of the various Networks
    namesLayers    = function(value) {
      if (missing(value))
        return(private$names_layers_)
      else {
        stopifnot(length(value)==private$model)
        private$names_layers_ <- value
      }
    }
  )
)
