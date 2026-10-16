#----------------------------------------------------------------------------------
myRepeat <- function(v,Qrow,Qcol){c(rep(v[1],Qrow),rep(v[2],Qcol))}
#----------------------------------------------------------------------------------
# Opacity of the block colors, growing with the value from `blockAlpha.min` to `blockAlpha.max`.
# The values are rescaled between `limits` (by default [0, max] for non negative values, the range otherwise)
blockAlpha <- function(m, plotOptions, limits = NULL){
  if (is.null(limits)) {limits <- if (all(m >= 0, na.rm = TRUE)) c(0, max(m, na.rm = TRUE)) else range(m, na.rm = TRUE)}
  intensity <- (pmin(pmax(m, limits[1]), limits[2]) - limits[1]) / diff(limits)
  intensity[!is.finite(intensity)] <- 0
  plotOptions$blockAlpha.min + (plotOptions$blockAlpha.max - plotOptions$blockAlpha.min) * intensity
}
# color `col` with opacity `a`, flattened on white
blendWhite <- function(col, a){
  grDevices::rgb(t(255 - a * (255 - grDevices::col2rgb(col))), maxColorValue = 255)
}
#----------------------------------------------------------------------------------
#' @importFrom rlang .data
#' @importFrom utils head modifyList
#' @importFrom prodlim row.match
#' @importFrom reshape2 melt
#----------------------------------------------------------------------------------
plotMatrix = function(Mat, dimLabels, clustering = NULL, plotOptions = list()){



  ###############################################
  # Core plotting options ----------------------------------------------------
  # `blockBackground` (default TRUE) adds a background to each block (row/col
  # groups defined by `clustering`), colored with `block.color`. Opacity grows with
  # the average value of `Mat` inside the block, from `blockAlpha.min` to `blockAlpha.max`.
  # Without clustering, the blocks have no meaning and no background is drawn.
  currentOptions = list(
    line.color = 'red',
    legend = FALSE,
    rowNames = FALSE,
    colNames = FALSE,
    names.size = NULL,
    title = NULL,
    blockBackground = TRUE,
    graphonColors = FALSE,
    block.color = 'red',
    blockAlpha.min = 0.1,
    blockAlpha.max = 0.6,
    legend.title = FALSE,
    legend.position = 'bottom'
  )
  # Overwrite defaults with any user-supplied options
  currentOptions[names(plotOptions)] <- plotOptions
  ################################################

  n1 <- dim(Mat)[1]
  n2 <- dim(Mat)[2]

  if (length(dimLabels) == 1){dimLabels = rep(dimLabels,2); names(dimLabels)  =c('row','col')}
  if (is.null(names(dimLabels))){names(dimLabels) = c('row','col')}



  if (is.null(rownames(Mat))){rownames(Mat) = as.factor(1:n1)}
  if (is.null(colnames(Mat))){colnames(Mat) = as.factor(1:n2)}


  rowFG <- dimLabels[1]
  colFG <- dimLabels[2]

  if ((!currentOptions$rowNames)){rownames(Mat)= 1:n1}
  if ((!currentOptions$colNames)){colnames(Mat)= 1:n2}

  u <- range(c(Mat))

  binary = FALSE
  val <- sort(unique(c(Mat)))
  if (setequal(val ,c(0,1))) {binary = TRUE}

  if (!is.null(clustering)) {
    l <- length(clustering)
    if (l == 1) {
      oRow <- oCol <- order(clustering$row)
      uRow <- cumsum(table(clustering$row)) + 0.5
      uRow <- uRow[-length(uRow)]
      sepRow <- as.data.frame(uRow)
      sepCol <- sepRow
    }
    if (l == 2) {
      oRow <- order(clustering$row)
      oCol <- order(clustering$col)
      uRow <- cumsum(table(clustering$row)) + 0.5
      uRow <- uRow[-length(uRow)]
      sepRow <- as.data.frame(uRow)
      uCol <- cumsum(table(clustering$col)) + 0.5
      uCol <- uCol[-length(uCol)]
      sepCol <- as.data.frame(uCol)
    }
    Mat <- Mat[oRow,oCol]
    names(sepCol) = names(sepRow) = 'sep'
    sepRow = n1 - sepRow

    if (currentOptions$blockBackground) {
      clRow <- clustering$row[oRow]
      clCol <- (if (l == 2) clustering$col else clustering$row)[oCol]
      endRow <- as.vector(cumsum(table(clRow))); startRow <- c(1, head(endRow, -1) + 1)
      endCol <- as.vector(cumsum(table(clCol))); startCol <- c(1, head(endCol, -1) + 1)
      blocks <- expand.grid(k = seq_along(endRow), q = seq_along(endCol))
      blocks$mean <- mapply(function(k, q) {
        mean(Mat[startRow[k]:endRow[k], startCol[q]:endCol[q]], na.rm = TRUE)
      }, blocks$k, blocks$q)
      # discrete axes: column j sits at x = j, row i at y = n1 - i + 1
      blocks$xmin <- startCol[blocks$q] - 0.5
      blocks$xmax <- endCol[blocks$q] + 0.5
      blocks$ymin <- n1 - endRow[blocks$k] + 0.5
      blocks$ymax <- n1 - startRow[blocks$k] + 1.5
      m <- blocks$mean
      blocks$alpha <- blockAlpha(blocks$mean, currentOptions)
    }
  }
  withBackground <- !is.null(clustering) && currentOptions$blockBackground
  if(is.null(currentOptions$line.width)){currentOptions$line.width = mean(dim(Mat))/300}



  FGRow <- rep(dimLabels[1], times = n1)
  FGCol <- rep(dimLabels[2], times = n2)

  melted_Mat <- reshape2::melt(t(Mat))

  names(melted_Mat) <- c('names_col','names_row','link')
  if(is.numeric(melted_Mat$names_col)){melted_Mat$names_col <- as.character(melted_Mat$names_col)}
  if(is.numeric(melted_Mat$names_row)){melted_Mat$names_row <- as.character(melted_Mat$names_row)}
  melted_Mat$index_col = rep(1:n2,n1)
  melted_Mat$index_row = n1 - rep(1:n1,each = n2) + 1
  melted_Mat$names_col <- as.factor(melted_Mat$names_col)
  melted_Mat$names_col <- factor(melted_Mat$names_col, levels = melted_Mat$names_col[melted_Mat$index_row == 1])
  melted_Mat$names_row <- as.factor(melted_Mat$names_row)
  melted_Mat$names_row <- factor(melted_Mat$names_row, levels = melted_Mat$names_row[melted_Mat$index_col == 1][n1:1])

  melted_Mat$FG_row <- as.factor(rep(FGRow, each = n2))
  melted_Mat$FG_col <- as.factor(rep(FGCol,  n1))
  if (binary) { melted_Mat$link <- as.factor(melted_Mat$link)}

  g <- ggplot(data = melted_Mat, aes(y = .data$names_row, x = .data$names_col, fill = .data$link))
  # with a block background, low values are transparent so that the background shows
  low.color <- if (withBackground) "transparent" else "white"
  g <- g + geom_tile()
  if (withBackground) {
    g <- g + geom_rect(data = blocks, inherit.aes = FALSE,
                       aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, alpha = .data$alpha),
                       fill = currentOptions$block.color) +
      scale_alpha_identity()
    # draw the background under the tiles (added after them so the axes stay discrete)
    g$layers <- g$layers[c(2, 1)]
  }
  if (!binary && currentOptions$graphonColors) {
    # same colors as the graphon (used for the expected network)
    limitsGraphon <- if (all(Mat >= 0, na.rm = TRUE)) c(0, max(Mat, na.rm = TRUE)) else u
    g <- g + scale_fill_gradient(low  = blendWhite(currentOptions$block.color, currentOptions$blockAlpha.min),
                                 high = blendWhite(currentOptions$block.color, currentOptions$blockAlpha.max),
                                 limits = limitsGraphon, na.value = "white")
  } else if (!binary) {g <-  g +  scale_fill_gradient(low = low.color, high = "black", limits = u,na.value = "white")}
  if (binary) {g <- g + scale_fill_manual(breaks = c("0", "1"),values = c(low.color, "black"),na.value = "white")}
  g <- g  +  scale_x_discrete(drop = FALSE, expand = c(0, 0)) + scale_y_discrete(drop = FALSE, expand = c(0, 0))
  # size of the names: by default, small enough for the names not to overlap on a
  # 5 inches plot (rows and columns have the same width since aspect.ratio = n1/n2)
  names.size <- currentOptions$names.size
  if (is.null(names.size)) {names.size <- min(9, 250 / max(n1, n2))}
  if (!currentOptions$colNames){
    g <- g + theme(axis.text.x = element_blank())
  }else{
    # vertical names, ending at the axis and centered on their column
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = names.size))
  }
  if (!currentOptions$rowNames){
    g <- g + theme(axis.text.y = element_blank())
  }else{
    g <- g + theme(axis.text.y = element_text(size = names.size))
  }

  # the labels of the dimensions are the axis titles
  g <- g +  labs(x = colFG, y = rowFG) +
    theme(aspect.ratio = n1/n2, axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
          plot.title = element_text(hjust = 0.5))
  if (!currentOptions$legend){g <- g +theme(legend.position = 'none')}else{
    g <- g +theme(legend.position = currentOptions$legend.position)
    if(!currentOptions$legend.title){g <- g+ theme(legend.title = element_blank())}}




  if (!is.null(clustering)) {
    g <- g + geom_vline(data = sepCol,mapping = aes(xintercept = .data$sep),linewidth = currentOptions$line.width, col=currentOptions$line.color)
    g <- g + geom_hline(data = sepRow + 1,mapping = aes(yintercept = .data$sep),linewidth = currentOptions$line.width, col=currentOptions$line.color)
  }
  if (!is.null(currentOptions$title)){g <- g + ggtitle(currentOptions$title) }
  g
  #if (!is.null(fileNameSave)) { ggsave(fileNameSave, width = 20, height = 20, units = "cm") }else{g}
}
##################################################################################
#----------------------------------------------------------------------------------
# Graphon of an SBM: the block connectivity drawn on the unit square, with blocks
# sized by the block proportions. The axes are in number of nodes, so that the aspect
# ratio of the plot reflects the dimensions of the network.
#----------------------------------------------------------------------------------
plotGraphon = function(thetaMean, pi, nbNodes, dimLabels, plotOptions = list()){

  currentOptions = list(
    line.color = 'red',
    line.width = NULL,
    legend = FALSE,
    legend.title = FALSE,
    legend.position = 'bottom',
    title = NULL,
    block.color = 'red',
    blockAlpha.min = 0.1,
    blockAlpha.max = 0.6,
    limits = NULL
  )
  currentOptions[names(plotOptions)] <- plotOptions

  if (!is.list(pi)) {pi <- list(pi, pi)}
  nbNodes <- rep_len(nbNodes, 2)
  dimLabels <- rep_len(dimLabels, 2)
  if (is.null(currentOptions$line.width)) {currentOptions$line.width = mean(nbNodes)/300}

  # block 1 in the top left corner, as in plotMatrix
  endRow <- nbNodes[1] * cumsum(pi[[1]]); startRow <- c(0, head(endRow, -1))
  endCol <- nbNodes[2] * cumsum(pi[[2]]); startCol <- c(0, head(endCol, -1))
  blocks <- expand.grid(k = seq_along(pi[[1]]), q = seq_along(pi[[2]]))
  blocks$value <- as.matrix(thetaMean)[cbind(blocks$k, blocks$q)]
  blocks$xmin <- startCol[blocks$q]
  blocks$xmax <- endCol[blocks$q]
  blocks$ymin <- -endRow[blocks$k]
  blocks$ymax <- -startRow[blocks$k]

  limits <- currentOptions$limits
  if (is.null(limits)) {limits <- if (all(blocks$value >= 0)) c(0, max(blocks$value)) else range(blocks$value)}
  blocks$value <- pmin(pmax(blocks$value, limits[1]), limits[2])

  g <- ggplot(blocks) +
    geom_rect(aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax, fill = .data$value)) +
    scale_fill_gradient(low  = blendWhite(currentOptions$block.color, currentOptions$blockAlpha.min),
                        high = blendWhite(currentOptions$block.color, currentOptions$blockAlpha.max),
                        limits = limits) +
    geom_vline(xintercept = startCol[-1], linewidth = currentOptions$line.width, col = currentOptions$line.color) +
    geom_hline(yintercept = -startRow[-1], linewidth = currentOptions$line.width, col = currentOptions$line.color) +
    coord_fixed(expand = FALSE) +
    labs(x = dimLabels[2], y = dimLabels[1]) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
          plot.title = element_text(hjust = 0.5))

  if (!currentOptions$legend){g <- g + theme(legend.position = 'none')}else{
    g <- g + theme(legend.position = currentOptions$legend.position)
    if (!currentOptions$legend.title){g <- g + theme(legend.title = element_blank())}}
  if (!is.null(currentOptions$title)){g <- g + ggtitle(currentOptions$title)}
  g
}
##################################################################################
#----------------------------------------------------------------------------------
##################################################################################

plotMultipartiteMatrix = function(listMat, E, nbNodes, namesFG,namesLayers, distrib, clustering, plotOptions = list()) {

  #-------------------------------------------
  list_Mat <- listMat
  nbNet <- length(list_Mat)

  ###---------------------- check for multiplex
  wE <- which(duplicated(E) == FALSE)
  if(length(wE) > 0){
    if (is.null(namesLayers)){namesLayers = paste("Layer", 1:nbNet, sep="")}
    for (i in wE){
      u.i <- rowSums(matrix(E[i,],nrow = nrow(E),2,byrow = T) == E)
      p.i <- which(u.i == 2)
      if(length(p.i) > 1){
      test_bipartite <- (E[p.i[1],1]!=E[p.i[1],2])
      FGi <- E[p.i[1],2]
      for (k in 1:length(p.i)){
        namesFG <- c(namesFG,paste(namesFG[FGi],namesLayers[k],sep='. '))
        nbNodes <- c(nbNodes,nbNodes[FGi])
        if (!is.null(clustering)){clustering[[length(clustering) + 1]] <- clustering[[FGi]]}
        if (!is.null(clustering)){names(clustering)<- namesFG}
        names(nbNodes) <- namesFG
        if(test_bipartite & (k==1)){E[p.i[k],2] <- max(E)}else{E[p.i[k],2] <- max(E)+1}
      }
      if(test_bipartite){
        namesFG <- namesFG[-FGi]
        nbNodes <- nbNodes[-FGi]
        if (!is.null(clustering)){
          rm_FGi <- c(1:(FGi-1),(FGi+1):length(clustering))
          clustering <- lapply(rm_FGi,function(l){clustering[[l]]})
          names(clustering) <- namesFG
        }
      }
    }
    }
    }

  nbFG <- length(unique(c(E)))





  #----------------------------------------
  currentOptions = list(line.color  = 'red',legend = FALSE,compact = TRUE, normalized = FALSE,title=NULL)
  currentOptions$legend.title = FALSE
  currentOptions$legend.position ='bottom'
  currentOptions$nodeNames = FALSE
  currentOptions[names(plotOptions)] = plotOptions
  normalized  <- currentOptions$normalized
  reordered <- !is.null(clustering)
  #----------------------------------------------

  if (reordered){
    if (!is.null(names(clustering))){clustering  <- lapply(1:nbFG,FUN  = function(i){clustering[[namesFG[i]]]})}
    nbBlocks <- sapply(clustering,function(s){length(unique(s))})
    list_Mat_reorder <- lapply(1:nbNet, function(l){
      clustering_row <- clustering[[E[l, 1]]]
      clustering_col <- clustering[[E[l, 2]]]
      oRow <- order(clustering_row)
      oCol <- order(clustering_col)
      list_Mat[[l]][oRow, oCol]})
    clustering <- lapply(1:nbFG,FUN = function(i){sort(clustering[[i]])})
    list_Mat <- list_Mat_reorder
  }
  if (normalized){
    list_Mat = lapply(1:nbNet,function(l){
      model <- distrib[l]
      mat  <- list_Mat[[l]]
      if (model == "poisson"){mat = mat/max(mat)}
      range_mat <- max(mat) - min(mat)
      if ( !(model %in% c("bernoulli","poisson"))){
        if (range_mat !=0){
          mat = (mat - min(mat))/(range_mat) + 0.1
        }
        else{if (min(mat)==0){
          mat = mat}
          else{
            mat = mat/mat}
        }
      }

      return(mat)
    }
    )
  }
  binary <- all(unlist(list_Mat) %in% c(0, 1, NA))
  ############## Optimize positions of matrices
  if (currentOptions$compact) {
    TranposColl <-
      lapply(
        0:(2 ^ nbNet / 2 - 1),
        FUN = function(x)
          head(as.integer(intToBits(x)), nbNet)
      )
    EColl <- lapply(
      1:(2 ^ nbNet / 2),
      FUN = function(i) {
        u <- which(TranposColl[[i]] == 1)
        Ex <- E
        if (length(u) > 0) {for (i in u) {Ex[i, 1:2] = E[i, 2:1]}}
        Ex
      }
    )
    GColl <- sapply(
      1:(2 ^ nbNet / 2),
      FUN = function(i) {
        G <- matrix(0, nbFG, nbFG)
        for (j in 1:nbNet) {G[EColl[[i]][j, 1], EColl[[i]][j, 2]] = 1}
        sum(rowSums(G) == 0) + sum(colSums(G) == 0)
      }
    )
    best <- which.max(GColl)
    E <- EColl[[best]]
    Tbest <- TranposColl[[best]]
    w = which(Tbest == 1)
    if (length(w) > 0) {
      for (i in w) {
        list_Mat[[i]] = t(list_Mat[[i]])
      }
    }
  }
  ############################  meta matrix
  # network of matrices
  G <- matrix(0, nbFG, nbFG)
  for (i in 1:nbNet){G[E[i, 1], E[i, 2]] = 1}

  # places in metamatrix
  GRow <- as.numeric(rowSums(G)>=1)
  GCol <- as.numeric(colSums(G)>=1)
  EndFG_row <- cumsum(nbNodes * GRow)
  BegFG_row <- c(0, cumsum(nbNodes * GRow )[-nbFG]) + 1
  EndFG_col <- cumsum(nbNodes * GCol)
  BegFG_col <- c(0, cumsum(nbNodes * GCol )[-nbFG]) + 1


  n1 <- sum(nbNodes*GRow)
  n2 <- sum(nbNodes*GCol)
  MetaMat <- matrix(NA,n1 ,n2 )
  row_Names <- rep(NA,n1)
  col_Names <- rep(NA,n2)
  for (i in 1:nbFG) {
    for (j in 1:nbFG) {
      if (G[i, j] == 1) {
        place_row <-  BegFG_row[i]:EndFG_row[i]
        place_col <-  BegFG_col[j]:EndFG_col[j]
        kij <- which(rowSums(E == matrix(rep(c(i, j), nbNet),ncol = 2,nrow = nbNet, byrow = T )) == 2)
        net_ij <- list_Mat[[kij]]
        MetaMat[place_row, place_col] <- net_ij
        if(currentOptions$nodeNames){
          row_Names[place_row] <-rownames(net_ij)
          col_Names[place_col] <-colnames(net_ij)
        }
      }
    }
  }
  rownames(MetaMat) <- row_Names
  colnames(MetaMat) <- col_Names


  if(is.null(currentOptions$line.width)){currentOptions$line.width = mean(dim(MetaMat))/300}
  FGCol <- rep(namesFG, times = nbNodes*GCol)
  FGRow <- rep(namesFG, times = nbNodes*GRow)

  ############# meltedMat
  melted_Mat = reshape2::melt(t(MetaMat))
  names(melted_Mat) <- c('names_col','names_row','link')
  if (all(is.na(melted_Mat$names_col))){melted_Mat$names_col <- rep(1:n2,n1)}
  if (all(is.na(melted_Mat$names_row))){melted_Mat$names_row <- n1 - rep(1:n1,each = n2) + 1}

  if(is.numeric(melted_Mat$names_col)){melted_Mat$names_col <- as.character(melted_Mat$names_col)}
  if(is.numeric(melted_Mat$names_row)){melted_Mat$names_row <- as.character(melted_Mat$names_row)}

  melted_Mat$index_col = rep(1:n2,n1)
  melted_Mat$index_row = n1 - rep(1:n1,each = n2) + 1
  melted_Mat$names_col <- as.factor(melted_Mat$names_col)
  melted_Mat$names_col <- factor(melted_Mat$names_col, levels = melted_Mat$names_col[melted_Mat$index_row == 1])
  melted_Mat$names_row <- as.factor(melted_Mat$names_row)
  melted_Mat$names_row <- factor(melted_Mat$names_row, levels = melted_Mat$names_row[melted_Mat$index_col == 1][n1:1])
  melted_Mat$FG_row <- as.factor(rep(FGRow, each = n2))
  melted_Mat$FG_col <- as.factor(rep(FGCol,  n1))
  range_melted <- range(c(melted_Mat$link))
  if (binary) { melted_Mat$link <- as.factor(melted_Mat$link)}


  ############# PLOT

  g <- ggplot2::ggplot(melted_Mat, aes(y = .data$names_row, x = .data$names_col, fill = .data$link))
  g <- g +  geom_tile()
  g <- g +  theme(axis.ticks = element_blank(),panel.background = element_rect(fill = "white"))
  g <- g +  labs(x = '', y = '')
  g <- g +  scale_x_discrete(drop = TRUE) + scale_y_discrete(drop = TRUE)

  name_legend <- ifelse(normalized,"Norm. Link","Link")
  if (!binary) {
    g <- g +  scale_fill_gradient(
      low = "white",
      high = "black",
      limits = range_melted,
      na.value = "transparent",
      name = name_legend
    )
  } else{
    g <-
      g + scale_fill_manual(
        breaks = c("0", "1"),
        values = c("white", "black"),
        na.value = "transparent"
      )
  }
  if (!currentOptions$nodeNames){
    g <- g + theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  }else{
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 0))
  }

  if(!currentOptions$legend){
    g <- g + theme(legend.position = 'none')}
  else{
    g <- g +theme(legend.position = currentOptions$legend.position)
    if(!currentOptions$legend.title){g <- g+ theme(legend.title = element_blank())}
  }


  #g <- g + facet_grid(FG_row~ FG_col, scales = 'free')#coord_equal()
  g <- g + facet_grid(FG_row~ FG_col, scales = 'free', space = 'free')

  #g <- g + facet_grid(FG_row~ FG_col)



  ########## separators
  if (reordered){
    separate <- unlist(lapply(1:nbFG, function(l){
      sepColl <- cumsum(table(clustering[[l]]))
      #sepColl <- table(clustering[[l]])
      sepColl <- sepColl[-length(sepColl)]
      sepColl}))
    nbSep <- length(separate)



    separCol <- data.frame(sepCol = c(separate + 0.5),#+ rep(c(0,cumsum(GCol*nbNodes)[-nbFG]),nbBlocks-1) + 0.5),
                           FG_col = rep(namesFG,nbBlocks-1),
                           FG_col_index = rep(1:nbFG,nbBlocks-1))
    separCol <- do.call("rbind", replicate(nbFG^2,separCol ,simplify = FALSE))
    separCol$FG_row <- rep(namesFG,each = nbSep)
    separCol$FG_row_index <-  rep(1:nbFG,each = nbSep)
    testCol <- vapply(1:(nbSep*nbFG),function(i){1*(G[separCol$FG_row_index[i],separCol$FG_col_index[i]]==1)},1)
    separCol <- separCol[testCol==1,]

    separCol$nameSepCol = c()
    if(nrow(separCol)>1){
    for (i in 1:nrow(separCol)){
      separCol$nameSepCol[i] <- melted_Mat$names_col[melted_Mat$index_row==1][separCol$sepCol[i]]}
    }

    separRow <- data.frame(sepRow = 0.5 + rep(nbNodes*GRow,nbBlocks - 1) -(separate),
                           FG_row = rep(namesFG,nbBlocks-1),
                           FG_row_index = rep(1:nbFG,nbBlocks-1))
    separRow <- do.call("rbind", replicate( nbFG,separRow ,simplify = FALSE))
    separRow$FG_col <- rep(namesFG,each = nbSep)
    separRow$FG_col_index <- rep(1:nbFG,each = nbSep)
    testRow <- vapply(1:(nbSep*nbFG),function(i){1*(G[separRow$FG_row_index[i],separRow$FG_col_index[i]]==1)},1)
    separRow <- separRow[testRow==1,]

    g <- g + geom_vline(data= separCol, aes(xintercept =  .data$sepCol), linewidth = currentOptions$line.width, col=currentOptions$line.color)
    g <- g + geom_hline(data= separRow, aes(yintercept = .data$sepRow), linewidth = currentOptions$line.width, col=currentOptions$line.color)
  }

  if (!is.null(currentOptions$title)){g <- g + ggtitle(currentOptions$title) }

  g
}

# default colors of ggplot2 for n discrete values (hue palette)
ggplotColors <- function(n){
  grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[1:n]
}
# colors of the types of nodes in the mesoscopic view: hue palette without the reds,
# so that the nodes are not confused with the connectivities (drawn in red by default)
mesoColors <- function(n){
  if (n == 1) return(grDevices::hcl(h = 230, l = 65, c = 100))
  grDevices::hcl(h = seq(285, 75, length.out = n), l = 65, c = 100)
}
#-----------------------------------------------------------------
# Mesoscopic view of a simple, bipartite or multipartite SBM: the graphon drawn as a network.
# Each block is a node whose area is proportional to the block proportion, colored according to
# its type of nodes (functional group). The connectivity between two blocks (and
# within a block, as a loop) is drawn as an edge whose width and color grow with the connectivity,
# with the same colors as the graphon, on a scale common to all the networks.
#   piList  : list of the block proportions, one element per functional group
#   netList : list of the networks, each one a list with fg (the two functional groups),
#             mean (the matrix of connectivities) and directed (logical, NA for bipartite)
#   defaultLayout : layout used when plotOptions$layout is NULL: 'groups' (one row of blocks per
#             type of nodes) or 'free' (blocks attracted by their connectivity, whatever their type)
#-----------------------------------------------------------------
#' @importFrom graphics plot
plotMeso <- function(piList, netList, nodeLabels, plotOptions, defaultLayout = 'groups'){

  nbFG <- length(piList)
  currentOptions <- list(
    seed = NULL,
    title = NULL,
    layout = NULL,
    vertex.color = mesoColors(nbFG),
    vertex.frame.color = NA,
    vertex.shape = 'circle',
    vertex.size = 1,
    vertex.label.name = NULL,
    vertex.label.color = 'black',
    vertex.label.size = 3.5,
    vertex.label.inside = 0.5,
    withinBlock = 'loop',
    edge.threshold = -Inf,
    edge.width = 1,
    edge.curved = 0.2,
    edge.arrow.size = 1,
    block.color = 'red',
    blockAlpha.min = 0.1,
    blockAlpha.max = 0.6,
    limits = NULL
  )
  currentOptions[names(plotOptions)] <- plotOptions
  for (opt in c('vertex.color', 'vertex.frame.color', 'vertex.shape', 'vertex.size')) {
    currentOptions[[opt]] <- rep_len(currentOptions[[opt]], nbFG)
  }

  #------------------ blocks (nodes of the mesoscopic graph)
  nbBlocks <- lengths(piList)
  offset <- c(0, cumsum(nbBlocks))[1:nbFG]
  nodes <- data.frame(fg = rep(1:nbFG, nbBlocks), q = unlist(lapply(nbBlocks, seq_len)), pi = unlist(piList))
  Q <- nrow(nodes)

  vertex.label <- currentOptions$vertex.label.name
  if (is.null(vertex.label)) {
    if (nbFG > 1) {
      # shortest prefixes of the labels of the dimensions which distinguish them
      dimLab <- unlist(nodeLabels)
      i <- 1
      prefix <- substr(dimLab, 1, i)
      while (anyDuplicated(prefix) && i < max(nchar(dimLab))) {i <- i + 1; prefix <- substr(dimLab, 1, i)}
      nodes$label <- paste0(prefix[nodes$fg], nodes$q)
    } else {
      nodes$label <- as.character(nodes$q)
    }
  } else {
    # prefix of the labels, one by type of nodes (as in the previous versions)
    nodes$label <- paste0(rep_len(unlist(vertex.label), nbFG)[nodes$fg], nodes$q)
  }

  #------------------ connectivities (edges of the mesoscopic graph)
  edges <- do.call(rbind, lapply(netList, function(net) {
    i <- net$fg[1]; j <- net$fg[2]
    m <- matrix(net$mean, nbBlocks[i], nbBlocks[j])
    e <- expand.grid(q = 1:nbBlocks[i], l = 1:nbBlocks[j])
    directed <- isTRUE(net$directed)
    if (i == j && !directed) {e <- e[e$q <= e$l, ]}
    data.frame(from = offset[i] + e$q, to = offset[j] + e$l,
               value = m[cbind(e$q, e$l)], directed = directed)
  }))
  # same colors as the graphon, on a scale common to all the networks
  allValues <- unlist(lapply(netList, function(net) net$mean))
  limits <- currentOptions$limits
  if (is.null(limits)) {limits <- if (all(allValues >= 0)) c(0, max(allValues)) else range(allValues)}
  alpha <- blockAlpha(edges$value, currentOptions, limits)
  edges$intensity <- (alpha - currentOptions$blockAlpha.min) / (currentOptions$blockAlpha.max - currentOptions$blockAlpha.min)
  edges$color <- vapply(alpha, function(a) blendWhite(currentOptions$block.color, a), character(1))
  edges$width <- currentOptions$edge.width * (0.2 + 2.8 * edges$intensity)
  edges <- edges[edges$value >= currentOptions$edge.threshold, ]
  edges <- edges[order(edges$intensity), ]

  #------------------ layout
  layout <- currentOptions$layout
  if (is.null(layout)) {layout <- if (nbFG == 1) 'free' else defaultLayout}
  if (is.character(layout)) {
    layout <- match.arg(layout, c('groups', 'free'))
    if (layout == 'groups' && nbFG > 1) {
      layout <- mesoLayoutFG(t(sapply(netList, function(net) net$fg)), nbBlocks)
    } else {
      # blocks attracted by their connectivity (the seed is set before the random layout)
      if (!is.null(currentOptions$seed)) {set.seed(currentOptions$seed)}
      inter <- edges[edges$from != edges$to, ]
      g <- igraph::make_empty_graph(n = Q, directed = FALSE)
      g <- igraph::add_edges(g, rbind(inter$from, inter$to), weight = 0.01 + inter$intensity)
      layout <- if (Q == 1) matrix(0, 1, 2) else igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
      # spread the blocks: the closest ones are at distance 1
      if (Q > 1) {layout <- layout / min(stats::dist(layout))}
    }
  }
  layout <- as.matrix(layout)
  nodes$x <- layout[, 1]
  nodes$y <- layout[, 2]

  # radius of the blocks: area proportional to the proportion, without overlap. The proportions
  # sum to 1 within each type of nodes: the largest block of each type has the maximal size
  dmin <- if (Q > 1) min(stats::dist(layout)) else 1
  piMax <- vapply(piList, max, numeric(1))[nodes$fg]
  nodes$r <- currentOptions$vertex.size[nodes$fg] * 0.4 * dmin * sqrt(nodes$pi / piMax)
  nodes$shape <- currentOptions$vertex.shape[nodes$fg]
  nodes$fill  <- currentOptions$vertex.color[nodes$fg]
  # simple SBM: the connectivity within the blocks can be given by the color of the nodes
  fillWithin <- nbFG == 1 && currentOptions$withinBlock == 'fill'
  if (fillWithin) {
    selfValue <- vapply(seq_len(Q), function(k) {
      v <- edges$value[edges$from == k & edges$to == k]
      if (length(v)) v else limits[1]}, numeric(1))
    nodes$fill <- vapply(blockAlpha(selfValue, currentOptions, limits),
                         function(a) blendWhite(currentOptions$block.color, a), character(1))
  }
  nodes$frame <- currentOptions$vertex.frame.color[nodes$fg]

  # free direction around each node: middle of the largest angular gap between the other blocks
  nodes$angle <- vapply(seq_len(Q), function(i) {
    others <- setdiff(seq_len(Q), i)
    if (length(others) == 0) return(base::pi / 2)
    angles <- sort(atan2(nodes$y[others] - nodes$y[i], nodes$x[others] - nodes$x[i]))
    gaps <- diff(c(angles, angles[1] + 2 * base::pi))
    angles[which.max(gaps)] + max(gaps) / 2
  }, numeric(1))

  # nodes as polygons (circles, squares, diamonds or triangles of the same area)
  polys <- do.call(rbind, lapply(seq_len(Q), function(k) {
    r <- nodes$r[k]
    xy <- switch(nodes$shape[k],
      square   = {h <- r * sqrt(base::pi) / 2; cbind(c(-h, h, h, -h), c(-h, -h, h, h))},
      diamond  = {h <- r * sqrt(base::pi / 2); cbind(c(0, h, 0, -h), c(-h, 0, h, 0))},
      triangle = {R <- r * sqrt(4 * base::pi / sqrt(27)); t <- base::pi / 2 + c(0, 2, 4) * base::pi / 3
                  cbind(R * cos(t), R * sin(t))},
      {t <- seq(0, 2 * base::pi, length.out = 61)[-61]; cbind(r * cos(t), r * sin(t))})
    data.frame(id = k, x = nodes$x[k] + xy[, 1], y = nodes$y[k] + xy[, 2])
  }))
  polys$fill  <- nodes$fill[polys$id]
  polys$frame <- nodes$frame[polys$id]

  #------------------ edges between blocks, from the border of the nodes
  inter <- edges[edges$from != edges$to, ]
  if (nrow(inter) > 0) {
    dx <- nodes$x[inter$to] - nodes$x[inter$from]
    dy <- nodes$y[inter$to] - nodes$y[inter$from]
    len <- sqrt(dx^2 + dy^2)
    # distance from the center to the border of the node, in the direction of the edge
    toBorder <- function(k) {
      ux <- abs(dx) / len; uy <- abs(dy) / len
      r <- nodes$r[k]
      ifelse(nodes$shape[k] == 'square', r * sqrt(base::pi) / 2 / pmax(ux, uy),
      ifelse(nodes$shape[k] == 'diamond', r * sqrt(base::pi / 2) / (ux + uy), r))
    }
    inter$x    <- nodes$x[inter$from] + dx / len * toBorder(inter$from)
    inter$y    <- nodes$y[inter$from] + dy / len * toBorder(inter$from)
    inter$xend <- nodes$x[inter$to]   - dx / len * toBorder(inter$to)
    inter$yend <- nodes$y[inter$to]   - dy / len * toBorder(inter$to)
    # undirected edges between aligned blocks would go through the blocks between them: curved
    inter$curved <- inter$directed | abs(dy) < 1e-8 * max(1, len)
  }
  # connectivity within a block: loop on the outer side of the node
  loops <- if (fillWithin) edges[0, ] else edges[edges$from == edges$to, ]
  loopPaths <- NULL
  if (nrow(loops) > 0) {
    loopPaths <- do.call(rbind, lapply(seq_len(nrow(loops)), function(k) {
      i <- loops$from[k]
      d <- c(cos(nodes$angle[i]), sin(nodes$angle[i]))
      rl <- 0.6 * nodes$r[i]
      cl <- c(nodes$x[i], nodes$y[i]) + d * (nodes$r[i] + 0.5 * rl)
      t <- seq(0, 2 * base::pi, length.out = 61)
      data.frame(id = k, x = cl[1] + rl * cos(t), y = cl[2] + rl * sin(t),
                 color = loops$color[k], width = loops$width[k])
    }))
  }

  #------------------ labels: inside the large blocks, outside the small ones (beyond their loop)
  # (minimal relative radius for a label of 2 characters, proportional to the number of characters)
  inside <- nodes$r >= pmin(1, currentOptions$vertex.label.inside * nchar(nodes$label) / 2) * max(nodes$r)
  hasLoop <- seq_len(Q) %in% loops$from
  ext <- nodes$r * ifelse(nodes$shape == 'circle', 1, 1.3) * ifelse(hasLoop, 1.9, 1) + 0.05 * dmin
  nodes$lx <- ifelse(inside, nodes$x, nodes$x + ext * cos(nodes$angle))
  nodes$ly <- ifelse(inside, nodes$y, nodes$y + ext * sin(nodes$angle))
  # text on the outer side of the anchor
  nodes$hjust <- ifelse(inside, 0.5, (1 - round(cos(nodes$angle), 10)) / 2)
  nodes$vjust <- ifelse(inside, 0.5, (1 - round(sin(nodes$angle), 10)) / 2)

  #------------------ plot
  g <- ggplot()
  if (!is.null(loopPaths)) {
    g <- g + geom_path(data = loopPaths, aes(x = .data$x, y = .data$y, group = .data$id,
                                              colour = .data$color, linewidth = .data$width),
                       lineend = 'round')
  }
  edgeAes <- aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                 colour = .data$color, linewidth = .data$width)
  if (nrow(inter) > 0) {
    if (any(!inter$curved)) {
      g <- g + geom_segment(data = inter[!inter$curved, ], edgeAes)
    }
    if (any(inter$curved & inter$directed)) {
      g <- g + geom_curve(data = inter[inter$curved & inter$directed, ], edgeAes,
                          curvature = currentOptions$edge.curved,
                          arrow = arrow(length = unit(2.5 * currentOptions$edge.arrow.size, "mm"), type = "closed"))
    }
    if (any(inter$curved & !inter$directed)) {
      g <- g + geom_curve(data = inter[inter$curved & !inter$directed, ], edgeAes,
                          curvature = max(abs(currentOptions$edge.curved), 0.3))
    }
  }
  g <- g +
    geom_polygon(data = polys, aes(x = .data$x, y = .data$y, group = .data$id,
                                   fill = .data$fill, colour = .data$frame)) +
    geom_text(data = nodes, aes(x = .data$lx, y = .data$ly, label = .data$label,
                                hjust = .data$hjust, vjust = .data$vjust),
              colour = currentOptions$vertex.label.color, size = currentOptions$vertex.label.size) +
    scale_colour_identity() + scale_fill_identity() + scale_linewidth_identity() +
    coord_fixed() +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(currentOptions$title)) {g <- g + ggtitle(currentOptions$title)}
  # coordinates of the blocks, to be modified and given back in plotOptions$layout
  attr(g, 'blockLayout') <- matrix(c(nodes$x, nodes$y), Q, 2, dimnames = list(nodes$label, c('x', 'y')))
  g
}

#-----------------------------------------------------------------
# Layout of the blocks of several functional groups (FG): the blocks of each FG are aligned on a row,
# spaced by 1. The FG with the most connections (the hub) is in the middle; the FG connected to it are
# placed on the rows above or below it (balancing the number of blocks), the FG further away beyond
# them, on the same side. FG on the same row are placed side by side.
# For a bipartite SBM: the row blocks above, the column blocks below.
#   arch     : matrix with 2 columns, the pairs of FG connected by a network
#   nbBlocks : number of blocks of each FG
#-----------------------------------------------------------------
mesoLayoutFG <- function(arch, nbBlocks){

  nbFG <- length(nbBlocks)
  A <- matrix(FALSE, nbFG, nbFG)
  A[arch] <- TRUE
  A <- A | t(A)
  diag(A) <- FALSE
  hub <- which.max(rowSums(A))

  # breadth first search from the hub: level = distance to the hub, side = +1 (above) or -1 (below)
  level <- rep(NA, nbFG); side <- rep(0, nbFG)
  level[hub] <- 0
  blocksBySide <- c(0, 0) # above, below
  queue <- hub
  repeat {
    if (length(queue) == 0) {
      # FG not connected to the previous ones: start again from the first one
      if (all(!is.na(level))) break
      f <- which(is.na(level))[1]
      s <- if (blocksBySide[1] < blocksBySide[2]) 1 else -1
      level[f] <- max(level, na.rm = TRUE) + 1; side[f] <- s
      blocksBySide[(3 - s) / 2] <- blocksBySide[(3 - s) / 2] + nbBlocks[f]
      queue <- f
    }
    f <- queue[1]; queue <- queue[-1]
    for (g in which(A[f, ] & is.na(level))) {
      level[g] <- level[f] + 1
      side[g] <- if (f == hub) {if (blocksBySide[1] < blocksBySide[2]) 1 else -1} else side[f]
      blocksBySide[(3 - side[g]) / 2] <- blocksBySide[(3 - side[g]) / 2] + nbBlocks[g]
      queue <- c(queue, g)
    }
  }
  row <- side * level

  # blocks of the FG of each row side by side, the FG being separated by a gap of 1.5
  x <- vector('list', nbFG)
  for (r in unique(row)) {
    fgs <- which(row == r)
    starts <- cumsum(c(0, nbBlocks[fgs] + 0.5))[seq_along(fgs)]
    width <- starts[length(fgs)] + nbBlocks[fgs[length(fgs)]] - 1
    for (k in seq_along(fgs)) {x[[fgs[k]]] <- starts[k] + seq_len(nbBlocks[fgs[k]]) - 1 - width / 2}
  }
  maxQ <- max(tapply(nbBlocks, row, sum))
  gap <- max(1.5, 0.6 * maxQ)
  cbind(x = unlist(x), y = rep(row * gap, nbBlocks))
}
