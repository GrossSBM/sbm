
#----------------------------------------------------------------------------------
plotMatrix = function(Mat, dimLabels, clustering = NULL){


  rowFG <- dimLabels$row
  colFG <- dimLabels$col


  n1 <- dim(Mat)[1]
  n2 <- dim(Mat)[2]
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
  }

  index_row = rep(1:dim(Mat)[1],each = dim(Mat)[2])
  index_col = rep(1:dim(Mat)[2],dim(Mat)[1])

  melted_Mat =  data.frame(n1 - index_row , index_col)
  link = rep(-10,dim(Mat)[2]*dim(Mat)[1])
  for (k in 1:(dim(Mat)[2] * dim(Mat)[1])) {link[k] = Mat[index_row[k],index_col[k]]}
  melted_Mat$link = link
  if (binary){melted_Mat$link <- as.factor(melted_Mat$link)}
  colnames(melted_Mat) <- c('index_row', 'index_col', 'link')

  g <- ggplot(data = melted_Mat, aes(y = index_row, x = index_col, fill = link))
  g <- g + geom_tile()
  if (!binary) {g <-  g +  scale_fill_gradient(low = "white", high = "black", limits = u,na.value = "transparent")}
  if (binary) {g <- g + scale_fill_manual(breaks = c("0", "1"),values = c("white", "black"),na.value = "transparent")}
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  g <- g + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  g <- g +  labs(x = colFG, y = rowFG) +  theme(aspect.ratio = n1/n2)

  if (!is.null(clustering)) {
    g <- g + geom_vline(data = sepCol,mapping = aes_string(xintercept = 'sep'),col = 'grey')
    g <- g + geom_hline(data = sepRow,mapping = aes_string(yintercept = 'sep'),col = 'grey')
  }
  g
  #if (!is.null(fileNameSave)) { ggsave(fileNameSave, width = 20, height = 20, units = "cm") }else{g}
}
#----------------------------------------------------------------------------------
plotMultipartiteMatrix = function(list_Mat, E, nbNodes, namesFG, normalized, clustering) {



  nbFG <- length(unique(c(E)))
  nbNet <- length(list_Mat)

  if (normalized){

    list_Mat = lapply(list_Mat,function(mat){
      val = c(mat)
      mat = (mat - min(val))/(max(val)-min(val))})
  }


  # meta matrix
  N <- sum(nbNodes)
  MetaMat <- matrix(NA, N, N)

  # network of matrices
  G <- matrix(0, nbFG, nbFG)
  G[E[, 1], E[, 2]] = 1
  uRow <- which(rowSums(G) > 0)
  uCol <- which(colSums(G) > 0)
  # places in metamatrix
  EndFG <- cumsum(nbNodes)
  BegFG <- c(0, cumsum(nbNodes)[-nbFG]) + 1
  #  reordering if clusterings

  #browser()
  if (!is.null(clustering)) {
    sepRow <- sepCol <- list()
    list_Mat_reorder <- list_Mat
    # reorder
    for (l in 1:nbNet) {
      clustering_row <- clustering[[E[l, 1]]]
      clustering_col <- clustering[[E[l, 2]]]
      oRow <- order(clustering_row)
      oCol <- order(clustering_col)
      list_Mat_reorder[[l]] <- list_Mat[[l]][oRow, oCol]
      uRowSep <- cumsum(table(clustering_row)) + 0.5
      uRowSep <- uRowSep[-length(uRowSep)]

      uColSep <- cumsum(table(clustering_col)) + 0.5
      uColSep <- uColSep[-length(uColSep)]
      sepCol[[l]] <- as.data.frame(BegFG[E[l, 2]] + uColSep)
      sepRow[[l]] <- as.data.frame(uRowSep)
      sepRow[[l]] = BegFG[E[l, 1]] + nbNodes[E[l, 1]] - sepRow[[l]]
      names(sepCol[[l]]) = names(sepRow[[l]]) = 'sep'
      sepRow[[l]]$FG <- namesFG[E[l, 1]]
      sepCol[[l]]$FG <- namesFG[E[l, 2]]
    }
    list_Mat <- list_Mat_reorder
    sepRow <- do.call(rbind,sepRow)
    sepCol <- do.call(rbind,sepCol)
  }



  #indFGCol <- indFGRow <- rep(1:nbFG, times = nbNodes)
  for (i in 1:nbFG) {
    for (j in 1:nbFG) {
      if (G[i, j] == 1) {
        place_row <-  BegFG[i]:EndFG[i]
        place_col <-  BegFG[j]:EndFG[j]
        kij <-
          which(rowSums(E == matrix(
            rep(c(i, j), nbNet),
            ncol = 2,
            nrow = nbNet,
            byrow = T
          )) == 2)

        MetaMat[place_row, place_col] <- list_Mat[[kij]]
      }
    }
  }
  FGCol <- FGRow <- rep(namesFG, times = nbNodes)
  uRow <- which(apply(MetaMat , 1, function(u){all(is.na(u))}))
  uCol <- which(apply(MetaMat , 2, function(u){all(is.na(u))}))

  # reduce matrix
  MetaMatSmall <- MetaMat
  FGColSmall <- FGCol
  FGRowSmall <- FGRow
  if (length(uRow) > 0 ){
    MetaMatSmall <- MetaMatSmall[-uRow,]
    FGRowSmall <- FGRowSmall[-uRow]
  }
  if (length(uCol) > 0 ){
    MetaMatSmall <- MetaMatSmall[,-uCol]
    FGColSmall <- FGColSmall[-uCol]
  }





  #indFGColSmall <- indFGCol[-uCol]
  #indFGRowSmall <- indFGRow[-uRow]
  binary = FALSE
  val <- sort(unique(c(MetaMatSmall)))
  if (setequal(val , c(0, 1))) {
    binary = TRUE
  }

  # changement en MELTED
  index_row = rep(1:dim(MetaMatSmall)[1], each = dim(MetaMatSmall)[2])
  index_col = rep(1:dim(MetaMatSmall)[2], dim(MetaMatSmall)[1])
  melted_Mat = data.frame(c(t(MetaMatSmall)))
  names(melted_Mat) <- 'link'
  melted_Mat$index_row <- dim(MetaMatSmall)[1] - index_row
  melted_Mat$index_col <- index_col
  melted_Mat$FG_row <- as.factor(rep(FGRowSmall, each = dim(MetaMatSmall)[2]))
  melted_Mat$FG_col <- as.factor(rep(FGColSmall,  dim(MetaMatSmall)[1]))
  if (binary) { melted_Mat$link <- as.factor(melted_Mat$link)}

  # plots

  g <- ggplot(melted_Mat, aes(y = index_row, x = index_col, fill = link))
  g <- g + geom_tile()
  name_legend <- ifelse(normalized,"Norm. Link","Link")
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  if (!binary) {
    g <-
      g +  scale_fill_gradient(
        low = "white",
        high = "black",
        limits = range(c(.data$link)),
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
  g <- g +  labs(x = '', y = '')
  g <- g + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  #browser()
  #if (!is.null(clustering)) {
  # g <- g + geom_vline(data = sepCol,aes(xintercept = 'sep'),col = 'orange')
  # g <- g + geom_hline(data = sepRow,aes(yintercept = 'sep'),col = 'orange')
  #}
  g <- g + facet_grid(FG_row ~ FG_col, scales = 'free', space = 'free')

  g
}

#-----------------------------------------------------------------
#' @importFrom graphics par
plotMeso <- function(thetaMean, pi,model,directed,bipartite,nbNodes,nodeLabels,plotOptions){


  #browser()
  currentOptions <- list(seed = NULL,
                         title = NULL,
                         layout = NULL,
                         vertex.color = 'salmon2',
                         vertex.frame.color = "black",#"white",                 # Node border color
                         vertex.shape = "circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                         vertex.size = 2,                               # Size of the node (default is 15)
                         vertex.size2 = NA,                              # The second size of the node (e.g. for a rectangle)
                         # === vertex label
                         vertex.label = NULL,                   # Character vector used to label the nodes
                         vertex.label.color =  "black",#"white",
                         vertex.label.font = 2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
                         vertex.label.cex = 0.9,                           # Font size (multiplication factor, device-dependent)
                         vertex.label.dist = 0,                          # Distance between the label and the vertex
                         vertex.label.degree = 0 ,                       # The position of the label in relation to the vertex (use pi)
                         # === Edge
                         edge.threshold = -Inf,
                         edge.color = "gray",#"white",                           # Edge color
                         edge.width = 10,                    # Edge width, defaults to 10
                         edge.arrow.size = 1,                            # Arrow size, defaults to 1
                         edge.arrow.width = 2,                           # Arrow width, defaults to 1
                         edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
                         edge.curved = 0.3)





  if (bipartite) {
    currentOptions$vertex.color = c('salmon2','darkolivegreen3')
    currentOptions$vertex.shape = c('circle','square')
    currentOptions$vertex.label = nodeLabels
    Qrow <- length(pi$row)
    Qcol <- length(pi$col)
    currentOptions$vertex.color = c(rep(currentOptions$vertex.color[1],Qrow),rep(currentOptions$vertex.color[2],Qcol))
    currentOptions$vertex.shape = c(rep(currentOptions$vertex.shape[1],Qrow),rep(currentOptions$vertex.shape[2],Qcol))
  }

  currentOptions[names(plotOptions)] <- plotOptions
  layout <- currentOptions$layout



  #------------------------------------------
  currentOptions[names(plotOptions)] <- plotOptions
  #-----------------------------------------------

  if (is.null(plotOptions$vertex.label)){
    vertex.label <- substr(currentOptions$vertex.label, 1, 1)
    i = 1
    while (length(unique(vertex.label)) < 2 & (i <= 2)) {i = i + 1; vertex.label <- substr( currentOptions$vertex.label, 1, i)}
  }else{
    vertex.label   = plotOptions$vertex.label
  }


  if (is.atomic(vertex.label)){vertex.label= as.list(vertex.label)}




  alpha <- thetaMean
  alpha[alpha < currentOptions$edge.threshold] <- 0

  if (model == "bernoulli"){alpha.norm  = alpha};
  if (model == "poisson"){alpha.norm = alpha/max(alpha)}
  if ( !(model %in% c("bernoulli","poisson"))){alpha.norm = (alpha - min(alpha))/(max(alpha) - min(alpha)) + 0.1}

  if (currentOptions$edge.threshold != -Inf) {
    cat(paste("Nota bene: threshold on connexions is",currentOptions$edge.threshold,sep = ' '))
  }

  if (bipartite){
    names(vertex.label) <- c('row','col')
    colnames(alpha.norm) <- paste(vertex.label$col,1:length(pi$col),sep='')
    rownames(alpha.norm) <- paste(vertex.label$row,1:length(pi$row),sep='')
    vlab <- c(rownames(alpha.norm),colnames(alpha.norm))
    g <- igraph::graph_from_incidence_matrix(alpha.norm, weighted = TRUE)
    u <- c(pi$row*nbNodes[1],pi$col*nbNodes[2])
    if (is.null(layout)){ layout <- igraph::layout_as_bipartite(g)}

  }else{
    vlab <- paste(vertex.label,1:length(pi),sep = "")
    u <- pi*nbNodes[1]
    mode <- ifelse(directed,'directed','undirected')
    g <- igraph::graph.adjacency(alpha.norm, mode = mode, weighted = TRUE)
    if (is.null(layout)){layout <- igraph::layout_with_fr(g)}

  }

  igraph::E(g)$width <- 1 + as.integer(igraph::E(g)$weight*10)

  set.seed(currentOptions$seed)

  old_par <- par(mar = rep(0.15,4))
  plot(g, layout = layout, ## see https://www.r-graph-gallery.com/248-igraph-plotting-parameters.html
       # === vertex
       vertex.color = currentOptions$vertex.color,
       vertex.frame.color = currentOptions$vertex.frame.color,
       vertex.shape = currentOptions$vertex.shape,
       vertex.size = currentOptions$vertex.size*u,
       vertex.size2 = currentOptions$vertex.size2,
       # === vertex label
       vertex.label = vlab,
       vertex.label.color = currentOptions$vertex.label.color,
       vertex.label.family = currentOptions$vertex.label.family,
       vertex.label.cex = currentOptions$vertex.label.cex,
       vertex.label.dist = currentOptions$vertex.label.dist,
       vertex.label.degree = currentOptions$vertex.label.degree,
       # === Edge
       edge.color = currentOptions$edge.color,
       edge.width = igraph::E(g)$weight*currentOptions$edge.width,
       edge.arrow.size = currentOptions$edge.arrow.size,
       edge.arrow.width = currentOptions$edge.arrow.width,
       edge.lty = currentOptions$edge.lty,
       edge.curved = currentOptions$edge.curved,
       main = currentOptions$title
  )
  par(old_par)
  list(g = g,layout = layout,plotOptions <- currentOptions)

}



#-----------------------------------------------------------------
plotMesoMultipartite <- function(E,list_thetaMean, list_pi,list_model,directed,bipartite,nbNodes,nodeLabels,plotOptions){
  #
  #
  #
  #   nbFG <- length(list_pi)
  #   nbNet <- nrow(E)
  #   nbBlocks <- sapply(list_pi,length)
  #   #---- les couleurs
  #
  #   color <- c('salmon2','darkolivegreen3','cadetblue2','darkgoldenrod2','lightsteelblue1','plum2','seagreen')
  #   if ( (nbFG > length(color)) & is.null(plotOptions$vertex.color)){stop('Too many FG. Define your own vertex.colors')}
  #   mycolor <- color[1:nbFG]
  #
  #
  #   currentOptions <- list(seed = NULL,
  #                          title = NULL,
  #                          layout = NULL,
  #                          vertex.color = mycolor,
  #                          vertex.frame.color = "black",                  # Node border color
  #                          vertex.shape = "circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
  #                          vertex.size = 2,                               # Size of the node (default is 15)
  #                          vertex.size2 = NA,                              # The second size of the node (e.g. for a rectangle)
  #                          # === vertex label
  #                          vertex.label = NULL,                   # Character vector used to label the nodes
  #                          vertex.label.color =  "black",
  #                          vertex.label.font = 2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
  #                          vertex.label.cex = 0.9,                           # Font size (multiplication factor, device-dependent)
  #                          vertex.label.dist = 0,                          # Distance between the label and the vertex
  #                          vertex.label.degree = 0 ,                       # The position of the label in relation to the vertex (use pi)
  #                          # === Edge
  #                          edge.threshold = -Inf,
  #                          edge.color = "gray",                           # Edge color
  #                          edge.width = 10,                    # Edge width, defaults to 10
  #                          edge.arrow.size = 1,                            # Arrow size, defaults to 1
  #                          edge.arrow.width = 2,                           # Arrow width, defaults to 1
  #                          edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
  #                          edge.curved = 0.3)
  #
  #
  #
  #   #### labels of the vertices depending on their FG
  #   vertex.Label   <- substr(nodeLabels,1,1)
  #   i = 1;
  #   while ( length(unique(vertex.Label))< nbFG & (i <= nbFG)){i = i + 1; vertex.Label <- substr(nodeLabels,1,i)}
  #
  #   #param <- resMBM$fittedModel[[whichModel]]$paramEstim
  #
  #   labelNode <- lapply(1:nbFG,function(q){paste(vertex.Label[q],1:nbBlocks[q],sep='')})
  #   colNode <- lapply(1:nbFG,function(q){rep(mycolor[q],nbBlocks[q])})
  #
  #   u <- c(pi$row*nbNodes[1],pi$col*nbNodes[2])
  #
  #   sizeNode <- lapply(1:nbFG,function(q){list_pi[[q]]*nbNodes[q]*currentOptions$vertex.size})
  #
  #   cumVK <-  c(0,cumsum(nbBlocks))
  #   codeNode <- lapply(2:(nbFG + 1),function(q){seq(cumVK[q - 1] + 1,cumVK[q],1)})
  #
  #   list_edges <- lapply(1:nbNet, function(i) {
  #     qRow <- E[i, 1]
  #     qCol <- E[i, 2]
  #     list_theta_i <- list_theta[[i]]
  #     c1 <- rep(codeNode[[qRow]], times = nbBlocks[qCol])
  #     c2 <- rep(codeNode[[qCol]], each = nbBlocks[qRow])
  #     edges_i <- cbind(c1, c2, c(list_theta_i))
  #     edges_i <- as.data.frame(edges_i)
  #     edges_i$type <- rep(dataR6$typeInter[i], length(c1))
  #     return(edges_i)
  #   })
  # allEdges <- do.call("rbind", list_edges)
  #
  #
  # allEdges$arrow_mode <- rep(0,nrow(allEdges))  # directed or nont directed
  # allEdges$arrow_mode[allEdges$type == "diradj"] = 2
  #
  # # allEdges <- rbind(allEdges,c(2,2,10))
  # # allEdges$type[nrow(allEdges)] = "diradj"
  # w <- which(allEdges[,3] > thres)
  # edges <- allEdges[w,c(1,2)]
  # curved <- 0 * (allEdges[w,4] == "diradj")
  #

}
