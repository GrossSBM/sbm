as_indicator <- function(clustering) {
  K <- length(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

as_clustering <- function(indicator) {
  if (is.null(indicator)) {
    cl <- numeric(0)
  } else {
    cl <- apply(indicator, 1, which.max)
  }
  cl
}

## Some utils function for math
.xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))
.logistic <- function(x) {1/(1 + exp(-x))}
.logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

.na2zero <- function(x) {
  x[is.na(x)] <- 0
  x
}

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}


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
plotMultipartiteMatrix = function(list_Mat, E, nbNodes, namesFG,normalizing,clustering) {


  nbFG <- length(unique(c(E)))
  nbNet <- length(list_Mat)

  if (normalizing){

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

  g <- ggplot(data = melted_Mat, aes(y = index_row, x = index_col, fill = link))
  g <- g + geom_tile()
  name_legend <- ifelse(normalizing,"Norm. Link","Link")
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  if (!binary) {
    u <- range(c(melted_Mat$link))
    g <-
      g +  scale_fill_gradient(
        low = "white",
        high = "black",
        limits = u,
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
plotMeso <- function(thetaMean, pi,directed,bipartite,nbNodes,nodeLabels,plotOptions){

  currentOptions <- list(seed = NULL,
                         title = NULL,
                         layout = NULL,
                         vertex.color = 'salmon2',
                         vertex.frame.color = "black",#"white",                 # Node border color
                         vertex.shape = "circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                         vertex.size = 2,                               # Size of the node (default is 15)
                         vertex.size2 = NA,                              # The second size of the node (e.g. for a rectangle)
                         # === vertex label
                         vertex.label = nodeLabels,                   # Character vector used to label the nodes
                         vertex.label.color =  "black",#"white",
                         vertex.label.font = 2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
                         vertex.label.cex = 0.9,                           # Font size (multiplication factor, device-dependent)
                         vertex.label.dist = 0,                          # Distance between the label and the vertex
                         vertex.label.degree = 0 ,                       # The position of the label in relation to the vertex (use pi)
                         # === Edge
                         edge.threshold = -Inf,
                         edge.color = "gray",#"white",                           # Edge color
                         edge.width = 10,                    # Edge width, defaults to 1
                         edge.arrow.size = 1,                            # Arrow size, defaults to 1
                         edge.arrow.width = 2,                           # Arrow width, defaults to 1
                         edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
                         edge.curved = 0.3)



  if (bipartite) {
    currentOptions$vertex.color = c('salmon2','darkolivegreen3')
    currentOptions$vertex.shape = c('circle','square')
    Qrow <- length(pi$row)
    Qcol <- length(pi$col)
    currentOptions$vertex.color = c(rep(currentOptions$vertex.color[1],Qrow),rep(currentOptions$vertex.color[2],Qcol))
    currentOptions$vertex.shape = c(rep(currentOptions$vertex.shape[1],Qrow),rep(currentOptions$vertex.shape[2],Qcol))
    currentOptions$layout = layout_as_bipartite
  }

  currentOptions[names(plotOptions)] <- plotOptions

  alpha <- thetaMean
  alpha[alpha < currentOptions$edge.threshold] <- 0
  alpha.norm = (alpha - min(alpha))/(max(alpha) - min(alpha))
  if (currentOptions$threshold != Inf) {
    cat(paste("Nota bene: threshold on connexions is",currentOptions$threshold,sep = ' '))
  }

  if (bipartite){
    colnames(alphanorm) <- paste(currentOptions$vertex.label$col,1:length(pi$col),sep='')
    rownames(alphanorm) <- paste(currentOptions$vertex.label$row,1:length(pi$row),sep='')
    vlab <- c(rownames(alpha.norm),colnames(alpha.norm))
    g <- graph_from_incidence_matrix(alpha.norm, weighted = TRUE)
     u <- c(pi$row*nbNodes[1],pi$col*nbNodes[2])
  }else{
    vlab <- paste(currentOptions$vertex.label,1:length(pi),sep = "")
    u <- pi*nbNodes
    mode <- ifelse(directed,'directed','undirected')
    g <- graph.adjacency(alpha.norm, mode = mode, weighted = TRUE)
  }

  E(g)$width <- 1 + as.integer(E(g)$weight*10)

  set.seed(currentOptions$seed)


  par(mar = rep(0.15,4))
  plot(g, layout = currentOptions$layout, ## see https://www.r-graph-gallery.com/248-igraph-plotting-parameters.html
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
       edge.width = E(g)$weight*currentOptions$edge.width,
       edge.arrow.size = currentOptions$edge.arrow.size,
       edge.arrow.width = currentOptions$edge.arrow.width,
       edge.lty = currentOptions$edge.lty,
       edge.curved = currentOptions$edge.curved,
       main = currentOptions$title

  )


}





