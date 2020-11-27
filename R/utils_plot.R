myRepeat <- function(v,Qrow,Qcol){c(rep(v[1],Qrow),rep(v[2],Qcol))}
#' @importFrom rlang .data
#' @importFrom utils head

#----------------------------------------------------------------------------------
plotMatrix = function(Mat, dimLabels, clustering = NULL,plotOptions = list()){




  currentOptions = list(line.color  = 'red',legend = FALSE,rowNames = FALSE, colNames = FALSE,title=NULL)
  currentOptions$legend.title = FALSE
  currentOptions$legend.position ='bottom'

  currentOptions[names(plotOptions)] = plotOptions

  n1 <- dim(Mat)[1]
  n2 <- dim(Mat)[2]

  if (is.null(rownames(Mat))){rownames(Mat) = as.factor(1:n1)}
  if (is.null(colnames(Mat))){colnames(Mat) = as.factor(1:n2)}


  rowFG <- dimLabels$row
  colFG <- dimLabels$col

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
  }
  if(is.null(currentOptions$line.width)){currentOptions$line.width = mean(dim(Mat))/300}



  FGCol <- rep(dimLabels$col, times = n2)
  FGRow <- rep(dimLabels$row, times = n1)

  melted_Mat = reshape::melt(t(Mat))

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
  g <- g + geom_tile()
  if (!binary) {g <-  g +  scale_fill_gradient(low = "white", high = "black", limits = u,na.value = "transparent")}
  if (binary) {g <- g + scale_fill_manual(breaks = c("0", "1"),values = c("white", "black"),na.value = "transparent")}
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  if (!currentOptions$colNames){
    g <- g + theme(axis.text.x = element_blank())
  }else{
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 0))
  }
  if (!currentOptions$rowNames){
    g <- g + theme(axis.text.y = element_blank())
  }

  g <- g +  labs(x = '', y = '') +  theme(aspect.ratio = n1/n2, axis.ticks = element_blank(), panel.background = element_rect(fill = "white"))
  if (!is.null(dimLabels$row) & !is.null(dimLabels$col)){
    g <- g+ facet_grid(FG_row ~ FG_col,scales = 'free', space = 'free')}
  if (!currentOptions$legend){g <- g +theme(legend.position = 'none')}else{
    g <- g +theme(legend.position = currentOptions$legend.position)
    if(!currentOptions$legend.title){g <- g+ theme(legend.title = element_blank())}}




  if (!is.null(clustering)) {
    g <- g + geom_vline(data = sepCol,mapping = aes_string(xintercept = 'sep'),size = currentOptions$line.width, col=currentOptions$line.color)
    g <- g + geom_hline(data = sepRow + 1,mapping = aes_string(yintercept = 'sep'),size = currentOptions$line.width, col=currentOptions$line.color)
  }
  if (!is.null(currentOptions$title)){g <- g + ggtitle(currentOptions$title) }
  g
  #if (!is.null(fileNameSave)) { ggsave(fileNameSave, width = 20, height = 20, units = "cm") }else{g}
}
##################################################################################
#----------------------------------------------------------------------------------
##################################################################################

plotMultipartiteMatrix = function(listMat, E, nbNodes, namesFG, distrib, clustering,plotOptions) {


  #----------------------------------------

  currentOptions = list(line.color  = 'red',legend = FALSE,compact = TRUE, normalized = FALSE,title=NULL)
  currentOptions$legend.title = FALSE
  currentOptions$legend.position ='bottom'
  currentOptions$nodeNames = FALSE

  currentOptions[names(plotOptions)] = plotOptions

  normalized  <- currentOptions$normalized
  reordered <- !is.null(clustering)
  #-------------------------------------------
  list_Mat <- listMat
  nbFG <- length(unique(c(E)))
  nbNet <- length(list_Mat)
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
  binary = ifelse(all(distrib == 'bernoulli'),TRUE,FALSE)

  ############## optimisation plot
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

        if (length(u) > 0) {
          for (i in u) {
            Ex[i, 1:2] = E[i, 2:1]
          }
        }
        Ex
      }
    )
    GColl <- sapply(
      1:(2 ^ nbNet / 2),
      FUN = function(i) {
        G <- matrix(0, nbFG, nbFG)
        for (j in 1:nbNet) {
          G[EColl[[i]][j, 1], EColl[[i]][j, 2]] = 1
        }
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
  if(is.null(currentOptions$line.width)){currentOptions$line.width = mean(dim(MetaMat))/300}
  FGCol <- rep(namesFG, times = nbNodes*GCol)
  FGRow <- rep(namesFG, times = nbNodes*GRow)

  ############# meltedMat




  melted_Mat = reshape::melt(t(MetaMat))
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


  g <- g + facet_grid(FG_row~ FG_col, scales = 'free', space = 'free')

  ########## separators
  if (reordered){
    separate <- unlist(lapply(1:nbFG, function(l){
      sepColl <- cumsum(table(clustering[[l]])) + 0.5
      sepColl <- sepColl[-length(sepColl)]
      sepColl}))
    nbSep <- length(separate)
    separCol <- data.frame(sepCol = (separate + rep(c(0,cumsum(GCol*nbNodes)[-nbFG]),nbBlocks-1)),
                           FG_col = rep(namesFG,nbBlocks-1),
                           FG_col_index = rep(1:nbFG,nbBlocks-1))
    separCol <- do.call("rbind", replicate( nbFG,separCol ,simplify = FALSE))
    separCol$FG_row <- rep(namesFG,each = nbSep)
    separCol$FG_row_index <- rep(1:nbFG,each = nbSep)
    testCol <- vapply(1:(nbSep*nbFG),function(i){1*(G[separCol$FG_row_index[i],separCol$FG_col_index[i]]==1)},1)
    separCol <- separCol[testCol==1,]

    separRow <- data.frame(sepRow = 1 + sum(nbNodes*GRow) -(separate + rep(c(0,cumsum(GRow*nbNodes)[-nbFG]),nbBlocks-1)),
                           FG_row = rep(namesFG,nbBlocks-1),
                           FG_row_index = rep(1:nbFG,nbBlocks-1))
    separRow <- do.call("rbind", replicate( nbFG,separRow ,simplify = FALSE))
    separRow$FG_col <- rep(namesFG,each = nbSep)
    separRow$FG_col_index <- rep(1:nbFG,each = nbSep)
    testRow <- vapply(1:(nbSep*nbFG),function(i){1*(G[separRow$FG_row_index[i],separRow$FG_col_index[i]]==1)},1)
    separRow <- separRow[testRow==1,]

    g <- g + geom_vline(data= separCol, aes(xintercept =  .data$sepCol),size= currentOptions$line.width, col=currentOptions$line.color)
    g <- g + geom_hline(data= separRow, aes(yintercept = .data$sepRow),size= currentOptions$line.width, col=currentOptions$line.color)
  }
  if (!is.null(currentOptions$title)){g <- g + ggtitle(currentOptions$title) }
  g
}

#-----------------------------------------------------------------
#' @importFrom graphics par
#' @importFrom graphics plot
plotMeso <- function(thetaMean, pi,model,directed,bipartite,nbNodes,nodeLabels,plotOptions){



  currentOptions <- list(seed = NULL,
                         title = NULL,
                         layout = NULL,
                         vertex.color = 'salmon2',
                         vertex.frame.color = "black",#"white",                 # Node border color
                         vertex.shape = "circle",                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
                         vertex.size = 1,                               # Size of the node (default is 15)
                         vertex.size2 = NA,                              # The second size of the node (e.g. for a rectangle)
                         # === vertex label
                         vertex.label.name = NULL,                   # Character vector used to label the nodes
                         vertex.label.color =  "black",#"white",
                         vertex.label.font = 2,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
                         vertex.label.cex = 0.9,                           # Font size (multiplication factor, device-dependent)
                         vertex.label.dist = 0,                          # Distance between the label and the vertex
                         vertex.label.degree = 0 ,                       # The position of the label in relation to the vertex (use pi)
                         # === Edge
                         edge.threshold = -Inf,
                         edge.color = "gray",#"white",                           # Edge color
                         edge.width = 1,                    # Edge width, defaults to 10
                         edge.arrow.size = 1,                            # Arrow size, defaults to 1
                         edge.arrow.width = 2,                           # Arrow width, defaults to 1
                         edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
                         edge.curved = 0.3
                         )


  if (bipartite) {
    currentOptions$vertex.color = c('salmon2','darkolivegreen3')
    currentOptions$vertex.shape = c('circle','square')
    currentOptions$vertex.label.name = nodeLabels
    currentOptions$vertex.size  = c(1,1)
    currentOptions$vertex.label.cex =c(0.9,0.9)
    currentOptions$vertex.label.color =  c("black","black")
    currentOptions$vertex.label.dist =  c(0,0)
    currentOptions$vertex.frame.color = c('black','black')
    currentOptions$vertex.size2  = c('NA','NA')
    currentOptions$vertex.label.font = c(2,2)                           # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
    currentOptions$vertex.label.degree = c(0,0)

    if ((length(plotOptions$vertex.size) == 1)){plotOptions$vertex.size = rep(plotOptions$vertex.size,2)}
    if ((length(plotOptions$vertex.label.cex) == 1)){plotOptions$vertex.label.cex = rep(plotOptions$vertex.label.cex,2)}
    if ((length(plotOptions$vertex.label.font) == 1)){plotOptions$vertex.label.font = rep(plotOptions$vertex.label.font,2)}
    if ((length(plotOptions$vertex.label.degree) == 1)){plotOptions$vertex.label.degree = rep(plotOptions$vertex.label.degree,2)}



  }
  #---------------------------------------"

  currentOptions[names(plotOptions)] <- plotOptions
  #---------------------------------------"

  layout <- currentOptions$layout

  if (is.null(plotOptions$vertex.label.name)){
    vertex.label <- substr(currentOptions$vertex.label.name, 1, 1)
    i = 1
    while (length(unique(vertex.label)) < 2 & (i <= 2)) {i = i + 1; vertex.label <- substr(currentOptions$vertex.label.name, 1, i)}
  }else{
    vertex.label   = plotOptions$vertex.label.name
  }


  myOptions <- currentOptions
  if (bipartite){
  Qrow <- length(pi$row)
  Qcol <- length(pi$col)
  w.vertex.options <- stringr::str_detect(names(currentOptions), "vertex.")
  myOptions <- lapply(1:length(currentOptions),
                      function(p){if (w.vertex.options[p]){myRepeat(currentOptions[[p]],Qrow,Qcol)}else{currentOptions[[p]]}})
  names(myOptions) <- names(currentOptions)
  }




  if (is.atomic(vertex.label)){vertex.label= as.list(vertex.label)}




  alpha <- thetaMean
  alpha[alpha < currentOptions$edge.threshold] <- 0

  if (model == "bernoulli"){alpha.norm  = alpha};
  if (model == "poisson"){alpha.norm = alpha/max(alpha)}
  if ( !(model %in% c("bernoulli","poisson"))){alpha.norm = (alpha - min(alpha))/(max(alpha) - min(alpha)) + 0.1}

  if (currentOptions$edge.threshold != -Inf) {
    cat(paste("Nota bene: threshold on connections is",currentOptions$edge.threshold,sep = ' '))
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
  plot(g, layout = layout,
       # === vertex
       vertex.color = myOptions$vertex.color,
       vertex.frame.color = myOptions$vertex.frame.color,
       vertex.shape = myOptions$vertex.shape,
       vertex.size = myOptions$vertex.size*u,
       vertex.size2 = myOptions$vertex.size2,
       # === vertex label
       vertex.label = vlab,
       vertex.label.color = myOptions$vertex.label.color,
       vertex.label.family = myOptions$vertex.label.family,
       vertex.label.cex = myOptions$vertex.label.cex,
       vertex.label.dist = myOptions$vertex.label.dist[1],
       vertex.label.degree = myOptions$vertex.label.degree,
       # === Edge
       edge.color = myOptions$edge.color,
       edge.width = 10*igraph::E(g)$weight*myOptions$edge.width,
       edge.arrow.size = myOptions$edge.arrow.size,
       edge.arrow.width = myOptions$edge.arrow.width,
       edge.lty = myOptions$edge.lty,
       edge.curved = myOptions$edge.curved,
       main = myOptions$title
  )
  par(old_par)
  list(g = g,layout = layout,plotOptions = currentOptions)

}


##################################################################################
#----------------------------------------------------------------------------------
##################################################################################
#-----------------------------------------------------------------
#' @importFrom graphics par
#' @importFrom graphics plot
plotMesoMultipartite <- function(E,theta, list_pi,v_distrib,directed,nbNodes,nodeLabels,plotOptions){


  directed[is.na(directed)] <- 'FALSE'
  nbFG <- length(list_pi)
  nbNet <- nrow(E)
  nbBlocks <- sapply(list_pi,length)

  list_thetaMean <- lapply(theta,function(u){u$mean})

  #------------------------------------------list of colors and forms  ---------
  colors <- c('salmon2','darkolivegreen3','cadetblue2','darkgoldenrod2','lightsteelblue1','plum2','seagreen')
  if ( (nbFG > length(colors)) & is.null(plotOptions$vertex.color)){stop('Too many FG. Define your own vertex.color')}
  forms <- c('circle', 'square', 'csquare', 'rectangle' , 'crectangle', 'vrectangle', 'pie', 'raster', 'sphere')
  if ( (nbFG > length(colors)) & is.null(plotOptions$vertexshape)){stop('Too many FG. Define your own vertex.shape')}

  #----------------------------------------------------------------------------------

  currentOptions <- list(seed = NULL,
                         title = NULL,
                         layout = NULL,
                         vertex.color = colors[1:nbFG],
                         vertex.frame.color = rep("black",nbFG),#"white",                 # Node border color
                         vertex.shape = forms[1:nbFG],                        #
                         vertex.size = rep(1,nbFG),                               # Size of the node (default is 15)
                         vertex.size2 = rep(NA,nbFG),                              # The second size of the node (e.g. for a rectangle)
                         # === vertex label
                         vertex.label.name = nodeLabels,                   # Character vector used to label the nodes
                         vertex.label.color =  rep("black",nbFG),#"white",
                         vertex.label.font = rep(2,nbFG),                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
                         vertex.label.cex = rep(0.9,nbFG),                           # Font size (multiplication factor, device-dependent)
                         vertex.label.dist = rep(0,nbFG),                          # Distance between the label and the vertex
                         vertex.label.degree = rep(0,nbFG) ,                       # The position of the label in relation to the vertex (use pi)
                         # === Edge
                         edge.threshold = -Inf,
                         edge.color = "gray",#"white",                           # Edge color
                         edge.width = 1,                    # Edge width, defaults to 10
                         edge.arrow.size = 1,                            # Arrow size, defaults to 1
                         edge.arrow.width = 1,                           # Arrow width, defaults to 1
                         edge.lty = "solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
                         edge.curved = 0.3)

  currentOptions[names(plotOptions)] <- plotOptions
  if (currentOptions$edge.threshold != -Inf) {
    cat(paste("Nota bene: threshold on connections is",currentOptions$edge.threshold,sep = ' '))
  }
  #---------------------------------------"
  if (is.null(plotOptions$vertex.label.name)){
    vertex.label <- substr(currentOptions$vertex.label.name, 1, 1)
    i = 1
    while (length(unique(vertex.label)) < 2 & (i <= 2)) {i = i + 1; vertex.label <- substr(currentOptions$vertex.label.name, 1, i)}
  }else{
    vertex.label   = plotOptions$vertex.label.name
  }
  #if (is.atomic(vertex.label)){vertex.label= as.list(vertex.label)}

  layout <- currentOptions$layout
  myOptions <- currentOptions
  w.vertex.options <- stringr::str_detect(names(currentOptions), "vertex.")
  myOptions <- lapply(1:length(currentOptions),
    function(p){if (w.vertex.options[p]){rep(currentOptions[[p]],nbFG)}else{currentOptions[[p]]}})
  names(myOptions) <- names(currentOptions)
  labelNode <- unlist(lapply(1:nbFG,function(q){paste(vertex.label[q],1:nbBlocks[q],sep='')}))

  #--------------------------- meatagraph g

  N <- sum(nbBlocks);
  DIR <- alpha <- matrix(0, N,N);
  alpha.norm <- alpha
  G <- matrix(0, nbFG, nbFG);
  G[E[, 1], E[, 2]] = 1
  EndFG <- cumsum(nbBlocks);  BegFG <- c(0, cumsum(nbBlocks)[-nbFG]) + 1
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
        alpha_rc <- list_thetaMean[[kij]]
        alpha_rc[alpha_rc < currentOptions$edge.threshold] <- 0
        alpha[place_row, place_col] <- alpha_rc
        model_kij <- v_distrib[kij]

        if (model_kij == "bernoulli"){alpha.norm_rc  = alpha_rc};
        if (model_kij == "poisson"){alpha.norm_rc = alpha_rc/max(alpha_rc)}
        if ( !(model_kij %in% c("bernoulli","poisson"))){alpha.norm_rc = (alpha_rc - min(alpha_rc))/(max(alpha_rc) - min(alpha_rc)) + 0.1}
        alpha.norm[place_row, place_col] <- alpha.norm_rc
        DIR[place_row,place_col] = directed[kij]
      }
    }
  }



  colnames(alpha.norm) <- rownames(alpha.norm) <- labelNode
  g <- igraph::graph.adjacency(alpha.norm, mode = 'directed', weighted = TRUE,diag = TRUE)
  if (is.null(layout)){layout <- igraph::layout_with_fr(g)}
  u <- unlist(lapply(1:nbFG, function(p){list_pi[[p]]*100}))

  igraph::E(g)$width <- 1 + as.integer(igraph::E(g)$weight*10)


  DIR <-c(c(t(DIR))[c(t(alpha.norm)) != 0] != 'FALSE')


  currentOptions$seed <- .Random.seed
  set.seed(currentOptions$seed)
  old_par <- par(mar = rep(0.15,4))
  plot(g, layout = layout,
     # === vertex
     vertex.color = myOptions$vertex.color,
     vertex.frame.color = myOptions$vertex.frame.color,
     vertex.shape = myOptions$vertex.shape,
     vertex.size = myOptions$vertex.size*u,
     vertex.size2 = myOptions$vertex.size2,
     # === vertex label
     vertex.label = labelNode,
     vertex.label.color = myOptions$vertex.label.color,
     vertex.label.family = myOptions$vertex.label.family,
     vertex.label.cex = myOptions$vertex.label.cex,
     vertex.label.dist = myOptions$vertex.label.dist[1],
     vertex.label.degree = myOptions$vertex.label.degree,
     # === Edge
     edge.color = myOptions$edge.color,
     edge.width = 10*igraph::E(g)$weight*myOptions$edge.width,
     edge.arrow.size = max(DIR)*myOptions$edge.arrow.size,
     edge.arrow.width = max(DIR)*myOptions$edge.arrow.width,
     edge.lty = myOptions$edge.lty,
     edge.curved = myOptions$edge.curved,
     main = myOptions$title
     )
    par(old_par)
    list(g = g,layout = layout,plotOptions = currentOptions)



}
