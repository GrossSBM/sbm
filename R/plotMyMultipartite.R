#' Plot the matrices corresonding to a Multipartite Network
#'
#' @param list_Net  : a list of objects representing the multipartite network (see)
#' @return a ggplot object corresponding to the plot
#' @export
#'
#' @examples
#' data(MPEcoNetwork, package = "GREMLIN")
#' PlantFlovis = sbm::defineNetwork(MPEcoNetwork$Inc_plant_flovis, type = "bipartite",model = 'bernoulli',dimLabels = list(row="Plants",col="Flovis"))
#' PlantAnt = sbm::defineNetwork(MPEcoNetwork$Inc_plant_ant,type = "bipartite",model = 'bernoulli',dimLabels =list(row = "Plants", col = "Ants"))
#' PlantBird = sbm::defineNetwork(MPEcoNetwork$Inc_plant_bird,type = "bipartite",model = 'bernoulli',dimLabels =list(row = "Plants",col = "Birds"))
#' plotMyMultipartiteMatrix(list(PlantFlovis,PlantAnt,PlantBird)
#'
#'

plotMyMultipartiteMatrix = function(list_Net){



  mydataR6 <- formattingData(list_Net)

  E <- mydataR6$E
  Q <- mydataR6$Q
  v_NQ <- mydataR6$v_NQ
  G <- matrix(0,Q,Q)
  G[E[,1],E[,2] ] = 1


  uRow <- which(rowSums(G)>0)
  uCol <- which(colSums(G)>0)

  cardE <- length(list_Net)

  EndBlocks <- cumsum(v_NQ)
  BegBlocks <- c(0,cumsum(v_NQ)[-Q]) +1

  N <- sum(v_NQ)
  MetaMat <- matrix(NA,N,N)
  dim(MetaMat)

  FGCol <- FGRow <- rep(mydataR6$namesFG, times = v_NQ)
  indFGCol <- indFGRow <- rep(1:Q, times = v_NQ)
  for (i in 1:Q){
    for (j in 1:Q){
      if (G[i,j] == 1){
      place_row <-  BegBlocks[i]: EndBlocks[i]
      place_col <-  BegBlocks[j]: EndBlocks[j]
      kij <- which(rowSums(E == matrix(rep(c(i,j),cardE),ncol=2,nrow= cardE,byrow = T)) == 2)

      MetaMat[place_row, place_col] <- mydataR6$mats[[kij]]
      }
    }
  }
  uRow <- which(apply(MetaMat , 1,function(u){all(is.na(u) )}))
  uCol <- which(apply(MetaMat , 2,function(u){all(is.na(u))}))
  MetaMatSmall <- MetaMat[-uRow, -uCol]
  FGColSmall <- FGCol[-uCol]
  FGRowSmall <- FGRow[-uRow]
  indFGColSmall <- indFGCol[-uCol]
  indFGRowSmall <- indFGRow[-uRow]
  binary = FALSE
  val <- sort(unique(c(MetaMatSmall)))
  if (setequal(val ,c(0,1))) {binary = TRUE}

  # changement en MELTED
  index_row = rep(1:dim(MetaMatSmall)[1],each = dim(MetaMatSmall)[2])
  index_col = rep(1:dim(MetaMatSmall)[2],dim(MetaMatSmall)[1])
  melted_Mat = data.frame(c(t(MetaMatSmall)))
  names(melted_Mat) <- 'link'
  melted_Mat$index_row <- dim(MetaMatSmall)[1] - index_row
  #melted_Mat$index_row <- index_row
  melted_Mat$index_col <- index_col
  melted_Mat$FG_row <- as.factor(rep(FGRowSmall,each = dim(MetaMatSmall)[2]))
  melted_Mat$FG_col <- as.factor(rep(FGColSmall,  dim(MetaMatSmall)[1]))
  if (binary){melted_Mat$link <- as.factor(melted_Mat$link)}

  g <- ggplot(data = melted_Mat, aes(y = index_row, x = index_col, fill = link))
  g <- g + geom_tile()
  g <- g  +  scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE)
  if (!binary) {
    u <- range(c(melted_Mat$link))
    g <-  g +  scale_fill_gradient(low = "white", high = "black", limits = u,na.value = "transparent")
  }else{
    g <- g + scale_fill_manual(breaks = c("0", "1"),values = c("white", "black"),na.value = "transparent")
  }
  g <- g +  labs(x = '', y = '')
  g <- g + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  g <- g + facet_grid(FG_row ~ FG_col,scales='free', space='free')
  g





}

