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
plotMatrix = function(Mat, rowFG, colFG, clustering = NULL){



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

  if (!is.null(clustering)){
    g <- g + geom_vline(data = sepCol,mapping = aes_string(xintercept = 'sep'),col = 'grey')
    g <- g + geom_hline(data = sepRow,mapping = aes_string(yintercept = 'sep'),col = 'grey')
  }
  #if (!is.null(fileNameSave)) { ggsave(fileNameSave, width = 20, height = 20, units = "cm") }else{g}
  g
}





