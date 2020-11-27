########################## SIMPLE
n1 = 25
n2 = 20
Mat <- matrix(sample(c(0,1),n1*n2,replace = TRUE),n1,n2);
colnames(Mat) <- sapply(1:ncol(Mat),function(i){
  li <- sample(2:7,1)
  paste(sample(letters,li,replace = TRUE),collapse = '')
})

rownames(Mat) <- sapply(1:nrow(Mat),function(i){
  li <- sample(2:7,1)
  paste(sample(letters,li,replace = TRUE),collapse = '')
}
)
dimLabels = list(row = 'truc',col='machin')

plotMatrix(Mat, dimLabels, clustering = NULL,plotOptions)

#------------- estimated

library(sbm)
resLBM <- estimateBipartiteSBM(Mat,model = 'bernoulli',dimLabels=dimLabels)
clustering = NULL
plotOptions = list(
  line.color  = 'red',
  legend = TRUE,
  legend.title = TRUE,
  rowNames = TRUE,
  colNames = TRUE,
  title=NULL)
plot(resLBM,plotOptions)


################# MULTI


U1 <- matrix(rbinom(1000,1,0.5),20,50)
colnames(U1) <- sapply(1:ncol(U1),function(i){li <- sample(2:7,1); paste(sample(letters,li,replace = TRUE),collapse = '')})
rownames(U1) <- sapply(1:nrow(U1),function(i){li <- sample(2:7,1);  paste(sample(letters,li,replace = TRUE),collapse = '')})
U2 <- matrix(rpois(20*30,8),30,20)
colnames(U2) <- sapply(1:ncol(U2),function(i){li <- sample(2:7,1); paste(sample(letters,li,replace = TRUE),collapse = '')})
rownames(U2) <- sapply(1:nrow(U2),function(i){li <- sample(2:7,1);  paste(sample(letters,li,replace = TRUE),collapse = '')})

listNet[[1]] <- defineSBM(U1,
                          model = 'bernoulli',
                          type  ='bipartite', directed = NA,
                          dimLabels = list(row="Questions",col="Students"))
listNet[[2]] <- defineSBM(U2,
                          model = 'poisson',
                          type  ='bipartite',directed = NA,
                          dimLabels = list(row="Competences",col="Questions"))
plotMyMultipartiteMatrix(listNet,plotOptions=list(legend = TRUE,compact = FALSE))
plotMyMultipartiteMatrix(listNet,plotOptions = list(legend = FALSE,normalized =  TRUE,nodeNames = TRUE))




