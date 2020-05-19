############################################################
############ Basic functions to use in ScoreSBM ###############


#-------------------------------------------------------------------------------
# Computes the number of (symmetric) pairs
#-------------------------------------------------------------------------------
n2nbPairs <- function(n, symmetric, diag=FALSE){
  if(symmetric){
    N <- ifelse(diag, n*(n+1)/2, n*(n-1)/2)
  }else{
    N <- ifelse(diag, N <- n^2, n*(n-1))
  }
  return(N)
}

#-------------------------------------------------------------------------------
# Computes the individual from the number of (symmetric) pairs
#-------------------------------------------------------------------------------
nbPairs2n <- function(N, symmetric, diag=FALSE){
  if(symmetric){
    n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
  }else{
    n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
  }
  return(n)
}

#-------------------------------------------------------------------------------
# Builds a (symmetric) matrix from a vector
#-------------------------------------------------------------------------------
vect2Mat <- function(V, symmetric, diag=FALSE)
{

  N <- length(V); n <- nbPairs2n(N, symmetric=symmetric, diag=diag);
  M <- matrix(0,n,n)

  if (symmetric == TRUE){
    # n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
    # M <- matrix(0,n,n)
    M[lower.tri(M, diag = diag)] <- V
    B <- t(M)
    diag(B) <- 0
    M <- B + M
  }else{
    # n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
    # M <- matrix(0,n,n)
    where <- (lower.tri(M, diag = diag) +  upper.tri(M, diag = diag)) > 0
    M[where] <- V
  }
  return(M)
}

#-------------------------------------------------------------------------------
# Builds a vector from a (symmetric) matrix
#-------------------------------------------------------------------------------
mat2Vect <- function(M, symmetric, diag=FALSE)
{
  if (symmetric) {
    where <- lower.tri(M, diag = diag)
  } else {
    where <- (lower.tri(M, diag=diag) +  upper.tri(M, diag = diag)) > 0
  }
  V <- M[where]
  return(V)
}

#-------------------------------------------------------------------------------
#- Gest the indices corresponding to a list of (symmetric) pairs
#-------------------------------------------------------------------------------
indices <- function(n, symmetric, diag=FALSE)
{
  N <- n2nbPairs(n, symmetric=symmetric, diag=diag)
  S <- vect2Mat(1:N, symmetric=symmetric, diag=diag)
  if (symmetric) {S[upper.tri(S)] = 0}
  res <- which(S!=0, arr.ind=TRUE)
  return(res)
}



#-------------------------------------------------------------------------------
#---- Transorm the list of the matrices into a unique matrix
#-------------------------------------------------------------------------------
scoreList2scoreMat <- function(listScores,symmetric){

  S <- listScores
  d <- length(S);
  #n <- nrow(S[[1]])
  mat_S <- sapply(1:d, function(q) {mat2Vect(S[[q]], symmetric = symmetric, diag = F)})
  # of dimension N  * d where N = n(n-1)/2 if symmetric, n(n-1) otherwise
  return(mat_S)
}



