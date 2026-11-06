# Define a network

Define a network

## Usage

``` r
defineSBM(
  netMat,
  model = "bernoulli",
  type = ifelse(ncol(netMat) == nrow(netMat), "simple", "bipartite"),
  directed = !isSymmetric(netMat),
  dimLabels = c(row = "row", col = "col"),
  covariates = list()
)
```

## Arguments

- netMat:

  a matrix describing the network: either an adjacency (square) or
  incidence matrix with possibly weighted entries.

- model:

  character describing the model for the relation between nodes
  (`'bernoulli'`, `'poisson'`, `'gaussian'`, ...). Default is
  `'bernoulli'`.

- type:

  Type of the matrix, choice between 'simple' and 'bipartite'

- directed:

  logical: is the network directed or not? Only relevant when `type` is
  `'Simple'`. Default is `TRUE` if `netMat` is symmetric, `FALSE`
  otherwise

- dimLabels:

  an optional vector of labels for each dimension (in row, in column).
  Default value = c('row' = row,'col'= col)

- covariates:

  a list of matrices with same dimension as mat describing covariates at
  the edge level. No covariate per Default.

## Value

an object SimpleSBM or BipartiteSBM with the informations required to
define a future multipartite network

## Examples

``` r
A <- matrix(rbinom(100,1,.2), 10, 10)
myNet <- defineSBM(A, "bernoulli", "simple", TRUE, "Actor")
```
