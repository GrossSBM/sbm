# R6 class for Bipartite SBM

R6 class for Bipartite SBM

## Super class

[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) -\>
`BipartiteSBM`

## Active bindings

- `dimLabels`:

  vector of two characters giving the label of each connected dimension
  (row, col)

- `blockProp`:

  list of two vectors of block proportions (aka prior probabilities of
  each block)

- `connectParam`:

  parameters associated to the connectivity of the SBM, e.g. matrix of
  inter/inter block probabilities when model is Bernoulli

- `probMemberships`:

  matrix of estimated probabilities for block memberships for all nodes

- `nbBlocks`:

  vector of size 2: number of blocks (rows, columns)

- `mask`:

  Mask for the (potential) NAs in the adjacency matrix

- `nbDyads`:

  number of dyads (potential edges in the network)

- `nbConnectParam`:

  number of parameter used for the connectivity

- `memberships`:

  list of size 2: vector of memberships in row, in column.

- `indMemberships`:

  matrix for clustering memberships

## Methods

### Public methods

- [`BipartiteSBM$new()`](#method-BipartiteSBM-initialize)

- [`BipartiteSBM$rMemberships()`](#method-BipartiteSBM-rMemberships)

- [`BipartiteSBM$rEdges()`](#method-BipartiteSBM-rEdges)

- [`BipartiteSBM$predict()`](#method-BipartiteSBM-predict)

- [`BipartiteSBM$show()`](#method-BipartiteSBM-show)

- [`BipartiteSBM$plot()`](#method-BipartiteSBM-plot)

- [`BipartiteSBM$clone()`](#method-BipartiteSBM-clone)

Inherited methods

- [`SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)

------------------------------------------------------------------------

### `BipartiteSBM$new()`

constructor for SBM

#### Usage

    BipartiteSBM$new(
      model,
      nbNodes,
      blockProp = vector("list", 2),
      connectParam = list(mean = matrix(0, 0, 0)),
      dimLabels = c(row = "row", col = "col"),
      covarParam = numeric(length(covarList)),
      covarList = list(),
      nodesCovarList = vector("list", 2),
      nodesCovarParam = vector("list", 2)
    )

#### Arguments

- `model`:

  character describing the type of model

- `nbNodes`:

  number of nodes in each dimension of the network

- `blockProp`:

  parameters for block proportions (vector of list of vectors)

- `connectParam`:

  list of parameters for connectivity with a matrix of means 'mean' and
  an optional scalar for the variance 'var'. The dimensions of mu must
  match `blockProp` lengths

- `dimLabels`:

  optional labels of each dimension (in row, in column)

- `covarParam`:

  optional vector of covariates effect

- `covarList`:

  optional list of covariates data

- `nodesCovarList`:

  optional list of two matrices with nodes covariates

- `nodesCovarParam`:

  optional list of two matrices of parameters of effects of nodes
  covariates on clustering

------------------------------------------------------------------------

### `BipartiteSBM$rMemberships()`

a method to sample new block memberships for the current SBM

#### Usage

    BipartiteSBM$rMemberships(store = FALSE)

#### Arguments

- `store`:

  should the sampled blocks be stored (and overwrite the existing data)?
  Default to FALSE

#### Returns

the sampled blocks

------------------------------------------------------------------------

### `BipartiteSBM$rEdges()`

a method to sample a network data (edges) for the current SBM

#### Usage

    BipartiteSBM$rEdges(store = FALSE)

#### Arguments

- `store`:

  should the sampled edges be stored (and overwrite the existing data)?
  Default to FALSE

#### Returns

the sampled network

------------------------------------------------------------------------

### `BipartiteSBM$predict()`

prediction under the current parameters

#### Usage

    BipartiteSBM$predict(covarList = self$covarList, theta_p0 = 0)

#### Arguments

- `covarList`:

  a list of covariates. By default, we use the covariates with which the
  model was estimated.

- `theta_p0`:

  double for thresholding...

------------------------------------------------------------------------

### `BipartiteSBM$show()`

show method

#### Usage

    BipartiteSBM$show(type = "Bipartite Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `BipartiteSBM$plot()`

basic matrix plot method for BipartiteSBM object or mesoscopic plot

#### Usage

    BipartiteSBM$plot(
      type = c("data", "expected", "graphon", "meso"),
      ordered = TRUE,
      plotOptions = list()
    )

#### Arguments

- `type`:

  character for the type of plot: either 'data' (true connection),
  'expected' (fitted connection), 'graphon' (block connectivity on the
  unit square, blocks sized by the block proportions) or 'meso'
  (mesoscopic view). Default to 'data'.

- `ordered`:

  logical: should the rows and columns be reordered according to the
  clustering? Default to `TRUE`.

- `plotOptions`:

  list with the parameters for the plot. See help of the corresponding
  S3 method for details.

#### Returns

a ggplot2 object

------------------------------------------------------------------------

### `BipartiteSBM$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BipartiteSBM$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
