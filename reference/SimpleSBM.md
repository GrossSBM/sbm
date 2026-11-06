# R6 class for Simple SBM

R6 class for Simple SBM

## Super class

[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) -\> `SimpleSBM`

## Active bindings

- `dimLabels`:

  a single character giving the label of the nodes

- `blockProp`:

  vector of block proportions (aka prior probabilities of each block)

- `connectParam`:

  parameters associated to the connectivity of the SBM, e.g. matrix of
  inter/inter block probabilities when model is Bernoulli

- `probMemberships`:

  matrix of estimated probabilities for block memberships for all nodes

- `mask`:

  Mask for the (potential) NAs in the adjacency matrix

- `nbBlocks`:

  number of blocks

- `nbDyads`:

  number of dyads (potential edges in the network)

- `nbConnectParam`:

  number of parameter used for the connectivity

- `nodesCovariates`:

  matrix of covariates on nodes

- `nodesCovariatesParam`:

  the matrix of nodes covariates parameters

- `nbNodesCovariates`:

  number of covariates on nodes

- `memberships`:

  vector of clustering

- `indMemberships`:

  matrix for clustering memberships

## Methods

### Public methods

- [`SimpleSBM$new()`](#method-SimpleSBM-initialize)

- [`SimpleSBM$rMemberships()`](#method-SimpleSBM-rMemberships)

- [`SimpleSBM$rEdges()`](#method-SimpleSBM-rEdges)

- [`SimpleSBM$predict()`](#method-SimpleSBM-predict)

- [`SimpleSBM$show()`](#method-SimpleSBM-show)

- [`SimpleSBM$plot()`](#method-SimpleSBM-plot)

- [`SimpleSBM$clone()`](#method-SimpleSBM-clone)

Inherited methods

- [`SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)

------------------------------------------------------------------------

### `SimpleSBM$new()`

constructor for SBM

#### Usage

    SimpleSBM$new(
      model,
      nbNodes,
      directed,
      blockProp,
      connectParam,
      dimLabels = c("node"),
      covarParam = numeric(length(covarList)),
      covarList = list(),
      nodesCovarList = list(node = matrix(0, 0, 0)),
      nodesCovarParam = list(node = matrix(0, 0, 0))
    )

#### Arguments

- `model`:

  character describing the type of model

- `nbNodes`:

  number of nodes in the network

- `directed`:

  logical, directed network or not.

- `blockProp`:

  parameters for block proportions (vector of list of vectors)

- `connectParam`:

  list of parameters for connectivity with a matrix of means 'mean' and
  an optional scalar for the variance 'var'. The size of mu must match
  `blockProp` length

- `dimLabels`:

  optional label for the node (default is "nodeName")

- `covarParam`:

  optional vector of covariates effect

- `covarList`:

  optional list of covariates data

- `nodesCovarList`:

  optional list of one matrix with nodes covariates

- `nodesCovarParam`:

  optional list of one matrix of parameters of effects of nodes
  covariates on clustering

------------------------------------------------------------------------

### `SimpleSBM$rMemberships()`

a method to sample new block memberships for the current SBM

#### Usage

    SimpleSBM$rMemberships(store = FALSE)

#### Arguments

- `store`:

  should the sampled blocks be stored (and overwrite the existing data)?
  Default to FALSE

#### Returns

the sampled blocks

------------------------------------------------------------------------

### `SimpleSBM$rEdges()`

a method to sample a network data (edges) for the current SBM

#### Usage

    SimpleSBM$rEdges(store = FALSE)

#### Arguments

- `store`:

  should the sampled edges be stored (and overwrite the existing data)?
  Default to FALSE

#### Returns

the sampled network

------------------------------------------------------------------------

### `SimpleSBM$predict()`

prediction under the currently parameters

#### Usage

    SimpleSBM$predict(covarList = self$covarList, theta_p0 = 0)

#### Arguments

- `covarList`:

  a list of covariates. By default, we use the covariates with which the
  model was estimated

- `theta_p0`:

  a threshold...

#### Returns

a matrix of expected values for each dyad

------------------------------------------------------------------------

### `SimpleSBM$show()`

show method

#### Usage

    SimpleSBM$show(type = "Simple Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `SimpleSBM$plot()`

basic matrix plot method for SimpleSBM object or mesoscopic plot

#### Usage

    SimpleSBM$plot(
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

### `SimpleSBM$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
