# R6 virtual class for SBM representation (mother class of SimpleSBM, BipartiteSBM, MultipartiteSBM)

R6 virtual class for SBM representation (mother class of SimpleSBM,
BipartiteSBM, MultipartiteSBM)

## Active bindings

- `modelName`:

  character, the family of model for the distribution of the edges

- `directed`:

  mode of the network data (directed or not or not applicable)

- `dimLabels`:

  vector or list of characters, the label of each dimension

- `nbNodes`:

  vector describing the number of the successive elements connecting the
  network

- `nbCovariates`:

  integer, the number of covariates

- `connectParam`:

  parameters associated to the connectivity of the SBM, e.g. matrix of
  inter/inter block probabilities when model is Bernoulli

- `covarParam`:

  vector of regression parameters associated with the covariates.

- `covarList`:

  list of matrices of covariates

- `covarArray`:

  the array of covariates

- `covarEffect`:

  effect of covariates

- `nbNodesCovariates`:

  the dimension of the nodes covariates

- `nodesCovariates`:

  the list of nodes covariates

- `nodesCovariatesParam`:

  the list of nodes covariates parameters

- `networkData`:

  the network data (adjacency or incidence matrix or list of such
  object)

- `expectation`:

  expected values of connection under the current model

## Methods

### Public methods

- [`SBM$new()`](#method-SBM-initialize)

- [`SBM$rNetwork()`](#method-SBM-rNetwork)

- [`SBM$show()`](#method-SBM-show)

- [`SBM$print()`](#method-SBM-print)

- [`SBM$clone()`](#method-SBM-clone)

------------------------------------------------------------------------

### `SBM$new()`

constructor for SBM

#### Usage

    SBM$new(
      model = vector("character", 0),
      directed = vector("logical", 0),
      dimension = vector("numeric", 0),
      dimLabels = vector("character", 0),
      blockProp = vector("numeric", 0),
      connectParam = vector("list", 0),
      covarParam = numeric(length(covarList)),
      covarList = list(),
      nodesCovarList = list(),
      nodesCovarParam = list()
    )

#### Arguments

- `model`:

  character describing the type of model

- `directed`:

  logical describing if the network data is directed or not

- `dimension`:

  dimension of the network data

- `dimLabels`:

  labels of each dimension

- `blockProp`:

  parameters for block proportions (vector or list of vectors)

- `connectParam`:

  list of parameters for connectivity

- `covarParam`:

  optional vector of covariates effect

- `covarList`:

  optional list of covariates data

- `nodesCovarList`:

  optional list of covariates on nodes (rows and cols),

- `nodesCovarParam`:

  optional list of vector of node covariates effects on the clustering

------------------------------------------------------------------------

### `SBM$rNetwork()`

a method to sample a network data for the current SBM (blocks and edges)

#### Usage

    SBM$rNetwork(store = FALSE)

#### Arguments

- `store`:

  should the sampled network be stored (and overwrite the existing
  data)? Default to FALSE

#### Returns

a list with the sampled block and network

------------------------------------------------------------------------

### `SBM$show()`

print method

#### Usage

    SBM$show(type = "Stochastic Block Model")

#### Arguments

- `type`:

  character to tune the displayed name

------------------------------------------------------------------------

### `SBM$print()`

print method

#### Usage

    SBM$print()

------------------------------------------------------------------------

### `SBM$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SBM$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
