# R6 Class definition of an Bipartite SBM fit

This class is designed to give a representation and adjust an LBM fitted
with blockmodels.

## Details

The list of parameters `estimOptions` essentially tunes the optimization
process and the variational EM algorithm, with the following parameters

- "nbCores integer for number of cores used. Default is 2 (1 under
  Positron, which does not allow forking)

- "verbosity" integer for verbosity (0, 1). Default is 1

- "plot" boolean, should the ICL be dynamically plotted or not. Default
  is TRUE, FALSE when "verbosity" is 0

- "exploreFactor" control the exploration of the number of groups

- "exploreMin" explore at least until exploreMin even if the exploration
  factor rule is achieved. Default 4. See the package blockmodels for
  details.

- "exploreMax" Stop exploration at exploreMax even if the exploration
  factor rule is not achieved. Default Inf. See the package blockmodels
  for details.

- "nbBlocksRange" minimal and maximal number or blocks explored

- "fast" logical: should approximation be used for Bernoulli model with
  covariates. Default to `TRUE`

## Super classes

[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) -\>
[`BipartiteSBM`](https://grosssbm.github.io/sbm/reference/BipartiteSBM.md)
-\> `BipartiteSBM_fit`

## Active bindings

- `loglik`:

  double: approximation of the log-likelihood (variational lower bound)
  reached

- `ICL`:

  double: value of the integrated classification log-likelihood

- `penalty`:

  double, value of the penalty term in ICL

- `entropy`:

  double, value of the entropy due to the clustering distribution

- `storedModels`:

  data.frame of all models fitted (and stored) during the optimization

## Methods

### Public methods

- [`BipartiteSBM_fit$new()`](#method-BipartiteSBM_fit-initialize)

- [`BipartiteSBM_fit$optimize()`](#method-BipartiteSBM_fit-optimize)

- [`BipartiteSBM_fit$setModel()`](#method-BipartiteSBM_fit-setModel)

- [`BipartiteSBM_fit$reorder()`](#method-BipartiteSBM_fit-reorder)

- [`BipartiteSBM_fit$show()`](#method-BipartiteSBM_fit-show)

- [`BipartiteSBM_fit$clone()`](#method-BipartiteSBM_fit-clone)

Inherited methods

- [`SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`BipartiteSBM$plot()`](https://grosssbm.github.io/sbm/reference/BipartiteSBM.html#method-plot)
- [`BipartiteSBM$predict()`](https://grosssbm.github.io/sbm/reference/BipartiteSBM.html#method-predict)
- [`BipartiteSBM$rEdges()`](https://grosssbm.github.io/sbm/reference/BipartiteSBM.html#method-rEdges)
- [`BipartiteSBM$rMemberships()`](https://grosssbm.github.io/sbm/reference/BipartiteSBM.html#method-rMemberships)

------------------------------------------------------------------------

### `BipartiteSBM_fit$new()`

constructor for a Bipartite SBM fit

#### Usage

    BipartiteSBM_fit$new(
      incidenceMatrix,
      model,
      dimLabels = c(row = "row", col = "col"),
      covarList = list(),
      nodesCovarList = vector("list", 2)
    )

#### Arguments

- `incidenceMatrix`:

  rectangular (weighted) matrix

- `model`:

  character (`'bernoulli'`, `'poisson'`, `'gaussian'`)

- `dimLabels`:

  labels of each dimension (in row, in columns)

- `covarList`:

  an optional list of covariates, each of whom must have the same
  dimension as `incidenceMatrix`

- `nodesCovarList`:

  an optional list of two matrices of covariates, each of whom must have
  the same dimension as `incidenceMatrix`

------------------------------------------------------------------------

### `BipartiteSBM_fit$optimize()`

function to perform optimization

#### Usage

    BipartiteSBM_fit$optimize(estimOptions = list())

#### Arguments

- `estimOptions`:

  a list of parameters controlling the inference algorithm and model
  selection. See details.

------------------------------------------------------------------------

### `BipartiteSBM_fit$setModel()`

method to select a specific model among the ones fitted during the
optimization. Fields of the current SBM_fit will be updated accordingly.

#### Usage

    BipartiteSBM_fit$setModel(index)

#### Arguments

- `index`:

  integer, the index of the model to be selected (row number in
  storedModels)

------------------------------------------------------------------------

### `BipartiteSBM_fit$reorder()`

permute group labels by order of decreasing probability

#### Usage

    BipartiteSBM_fit$reorder()

------------------------------------------------------------------------

### `BipartiteSBM_fit$show()`

show method

#### Usage

    BipartiteSBM_fit$show(type = "Fit of a Bipartite Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `BipartiteSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BipartiteSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
