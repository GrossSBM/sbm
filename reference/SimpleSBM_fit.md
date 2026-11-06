# R6 Class definition of a Simple SBM fit

This class is designed to give a representation and adjust an SBM fitted
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
[`SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.md) -\>
`SimpleSBM_fit`

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

- [`SimpleSBM_fit$new()`](#method-SimpleSBM_fit-initialize)

- [`SimpleSBM_fit$optimize()`](#method-SimpleSBM_fit-optimize)

- [`SimpleSBM_fit$setModel()`](#method-SimpleSBM_fit-setModel)

- [`SimpleSBM_fit$reorder()`](#method-SimpleSBM_fit-reorder)

- [`SimpleSBM_fit$show()`](#method-SimpleSBM_fit-show)

- [`SimpleSBM_fit$clone()`](#method-SimpleSBM_fit-clone)

Inherited methods

- [`SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`SimpleSBM$plot()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-plot)
- [`SimpleSBM$predict()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-predict)
- [`SimpleSBM$rEdges()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rEdges)
- [`SimpleSBM$rMemberships()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rMemberships)

------------------------------------------------------------------------

### `SimpleSBM_fit$new()`

constructor for a Simple SBM fit

#### Usage

    SimpleSBM_fit$new(
      adjacencyMatrix,
      model,
      directed,
      dimLabels = c(node = "nodeName"),
      covarList = list(),
      nodesCovarList = list(matrix(nrow = 0, ncol = 0))
    )

#### Arguments

- `adjacencyMatrix`:

  square (weighted) matrix

- `model`:

  character (`'bernoulli'`, `'poisson'`, `'gaussian'`)

- `directed`:

  logical, directed network or not. In not, `adjacencyMatrix` must be
  symmetric.

- `dimLabels`:

  list of labels of each dimension (in row, in columns)

- `covarList`:

  and optional list of covariates, each of whom must have the same
  dimension as `adjacencyMatrix`

- `nodesCovarList`:

  optional list of one matrix with nodes covariates

------------------------------------------------------------------------

### `SimpleSBM_fit$optimize()`

function to perform optimization

#### Usage

    SimpleSBM_fit$optimize(estimOptions = list())

#### Arguments

- `estimOptions`:

  a list of parameters controlling the inference algorithm and model
  selection. See details.

------------------------------------------------------------------------

### `SimpleSBM_fit$setModel()`

method to select a specific model among the ones fitted during the
optimization. Fields of the current SBM_fit will be updated accordingly.

#### Usage

    SimpleSBM_fit$setModel(index)

#### Arguments

- `index`:

  integer, the index of the model to be selected (row number in
  storedModels)

------------------------------------------------------------------------

### `SimpleSBM_fit$reorder()`

permute group labels by order of decreasing probability

#### Usage

    SimpleSBM_fit$reorder()

------------------------------------------------------------------------

### `SimpleSBM_fit$show()`

show method

#### Usage

    SimpleSBM_fit$show(type = "Fit of a Simple Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `SimpleSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
