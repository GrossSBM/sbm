# R6 Class definition of a Multipartite SBM fit

This class is designed to give a representation and adjust a
Multipartite SBM fitted with GREMLIN.

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
[`MultipartiteSBM`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM.md)
-\> `MultipartiteSBM_fit`

## Active bindings

- `loglik`:

  double: approximation of the log-likelihood (variational lower bound)
  reached

- `ICL`:

  double: value of the integrated classification log-likelihood

- `storedModels`:

  data.frame of all models fitted (and stored) during the optimization

## Methods

### Public methods

- [`MultipartiteSBM_fit$new()`](#method-MultipartiteSBM_fit-initialize)

- [`MultipartiteSBM_fit$optimize()`](#method-MultipartiteSBM_fit-optimize)

- [`MultipartiteSBM_fit$predict()`](#method-MultipartiteSBM_fit-predict)

- [`MultipartiteSBM_fit$setModel()`](#method-MultipartiteSBM_fit-setModel)

- [`MultipartiteSBM_fit$show()`](#method-MultipartiteSBM_fit-show)

- [`MultipartiteSBM_fit$clone()`](#method-MultipartiteSBM_fit-clone)

Inherited methods

- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`MultipartiteSBM$plot()`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM.html#method-plot)
- [`MultipartiteSBM$print()`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM.html#method-print)

------------------------------------------------------------------------

### `MultipartiteSBM_fit$new()`

constructor for Multipartite SBM

#### Usage

    MultipartiteSBM_fit$new(netList)

#### Arguments

- `netList`:

  list of SBM objects

------------------------------------------------------------------------

### `MultipartiteSBM_fit$optimize()`

estimation of multipartiteSBM via GREMLINS

#### Usage

    MultipartiteSBM_fit$optimize(estimOptions)

#### Arguments

- `estimOptions`:

  options for MultipartiteBM

------------------------------------------------------------------------

### `MultipartiteSBM_fit$predict()`

prediction under the currently estimated model

#### Usage

    MultipartiteSBM_fit$predict()

#### Returns

a list of matrices matrix of expected values for each dyad

------------------------------------------------------------------------

### `MultipartiteSBM_fit$setModel()`

method to select a specific model among the ones fitted during the
optimization. Fields of the current MultipartiteSBM_fit will be updated
accordingly.

#### Usage

    MultipartiteSBM_fit$setModel(index)

#### Arguments

- `index`:

  integer, the index of the model to be selected (row number in
  storedModels)

------------------------------------------------------------------------

### `MultipartiteSBM_fit$show()`

show method

#### Usage

    MultipartiteSBM_fit$show(type = "Fit of a Multipartite Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `MultipartiteSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultipartiteSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
