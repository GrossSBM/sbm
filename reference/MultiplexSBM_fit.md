# R6 Class definition of a Multiplex SBM fit

This class is designed to give a representation and adjust a Multiplex
SBM fitted with GREMLIN.

The list of parameters `estimOptions` essentially tunes the optimization
process and the variational EM algorithm, with the following parameters

- "nbCores" integer for number of cores used. Default is 2 (1 under
  Positron, which does not allow forking)

- "verbosity" integer for verbosity (0, 1). Default is 1

- "nbBlocksRange" List of length the number of functional groups, each
  element supplying the minimal and maximal number of blocks to be
  explored. The names of the list must be the names of the functional
  groups. Default value is from 1 to 10)

- "initBM" Boolean. True if using simple and bipartite SBM as
  initialisations. Default value = TRUE

- "maxiterVEM" Number of max. number of iterations in the VEM. Default
  value = 100

- "maxiterVE" Number of max. number of iterations in the VE. Default
  value = 100

## Super classes

[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) -\>
[`MultipartiteSBM`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM.md)
-\>
[`MultipartiteSBM_fit`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM_fit.md)
-\> `MultiplexSBM_fit`

## Active bindings

- `nbBlocks`:

  vector of size 2: number of blocks (rows, columns)

- `dependentNetwork`:

  : connection parameters in each network

- `storedModels`:

  data.frame of all models fitted (and stored) during the optimization

- `namesLayers`:

  : names of the various Networks

## Methods

### Public methods

- [`MultiplexSBM_fit$new()`](#method-MultiplexSBM_fit-initialize)

- [`MultiplexSBM_fit$optimize()`](#method-MultiplexSBM_fit-optimize)

- [`MultiplexSBM_fit$plot()`](#method-MultiplexSBM_fit-plot)

- [`MultiplexSBM_fit$show()`](#method-MultiplexSBM_fit-show)

- [`MultiplexSBM_fit$predict()`](#method-MultiplexSBM_fit-predict)

- [`MultiplexSBM_fit$clone()`](#method-MultiplexSBM_fit-clone)

Inherited methods

- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`MultipartiteSBM$print()`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM.html#method-print)
- [`MultipartiteSBM_fit$setModel()`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM_fit.html#method-setModel)

------------------------------------------------------------------------

### `MultiplexSBM_fit$new()`

constructor for Multiplex SBM

#### Usage

    MultiplexSBM_fit$new(netList, dependentNet = FALSE)

#### Arguments

- `netList`:

  list of SBM object with

- `dependentNet`:

  boolean indicating whether dependence is assumed between networks
  beyond the common dependence on the latent variables

------------------------------------------------------------------------

### `MultiplexSBM_fit$optimize()`

estimation of multipartiteSBM via GREMLINS

#### Usage

    MultiplexSBM_fit$optimize(estimOptions)

#### Arguments

- `estimOptions`:

  options for MultipartiteBM

------------------------------------------------------------------------

### `MultiplexSBM_fit$plot()`

plot Multiplex Network

#### Usage

    MultiplexSBM_fit$plot(
      type = c("data", "expected", "graphon"),
      ordered = TRUE,
      plotOptions = list()
    )

#### Arguments

- `type`:

  character for the type of plot: either 'data' (true connection),
  'expected' (fitted connection) or 'graphon' (block connectivity of
  each layer on the unit square, blocks sized by the block proportions;
  the marginal connectivity of each layer in the dependent case).
  Default to 'data'.

- `ordered`:

  TRUE is the matrices are plotted after reorganization with the blocks.
  Default value = TRUE

- `plotOptions`:

  list of plot options for the matrix view

------------------------------------------------------------------------

### `MultiplexSBM_fit$show()`

show method

#### Usage

    MultiplexSBM_fit$show(type = "Fit of a Multiplex Stochastic Block Model")

#### Arguments

- `type`:

  character used to specify the type of SBM

------------------------------------------------------------------------

### `MultiplexSBM_fit$predict()`

prediction under the currently estimated model

#### Usage

    MultiplexSBM_fit$predict()

#### Returns

a list of matrices matrix of expected values for each dyad

------------------------------------------------------------------------

### `MultiplexSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultiplexSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
