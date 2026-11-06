# R6 Class definition of a Multipartite SBM

R6 Class definition of a Multipartite SBM

## Super class

[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) -\>
`MultipartiteSBM`

## Active bindings

- `dimLabels`:

  vector of characters giving the label of each connected dimension

- `blockProp`:

  list of two vectors of block proportions (aka prior probabilities of
  each block)

- `connectParam`:

  parameters associated to the connectivity of the SBM, e.g. matrix of
  inter/inter block probabilities when model is Bernoulli

- `probMemberships`:

  matrix of estimated probabilities for block memberships for all nodes

- `nbBlocks`:

  : vector with the number of blocks in each FG

- `nbConnectParam`:

  number of parameter used for the connectivity

- `architecture`:

  organization of the multipartite network

- `nbNetworks`:

  number of networks in the multipartite network

- `memberships`:

  list of size 2: vector of memberships in all parts of the network

- `indMemberships`:

  matrix for clustering memberships

## Methods

### Public methods

- [`MultipartiteSBM$new()`](#method-MultipartiteSBM-initialize)

- [`MultipartiteSBM$show()`](#method-MultipartiteSBM-show)

- [`MultipartiteSBM$print()`](#method-MultipartiteSBM-print)

- [`MultipartiteSBM$plot()`](#method-MultipartiteSBM-plot)

- [`MultipartiteSBM$clone()`](#method-MultipartiteSBM-clone)

Inherited methods

- [`SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)

------------------------------------------------------------------------

### `MultipartiteSBM$new()`

constructor for Multipartite SBM

#### Usage

    MultipartiteSBM$new(
      model = character(0),
      architecture = matrix(NA, 0, 2),
      directed = logical(0),
      nbNodes = numeric(0),
      dimLabels = character(0),
      blockProp = list(),
      connectParam = list()
    )

#### Arguments

- `model`:

  character describing the type of model

- `architecture`:

  a 2-column matrix describing interactions between the networks

- `directed`:

  vector of logical: are the network directed or not?

- `nbNodes`:

  number of nodes in each dimension/part of the network

- `dimLabels`:

  labels of each par of the network

- `blockProp`:

  parameters for block proportions (vector of list of vectors)

- `connectParam`:

  parameters of connectivity (vector of list of vectors)

------------------------------------------------------------------------

### `MultipartiteSBM$show()`

print method

#### Usage

    MultipartiteSBM$show(type = "Multipartite Stochastic Block Model")

#### Arguments

- `type`:

  character to tune the displayed name

------------------------------------------------------------------------

### `MultipartiteSBM$print()`

print method

#### Usage

    MultipartiteSBM$print()

------------------------------------------------------------------------

### `MultipartiteSBM$plot()`

plot Multipartite Network

#### Usage

    MultipartiteSBM$plot(
      type = c("data", "expected", "meso"),
      ordered = TRUE,
      plotOptions = list()
    )

#### Arguments

- `type`:

  character for the type of plot: either 'data' (true connection),
  'expected' (fitted connection) or 'meso' (mesoscopic view). Default to
  'data'.

- `ordered`:

  TRUE is the matrices are plotted after reorganization with the blocks.
  Default value = TRUE

- `plotOptions`:

  list of plot options for the mesoscopic view or matrix view

------------------------------------------------------------------------

### `MultipartiteSBM$clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultipartiteSBM$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
