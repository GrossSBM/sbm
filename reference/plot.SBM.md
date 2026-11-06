# SBM Plot

Basic matrix plot method for SBM object or mesoscopic view

## Usage

``` r
# S3 method for class 'SBM'
plot(
  x,
  type = c("data", "expected", "graphon", "meso"),
  ordered = TRUE,
  plotOptions = list(),
  ...
)
```

## Arguments

- x:

  an object inheriting from class SBM

- type:

  character for the type of plot: either 'data' (true connection),
  'expected' (fitted connection), 'graphon' (block connectivity on the
  unit square, blocks sized by the block proportions) or 'meso'
  (mesoscopic). Default to 'data'. 'graphon' is available for simple,
  bipartite and multiplex SBM (one panel per layer).

- ordered:

  logical: should the rows and columns be ordered according to the
  clustering? Default to `TRUE` (not taken into account for 'meso').

- plotOptions:

  list with parameters for 'meso' type plot and data type plot. Details
  are given below

- ...:

  additional parameters for S3 compatibility. Not used

## Value

a ggplot2 object for the `'data'`, `'expected'` and `'graphon'`, and
`'meso'`. For `'meso'`, the coordinates of the blocks are in the
attribute `"blockLayout"` of the plot (a matrix with one row per block,
named by the labels of the blocks): they can be modified and given back
in `plotOptions$layout`, for instance
`lay <- attr(plot(mySBM, type = "meso"), "blockLayout")`, then
`lay["P1", ] <- c(0, 2)` and
`plot(mySBM, type = "meso", plotOptions = list(layout = lay))`

## Details

For simple and bipartite SBM, the mesoscopic view draws the graphon as a
network: each block is a node whose area is proportional to the block
proportion (the largest block of each type of nodes having the same
size), and the connectivity between two blocks is an edge whose width
and color grow with the connectivity, with the colors of the graphon.
The connectivity within a block is drawn as a loop. Edges are curved
arrows for directed networks. The list of parameters `plotOptions` for
this mesoscopic plot is:

- "seed": :

  seed to control the layout (free layout). Default value is NULL

- "title": :

  character string for the title. Default value is NULL

- "layout": :

  "free" (force-directed: the blocks are placed according to their
  connectivities), "groups" (one line of blocks per type of nodes, for
  bipartite SBM) or a matrix with the coordinates of the blocks (one row
  per block, rows then columns for the bipartite SBM). Default value is
  NULL: "free" for simple SBM, "groups" for bipartite SBM

- "vertex.color": :

  Color of the blocks, one per type of nodes. Default values are taken
  in a hue palette without the reds of the connectivities

- "vertex.frame.color": :

  Color of the border of the blocks. Default value is NA (no border)

- "vertex.shape": :

  "circle", "square", "diamond" or "triangle", one per type of nodes.
  Default value is "circle": the types of nodes are distinguished by
  their colors

- "vertex.size": :

  Factor multiplying the size of the blocks, one per type of nodes.
  Default value is 1

- "vertex.label.name": :

  Prefix of the labels of the blocks, one per type of nodes. Default
  value is NULL: the number of the block for a simple SBM, the shortest
  prefixes of the labels of the dimensions which distinguish them for a
  bipartite SBM

- "vertex.label.color": :

  Default value is "black"

- "vertex.label.size": :

  Size of the labels. Default value is 3.5

- "vertex.label.inside": :

  The label of a block is drawn inside it when its radius is at least
  this fraction of the radius of the largest block (for a label of 2
  characters, the fraction being proportional to the number of
  characters), outside it otherwise. Default value is 0.5. Set it to 0
  to draw all the labels inside the blocks

- "withinBlock": :

  Representation of the connectivity within the blocks for simple SBM:
  "loop" (a loop, as the other connections) or "fill" (the color of the
  block). Default value is "loop"

- "edge.threshold": :

  Threshold under which the connectivity is not plotted. Default value
  is -Inf

- "edge.width": :

  Factor multiplying the width of the edges. Default value is 1

- "edge.curved": :

  Curvature of the edges for directed networks. Default value is 0.2

- "edge.arrow.size": :

  Factor multiplying the size of the arrows for directed networks.
  Default value is 1

- "block.color", "blockAlpha.min", "blockAlpha.max", "limits": :

  Colors of the connectivities, as for the graphon (see below)

For multipartite SBM, the list of parameters `plotOptions` for the
mesoscopic plot (drawn with ggplot2, the free layout being computed with
igraph) is:

- "seed": :

  seed to control the free layout. Default value is NULL

- "title": :

  character string for the title. Default value is NULL

- "layout": :

  "free" (the blocks are placed according to their connectivities,
  whatever their functional group), "groups" (one row of blocks per
  functional group) or a matrix with the coordinates of the blocks (one
  row per block, in the order of the functional groups). Default value
  is "free"

- "vertex.color": :

  Colors of the functional groups. Default values are taken in a hue
  palette without the reds of the connectivities

- "vertex.frame.color": :

  Color of the border of the blocks. Default value is NA (no border)

- "vertex.shape": :

  "circle", "square", "diamond" or "triangle", one per functional group.
  Default value is "circle": the functional groups are distinguished by
  their colors

- "vertex.size": :

  Size of the node (default is 1)

- "vertex.size2": :

  The second size of the node (e.g. for a rectangle)

- "vertex.label.name": :

  Names of the vertices. Default value is the label of the nodes

- "vertex.label.color": :

  Default value is "black"

- "vertex.label.font": :

  Default value is 2. Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5
  symbol

- "vertex.label.cex": :

  Font size (multiplication factor, device-dependent).Default value is
  0.9.

- "vertex.label.dist": :

  Distance between the label and the vertex. Default value is 0

- "vertex.label.degree": :

  The position of the label in relation to the vertex. default value is
  0

- "vertex.label.inside": :

  The label of a block is drawn inside it when its radius is at least
  this fraction of the radius of the largest block (for a label of 2
  characters, the fraction being proportional to the number of
  characters), outside it otherwise. Default value is 0.5

- "edge.threshold": :

  Threshold under which the edge is not plotted. Default value is = -Inf

- "edge.color": :

  Default value is "gray"

- "edge.width": :

  Factor parameter. Default value is 1

- "edge.arrow.size": :

  Default value is 1

- "edge.arrow.width": :

  Default value is 1

- "edge.lty": :

  Line type, could be 0 or "blank", 1 or "solid", 2 or "dashed", 3 or
  "dotted", 4 or "dotdash", 5 or "longdash", 6 or "twodash". Default
  value is "solid"

- "edge.curved": :

  Default value is = 0.3.

For type = 'data' or 'expected plot', the list of parameters
`plotOptions` is

- "legend": :

  Boolean. Set TRUE if you want to see the legend. Default value is
  FALSE

- "legend.title"::

  Boolean. Set TRUE if you want to print the title of the legend.
  Default value is FALSE

- "legend.position"::

  Position of the legend. Possible values are 'bottom',
  'top','left,'right'. Default value is 'bottom'

- "rowNames"::

  Set true if the rownames must be plotted. Default value is FALSE

- "colNames"::

  Set true if the colNames must be plotted. Default value is FALSE

- "names.size"::

  Numeric. Size (in pts) of the row and column names. Default value is
  NULL: the size decreases with the number of nodes, so that the names
  do not overlap on a 5 inches plot

- "line.color": :

  Chain of character. The color of the lines to separate groups if a
  clustering is provided. Default value is red

- "line.width": :

  Numeric. Width of the lines to separate groups. Default value is NULL,
  automatically chosen

- "title": :

  Chain of character. Title of the plot. Default value is NULL

- "blockBackground": :

  Boolean. Color the background of each block (when `ordered = TRUE`),
  with an opacity growing with the mean value of the entries in the
  block. Default value is TRUE for `type = 'data'` and FALSE for
  `type = 'expected'`

- "graphonColors": :

  Boolean. Color the entries with the colors of the graphon (see below)
  instead of gray levels. Default value is TRUE for `type = 'expected'`
  and FALSE for `type = 'data'`

- "block.color": :

  Chain of character. The color of the block background. Default value
  is red

- "blockAlpha.min", "blockAlpha.max": :

  Numeric. Opacity of the least and most dense blocks. Default values
  are 0.1 and 0.6

For type = 'graphon', the blocks are drawn with the colors of the block
background of the data plot. The axes are in number of nodes, so that
the aspect ratio of the plot reflects the dimensions of the network. The
list of parameters `plotOptions` is

- "legend", "legend.title", "legend.position", "title", "line.color",
  "line.width": :

  As for the 'data' plot

- "block.color", "blockAlpha.min", "blockAlpha.max": :

  As for the 'data' plot

- "limits": :

  Numeric vector of length 2. Range of the connectivity values mapped to
  the colors. Default value is NULL: from 0 to the maximal value (or the
  range for negative values). Set it to compare several graphons on the
  same color scale, for instance c(0, 1) for Bernoulli models

The graphon is drawn from the connectivity parameters only, without the
effect of the covariates.
