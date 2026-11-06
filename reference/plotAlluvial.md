# Plot an alluvial plot between clusterings

Plot an alluvial plot between clusterings

## Usage

``` r
plotAlluvial(listMemberships, plotOptions = list())
```

## Arguments

- listMemberships:

  : a list vectors containing the memberships

- plotOptions:

  : a list containing the options for Alluvial plots

## Value

display the alluvial plot, returns invisibly a list with the plotOptions
and the table of frequencies

## Details

The list of parameters `plotOptions` provides the following options

- "curvy" numeric, controls the curvature of the alluvial. Default value
  = 0.3

- "alpha numeric, vector of transparency of the stripes. Default value =
  0.8

- "gap.width" numeric, relative width of inter-category gaps. Default
  value = 0.1

- "col" vector of colors of the stripes. Default value is NULL: the
  stripes are colored according to the clusters of the first clustering

- "border" vector of border colors for the stripes. Default is white

- "reorder" boolean. If TRUE, the clusters of each clustering are
  stacked so as to follow the clusters of the previous one, which limits
  the crossings of the stripes. The labels of the clusters are not
  changed. Default value = TRUE

## Examples

``` r
listMemberships <- list(C1 = rep(c('A','B','C'),each=10),C2 = rep(c(1,2,4),10))
plotAlluvial(listMemberships)
```
