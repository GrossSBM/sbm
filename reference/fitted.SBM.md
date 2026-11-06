# Extract model fitted values

Extracts fitted values for object with class
([`SimpleSBM_fit`](https://grosssbm.github.io/sbm/reference/SimpleSBM_fit.md),
[`BipartiteSBM_fit`](https://grosssbm.github.io/sbm/reference/BipartiteSBM_fit.md))
or
[`multipartitepartiteSBM_fit`](https://grosssbm.github.io/sbm/reference/MultipartiteSBM_fit.md))

## Usage

``` r
# S3 method for class 'SBM'
fitted(object, ...)
```

## Arguments

- object:

  an R6 object inheriting from SimpleSBM_fit, BipartiteSBM_fit or
  MultipartiteSBM_fit

- ...:

  additional parameters for S3 compatibility. Not used

## Value

a matrix of expected fitted values for each dyad
