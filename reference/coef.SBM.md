# Extract model coefficients

Extracts model coefficients from objects with class
[`SBM`](https://grosssbm.github.io/sbm/reference/SBM.md) and children
([`SimpleSBM_fit`](https://grosssbm.github.io/sbm/reference/SimpleSBM_fit.md),
[`BipartiteSBM_fit`](https://grosssbm.github.io/sbm/reference/BipartiteSBM_fit.md))

## Usage

``` r
# S3 method for class 'SBM'
coef(object, type = c("connectivity", "block", "covariates"), ...)
```

## Arguments

- object:

  an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or
  BipartiteSBM_fit)

- type:

  type of parameter that should be extracted. Either 'block' for
  \$\$\pi\$\$, 'connectivity' for \$\$\theta\$\$, or "covariates" for
  \$\$\beta\$\$. Default is 'connectivity'.

- ...:

  additional parameters for S3 compatibility. Not used

## Value

vector or list of parameters.
