# Model Predictions

Make predictions from an SBM.

## Usage

``` r
# S3 method for class 'SBM'
predict(object, covarList = object$covarList, theta_p0 = 0, ...)
```

## Arguments

- object:

  an R6 object inheriting from class SBM_fit (like SimpleSBM_fit or
  BipartiteSBM_fit)

- covarList:

  a list of covariates. By default, we use the covariates associated
  with the model.

- theta_p0:

  a threshold...

- ...:

  additional parameters for S3 compatibility. Not used

## Value

a matrix of expected values for each dyad
