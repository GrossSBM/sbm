library(testthat)
library(sbm)
library(aricode)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

test_check("sbm")
