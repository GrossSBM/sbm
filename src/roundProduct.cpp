#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix roundProduct(arma::cube phi, arma::vec beta) {

  int N = phi.n_rows;
  arma::mat M = arma::zeros<arma::mat>(N,N);

  for (unsigned int k = 0; k < beta.size(); k++) {
    M += phi.slice(k) * beta[k];
  }

  return Rcpp::wrap(M);
}