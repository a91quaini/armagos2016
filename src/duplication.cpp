#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

mat duplication_cpp(int n) {

  mat A(n, n); A.ones(); A = trimatl(A);
  int m = n * (n + 1) / 2;
  A(find(A)) = linspace(1, m, m);
  A = symmatl(A);
  colvec j = vectorise(A);

  mat D = zeros(n * n, m);

  for (int r=0; r<n*n; ++r ) {
    D(r, j(r) - 1) = 1;
  }

  return D;

}
