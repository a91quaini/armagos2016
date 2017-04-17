#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

mat commutation_cpp(int n, int m) {

  colvec vec_eye_n = vectorise(eye(n, n));
  mat k = reshape(kron(vec_eye_n, eye(m, m)), n * m, n * m);
  return k;

}
