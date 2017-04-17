#include <RcppArmadillo.h>
#include "duplication.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
colvec UnbiasedEst_cpp(colvec nu, mat b, mat V, mat Vxi, colvec tau, int n, int T_obs, int K) {

  colvec c = ones(K + 1);
  c(span(1, K)) = -nu;
  colvec Bs1 = zeros(K);
  colvec v = zeros(n);
  colvec w = zeros(n);
  mat E2 = join_cols(zeros(1, K), eye(K, K));
  mat D = duplication_cpp(K + 1);

  for(int i=0; i<n; ++i) {

    mat W = reshape(D * V.row(i).t(), K + 1, K + 1);
    v(i) = conv_to<double>::from(tau(i) * c.t() * W * c);
    w(i) = 1 / v(i);

    mat Wxi = reshape(D * Vxi.row(i).t(), K + 1, K + 1);
    Bs1 = Bs1 + w(i) * tau(i) * E2.t() * Wxi * c;

  }

  mat Wn = diagmat(w);
  mat Qb = b.t() * Wn * b / n;
  colvec Bnu = solve(Qb, Bs1 / n);
  colvec nuB = nu - Bnu / T_obs;

  return nuB;

}
