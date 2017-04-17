#include <RcppArmadillo.h>
#include <string.h>
#include "duplication.h"
#include "UnbiasedEst.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List CR_reg_cpp(colvec a, mat b, mat V, mat Vxi, colvec tau, mat F_mat,
  int n, int nchi, int T_obs, int K, std::string dataset) {

  colvec nu1 = solve(b.t() * b, b.t() * a);
  colvec c = ones(K + 1);
  c(span(1, K)) = -nu1;
  mat sum1 = zeros(K, K);
  colvec sum2 = zeros(K);
  mat D = duplication_cpp(K + 1);

  colvec vecV;
  double w;
  for (int i=0; i<nchi; ++i) {

    vecV = D * V.row(i).t();
    mat W = reshape(vecV, K + 1, K + 1);
    double v = conv_to<double>::from(tau(i) * c.t() * W * c);
    w = 1 / v;
    sum1 = sum1 + w * b.row(i).t() * b.row(i);
    sum2 = sum2 + w * b.row(i).t() * a(i);

  }

  colvec nu = solve(sum1, sum2);
  colvec lambda(nu.n_elem);

  if ( dataset.compare("individual") != 0) {

    lambda = nu + mean(F_mat).t();

  } else {

  nu = UnbiasedEst_cpp(nu, b, V, Vxi, tau, nchi, T_obs, K);
  lambda = nu + mean(F_mat).t();

  }

  colvec RiskPremia = lambda * 12 * 100 / 10;
  colvec Yhat(nchi); Yhat.zeros();

  for (int i=0; i<nchi; ++i) {
      Yhat(i) = conv_to<double>::from(b.row(i) * RiskPremia);
  }

  return List::create(
    Named("RiskPremia") = RiskPremia
    , Named("Yhat") = Yhat
  );
}
