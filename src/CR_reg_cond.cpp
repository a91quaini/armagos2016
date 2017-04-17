#include <RcppArmadillo.h>
#include "duplication.h"
#include "commutation.h"
#include "Cmatrix.h"
#include "UnbiasedEst_cond.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List CR_reg_cond_cpp(int n, int nchi, int T_obs, int K, int p, int q, int d1, int d2, int d, mat beta1,
  mat beta2, mat V, colvec tau, mat Z, cube Zi, mat F_mat, std::string dataset) {

  mat D = duplication_cpp(p);
  mat Dplus = pinv(D);
  mat Wpq = commutation_cpp(p, q);

  mat sum1 = zeros(K * p, K * p);
  colvec sum2 = zeros(K * p);
  mat vecbeta3 = zeros(d1 * K * p, nchi);

  cube Bi(K, p ,nchi);
  cube Ci(K, q, nchi);
  cube b_timevar(T_obs - 1, K, nchi);

    for(int i=0; i<nchi; ++i) {

    Bi.slice(i) = reshape(beta2(span(0, K * p - 1), i), p, K).t();
    Ci.slice(i) = reshape(beta2(span(K * p, d2 - 1), i), q, K).t();
    b_timevar.slice(i) = Z.rows(0, T_obs - 2) * Bi.slice(i).t() + Zi.slice(i).rows(0, T_obs - 2) * Ci.slice(i).t();

    mat beta3 = join_cols(
      Dplus * kron(Bi.slice(i).t(), eye(p, p)),
      Wpq * kron(Ci.slice(i).t(), eye(p, p))
    );
    sum1 = sum1 + beta3.t() * beta3;
    sum2 = sum2 + beta3.t() * beta1.col(i);
    vecbeta3.col(i) = vectorise(beta3.t());

  }

  colvec nu1 = solve(sum1, sum2);
  List Cmatrix_out = Cmatrix_cpp(nu1, p, q, K, d1, d2);
  mat C = Cmatrix_out["C"];

  sum1 = zeros(K * p, K * p);
  sum2 = zeros(K * p);
  D = duplication_cpp(d);

  for(int i=0; i<nchi; ++i) {

    mat W = reshape(D * V.row(i).t(), d, d);
    mat v = tau(i) * C.t() * W * C;
    mat w = diagmat(ones(d1) / diagvec(v));
    mat beta3 = reshape(vecbeta3.col(i), K * p, d1).t();
    sum1 = sum1 + beta3.t() * w * beta3;
    sum2 = sum2 + beta3.t() * w * beta1.col(i);

  }

  colvec nu = solve(sum1, sum2);

  mat Fhat = solve(
    (Z.rows(0, T_obs - 2).t() * Z.rows(0, T_obs - 2)),
    Z.rows(0, T_obs - 2).t() * F_mat.rows(1, T_obs - 1)
  ).t();

  mat lambda(T_obs - 1, K);
  mat nuT(T_obs - 1, K);
  colvec Fhat_vec = vectorise(Fhat.t());
  mat LAMBDA(K, p);

  if ( dataset.compare("individual") != 0) {

    colvec LAMBDAv = nu + Fhat_vec;
    LAMBDA = reshape(LAMBDAv, p, K).t();
    lambda = Z.rows(0,T_obs - 2) * LAMBDA.t();
    mat NU = reshape(nu, p, K).t();
    nuT = Z.rows(0,T_obs - 2) * NU.t();

  } else {

    colvec nuB = UnbiasedEst_cond_cpp(nu, vecbeta3, V, K, p, q, nchi, T_obs, tau, d1, d2, d);

    colvec LAMBDAv = nuB + Fhat_vec;
    LAMBDA = (reshape(LAMBDAv, p, K)).t();
    lambda = Z.rows(0,T_obs - 2) * LAMBDA.t();
    mat NU = (reshape(nuB, p, K).t());
    nuT = Z.rows(0,T_obs - 2) * NU.t();

  }

  mat beta = join_cols(beta1, beta2);
  mat RiskPremia = lambda * 12 * 100 / 10;
  // nuT = nuT * 12 * 100 / 10;
  mat a_timevar(T_obs - 1, nchi); a_timevar.zeros();
  mat Yhat(T_obs - 1, nchi); Yhat.zeros();

  for (int i=0; i<nchi; ++i) {
    mat bi_timevar = b_timevar.slice(i);
    for (int t_obs=0; t_obs<(T_obs - 1); ++t_obs) {
      a_timevar(t_obs, i) = conv_to<double>::from(bi_timevar.row(t_obs) * nuT.row(t_obs).t());
      Yhat(t_obs, i) = conv_to<double>::from(bi_timevar.row(t_obs) * RiskPremia.row(t_obs).t());
    }
  }


  return List::create(
    Named("beta") = beta
    , Named("a_timevar") = a_timevar
    , Named("b_timevar") = b_timevar
    , Named("RiskPremia") = RiskPremia
    , Named("Yhat") = Yhat
  );
}
