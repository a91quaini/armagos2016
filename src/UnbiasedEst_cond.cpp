#include <RcppArmadillo.h>
#include "duplication.h"
#include "Cmatrix.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
colvec UnbiasedEst_cond_cpp(colvec nu, mat vecbeta3, mat V, int K, int p, int q, int n,
  int T_obs, colvec tau, int d1, int d2, int d) {

  List Cmatrix_out = Cmatrix_cpp(nu, p, q, K, d1, d2);
  mat C = Cmatrix_out["C"];
  mat Ja = Cmatrix_out["Ja"];
  mat Jb = kron(vectorise(eye(d1, d1)).t(),
    eye(K * p, K * p)) * kron(eye(d1, d1), Ja);
  mat E2 = join_cols(zeros(d1, d2), eye(d2, d2));
  mat w(d1, n); w.zeros();

  mat D = duplication_cpp(d);
  colvec Bs1 = zeros(d1 * d2);
  mat Qb3W = zeros(K * p, K * p);

  for (int i=0; i<n; ++i) {

    colvec vecV = D * V.row(i).t();
    mat W = reshape(vecV, d, d);
    mat v = tau(i) * C.t() * W * C;
    w.col(i) = diagvec(v);
    mat wi = diagmat(ones(d1) / w.col(i));
    Bs1 = Bs1 + tau(i) * vectorise(E2.t() * W * C * wi);
    mat beta3 = (reshape(vecbeta3.col(i), K * p, d1)).t();
    Qb3W = Qb3W + beta3.t() * wi * beta3;

  }

  Qb3W = Qb3W / n;
  mat Bnu = solve(Qb3W, Jb) * Bs1 / n;
  mat nuB = nu - Bnu / T_obs;

  return nuB;
}

