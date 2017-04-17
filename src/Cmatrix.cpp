#include <RcppArmadillo.h>
#include "duplication.h"
#include "commutation.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

List Cmatrix_cpp(colvec nu, int p, int q, int K, int d1, int d2) {

  mat D = duplication_cpp(p);
  mat Dplus = solve(D.t() * D, D.t());
  mat W1 = commutation_cpp(p * (p + 1) / 2, p * K);
  mat W2 = commutation_cpp(p, p);
  mat W3 = commutation_cpp(p * q, p * K);
  mat W4 = commutation_cpp(p, q);
  mat E1 = join_cols(eye(d1, d1), zeros(d2, d1));
  mat E2 = join_cols(zeros(d1, d2), eye(d2, d2));
  mat eye_p = eye(p, p);
  mat vec_eye_p = vectorise(eye_p);
  mat eye_K = eye(K, K);

  mat J11 = W1 * (kron(eye_K,
    (kron(eye_p, Dplus) * kron(W2, eye_p) * kron(eye_p, vec_eye_p))));

  mat J22 = W3 * (kron(eye_K,
    (kron(eye_p, W4) * kron(W4, eye_p) * kron(eye(q, q), vec_eye_p))));

  mat Ja = join_cols(
    join_rows(J11, zeros(pow(p, 2) * K * (p + 1) / 2 , q * K)),
    join_rows(zeros(pow(p, 2) * q * K, p * K), J22)
  );

  mat C = (E1.t() - kron(eye(d1, d1), nu.t()) * Ja * E2.t()).t();


  return List::create(
    Named("C") = C
    , Named("Ja") = Ja
  );

}

