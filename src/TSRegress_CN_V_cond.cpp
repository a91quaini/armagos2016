#include <RcppArmadillo.h>
#include "CondReg.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List TSRegress_CN_V_cond_cpp(mat R, colvec Rf, mat F_mat, mat Z, cube Zi) {

  std::ostream nullstream(0);
  arma::set_stream_err2(nullstream);

  int n = R.n_cols;
  int T_obs = R.n_rows ;
  int K = F_mat.n_cols;
  int p = Z.n_cols;
  int q = Zi.n_cols;
  int d1 = p * (p + 1) / 2 + p * q;
  int d2 = K * (p + q);
  int d = d1 + d2;

  cube X = CondReg_cpp(n, T_obs, d1, d2, Z, Zi, F_mat);

  mat I(size(R)); I.zeros(); I(find(R)).ones();
  colvec Ti = sum(I).t();

  mat ExcR = I % (R - repmat(Rf, 1, n));
  ExcR = ExcR.rows(1, T_obs - 1); I = I.rows(1, T_obs - 1); Ti = sum(I).t();

  mat beta1(d1, n); beta1.zeros();
  mat beta2(d2, n); beta2.zeros();
  colvec tau = (T_obs - 1) / Ti;
  colvec CN(n); CN.zeros();
  mat V(n, d * (d + 1) / 2); V.zeros();
  mat res(T_obs - 1, n); res.zeros();
  colvec RSS(n);
  cube QQxi(d, d, n);
  mat beta(d, n);

  for (int i=0; i<n; ++i) {

    mat Ii = repmat(I.col(i), 1, d);
    mat Xi = Ii % X.slice(i).t();
    mat Qxi = (Xi.t() * Xi) / Ti(i);
    QQxi.slice(i) = Qxi;

    mat beta = solve(Qxi, ((Xi.t() * ExcR.col(i)) / Ti(i)));
    beta1.col(i) = beta.rows(0, d1 - 1).col(0);
    beta2.col(i) = beta.rows(d1, beta.n_rows - 1).col(0);
    res.col(i) = ExcR.col(i) - Xi * beta;

    RSS(i) = conv_to<double>::from(sum(pow(res.col(i), 2)));

    mat Xtildei(T_obs - 1, d); Xtildei.zeros();
    for (int s = 0; s < d; ++s) {
      Xtildei.col(s) = Xi.col(s) / norm(Xi.col(s), 2);
    }

    mat Qxtildei = (Xtildei.t() * Xtildei) / Ti(i);
    // colvec l = eig_sym(Qxtildei);
    //
    // if (l.max() > 0) {
    //   CN(i) = sqrt(l.max() / l.min());
    // } else {
    //   CN(i) = 100;
    // }

    double l = rcond(Qxtildei);

    if (l > 1 / pow(15, 2)) {
      CN(i) = l;
    } else {
      CN(i) = 100;
    }

    mat ri = repmat(res.col(i), 1, d) % Xi;
    mat Sii = ri.t() * ri / Ti(i);
    mat Vmat = solve(Qxi.t(), solve(Qxi, Sii).t()).t();

    for(int j=0; j<(d * (d + 1) / 2); ++j) {
      for(int k=0; k<(K + 1); ++k) {
        for(int r=0; r<(K + 1); ++r) {
          if(k<=r) {
            V(i, j) = Vmat(k, r);
          }
        }
      }
    }

  }

  int chi1 = 15;
  int chi2 = (T_obs - 1) / 60;
  uvec ind = find(CN < chi1);

  beta1 = beta1.cols(ind); beta2 = beta2.cols(ind); V = V.rows(ind); I = I.cols(ind);
  res = res.cols(ind); tau = tau(ind); Ti = Ti.rows(ind); RSS = RSS(ind); ExcR = ExcR.cols(ind);
  cube X_ind(X.n_rows, X.n_cols, ind.n_elem);
  cube QQxi_ind(d, d, ind.n_elem);

  for (int i=0; i<ind.n_elem; ++i) {

    int j = ind(i);
    X_ind.slice(i) = X.slice(j);
    QQxi_ind.slice(i) = QQxi.slice(j);

  }

  uvec ind1 = find(tau <= chi2);
  beta1 = beta1.cols(ind1); beta2 = beta2.cols(ind1); V = V.rows(ind1); I = I.cols(ind1);
  res = res.cols(ind1); tau = tau(ind1); Ti = Ti.rows(ind1);
  cube X_ind1(X.n_rows, X.n_cols, ind1.n_elem); RSS = RSS(ind1); ExcR = ExcR.cols(ind1);

  int nchi = Ti.n_elem;
  colvec sigmai_hat = sqrt(RSS / (Ti - d));
  mat beta1_sderr(ind1.n_elem, d1);
  mat beta2_sderr(ind1.n_elem, d2);

  for (int i=0; i<ind1.n_elem; ++i) {

    int j = ind1(i);
    X_ind1.slice(i) = X_ind.slice(j);
    mat Qxi_ind1 = QQxi_ind.slice(j);
    colvec Qxi_inv = diagvec(inv_sympd(Qxi_ind1));
    Qxi_inv = sqrt(Qxi_inv);

    beta1_sderr.row(i) = sigmai_hat(i) * Qxi_inv(span(0, d1 - 1)).t();
    beta2_sderr.row(i) = sigmai_hat(i) * Qxi_inv(span(d1, d - 1)).t();

  }

  vec ind_last = linspace(0, n - 1, n);
  X = X_ind1;
  mat beta_sderr = join_rows(beta1_sderr, beta2_sderr);

  return List::create(
    Named("n") = n
    , Named("nchi") = nchi
    , Named("ExcR") = ExcR
    , Named("T_obs") = T_obs
    , Named("K") = K
    , Named("p") = p
    , Named("q") = q
    , Named("d1") = d1
    , Named("d2") = d2
    , Named("d") = d
    , Named("beta1") = beta1
    , Named("beta2") = beta2
    , Named("beta_sderr") = beta2_sderr
    , Named("V") = V
    , Named("tau") = tau
  );
}

