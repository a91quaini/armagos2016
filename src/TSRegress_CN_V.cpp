#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List TSRegress_CN_V_cpp(mat R, colvec Rf, mat F_mat) {

  std::ostream nullstream(0);
  arma::set_stream_err2(nullstream);

  int n = R.n_cols;
  int T_obs = R.n_rows;
  int K = F_mat.n_cols;

  mat I(size(R)); I.zeros(); I(find(R)).ones();
  colvec Ti = sum(I).t();

  mat ExcR = I % (R - repmat(Rf, 1, n));
  mat X(T_obs, K + 1); X.col(0).ones(); X.cols(1, K) = F_mat;
  mat Qx(K + 1, K + 1); Qx = X.t() * X / T_obs;

  colvec tau = T_obs / Ti;
  colvec a(n); a.zeros();
  mat b(n, K); b.zeros();
  colvec CN(n); CN.zeros();
  int Kk = (K + 1) * (K + 1 + 1) / 2;
  mat V(n, Kk); V.zeros();
  mat Vxi(n, Kk); Vxi.zeros();
  mat res(T_obs, n); res.zeros();
  colvec RSS(n);
  cube QQxi(K + 1, K + 1, n);
  mat Vmat(K + 1, K + 1);

  for (int i=0; i<n; ++i) {

    mat Ii = repmat(I.col(i), 1, K + 1);
    mat Xi = Ii % X;
    mat Qxi = (Xi.t() * Xi) / Ti(i);
    QQxi.slice(i) = Qxi;

    colvec beta = solve(Qxi, ((Xi.t() * ExcR.col(i)) / Ti(i)));
    a(i) = beta(0);
    b.row(i) = beta(span(1, K)).t();
    res.col(i) = ExcR.col(i) - Xi * beta;

    RSS(i) = conv_to<double>::from(sum(pow(res.col(i), 2)));

    mat Xtildei(T_obs, K + 1); Xtildei.zeros();
    for (int s = 0; s < K + 1; ++s) {
      Xtildei.col(s) = Xi.col(s) / norm(Xi.col(s), 2);
    }
    mat Qxtildei = (Xtildei.t() * Xtildei) / Ti(i);
    // colvec l = eig_sym(Qxtildei);

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

    mat ri = repmat(res.col(i), 1, K + 1) % Xi;
    mat Sii = ri.t() * ri / Ti(i);
    Vmat = solve(Qx.t(), solve(Qx, Sii).t()).t();
    mat Vmati = solve(Qxi.t(), solve(Qxi, Sii).t()).t();

    for(int j=0; j<Kk; ++j) {
      for(int k=0; k<(K + 1); ++k) {
        for(int r=0; r<(K + 1); ++r) {
          if(k<=r) {
            V(i, j) = Vmat(k, r);
            Vxi(i, j) = Vmati(k, r);
          }
        }
      }
    }

  }

  int chi1 = 15;
  int chi2 = T_obs / 12;
  uvec ind = find(CN < chi1);

  a = a(ind); b = b.rows(ind); V = V.rows(ind); Vxi = Vxi.rows(ind); I = I.cols(ind);
  res = res.cols(ind); tau = tau(ind); Ti = Ti.rows(ind); RSS = RSS(ind); ExcR = ExcR.cols(ind);
  cube QQxi_ind(K + 1, K + 1, ind.n_elem);

  for (int i=0; i<ind.n_elem; ++i) {

    int j = ind(i);
    QQxi_ind.slice(i) = QQxi.slice(j);

  }

  uvec ind1 = find(tau <= chi2);
  a = a(ind1); b = b.rows(ind1); V = V.rows(ind1); Vxi = Vxi.rows(ind1); I = I.cols(ind1);
  res = res.cols(ind1); tau = tau(ind1); Ti = Ti.rows(ind1); RSS = RSS(ind1); ExcR = ExcR.cols(ind1);

  int nchi = Ti.n_elem;
  colvec sigmai_hat = sqrt(RSS / (Ti - (K + 1)));
  colvec a_sderr(ind1.n_elem);
  mat b_sderr(ind1.n_elem, K);

  for (int i=0; i<ind1.n_elem; ++i) {

    int j = ind1(i);
    mat Qxi_ind1 = QQxi_ind.slice(j);
    colvec Qxi_inv = diagvec(inv_sympd(Qxi_ind1));
    Qxi_inv = sqrt(Qxi_inv);

    a_sderr(i) = sigmai_hat(i) * Qxi_inv(0);
    b_sderr.row(i) = sigmai_hat(i) * Qxi_inv(span(1, K)).t();

  }

  vec ind_last = linspace(0, n - 1, n);
  ind_last = ind_last(ind); ind_last = ind_last(ind1);

  return List::create(
    Named("n") = n
    , Named("ind_last") = ind_last
    , Named("nchi") = nchi
    , Named("ExcR") = ExcR
    , Named("T_obs") = T_obs
    , Named("K") = K
    , Named("a") = a
    , Named("a_sderr") = a_sderr
    , Named("b") = b
    , Named("b_sderr") = b_sderr
    , Named("V") = V
    , Named("Vxi") = Vxi
    , Named("tau") = tau
  );
}


