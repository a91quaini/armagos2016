#ifndef __UNBIASEDEST_H_INCLUDED__
#define __UNBIASEDEST_H_INCLUDED__

arma::colvec UnbiasedEst_cpp(arma::colvec nu, arma::mat b, arma::mat V,
  arma::mat Vxi, arma::colvec tau, int n, int T_obs, int K);

#endif
