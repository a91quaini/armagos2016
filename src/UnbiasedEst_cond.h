#ifndef __UNBIASEDESTCOND_H_INCLUDED__
#define __UNBIASEDESTCOND_H_INCLUDED__

arma::colvec UnbiasedEst_cond_cpp(arma::colvec nu, arma::mat vecbeta3, arma::mat V, int K, int p, int q, int n,
  int T_obs, arma::colvec tau, int d1, int d2, int d);

#endif
