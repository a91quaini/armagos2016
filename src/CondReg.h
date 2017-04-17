#ifndef __CONDREG_H_INCLUDED__
#define __CONDREG_H_INCLUDED__

arma::cube CondReg_cpp(int n, int T_obs, int d1, int d2, arma::mat Z, arma::cube Zi, arma::mat F_mat);

#endif
