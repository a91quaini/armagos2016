#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

cube CondReg_cpp(int n, int T_obs, int d1, int d2, mat Z, cube Zi, mat F_mat) {

  mat tri_l, tri_u;
  mat Xt, Xit;
  cube X(d1 + d2, T_obs - 1, n);
  X.fill(0);

  for (int i = 0; i < n; ++i) {

    mat x1(d1, T_obs); x1.zeros();
    mat x2(d2, T_obs); x2.zeros();

    for (int t_obs = 1; t_obs < T_obs; ++t_obs) {

      tri_l = Z.row(t_obs - 1).t() * Z.row(t_obs - 1);
      tri_u = Z.row(t_obs - 1).t() * Z.row(t_obs - 1);
      tri_l = trimatl(tri_l); tri_l.diag().zeros();
      tri_u = trimatu(tri_u); tri_u.diag().zeros();

      Xt = Z.row(t_obs - 1).t() * Z.row(t_obs - 1) + tri_l + tri_u;
      mat vech_Xt = trimatl(Xt); vech_Xt = nonzeros(vech_Xt);
      mat kron_F_Z = kron(F_mat.row(t_obs), Z.row(t_obs - 1));
      Xit = Zi.slice(i).row(t_obs - 1).t() * Z.row(t_obs - 1);
      mat vec_Xti = Xit.t();
      mat kron_F_Zi = kron(F_mat.row(t_obs), Zi.slice(i).row(t_obs - 1));

      x1.col(t_obs) = join_cols(vech_Xt, vec_Xti);
      x2.col(t_obs) = join_cols(kron_F_Z.t(), kron_F_Zi.t());

    }

    x1 = x1.cols(1, T_obs - 1);
    x2 = x2.cols(1, T_obs - 1);
    X.slice(i) = join_cols(x1, x2);

  }

  return X;
}
