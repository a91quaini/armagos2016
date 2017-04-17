#' This function estimates time-varying risk premia.
#'
#' @param R = matrix containing a possibly unbalanced panel of asset returns.
#' @param Rf = vector containing risk free returns.
#' @param F_mat = matrix containing factors (without constant).
#' @param Z = matrix containing the common instruments.
#' @param Zi = matrix containing the asset-specific instruments.
#' @param dataset = character vector indicating if ExcR are "portfolio" ExcR or "individual" stocks ExcR.
#'
#' @return out = list of integer n, integer nchi(trimmed cross-sectional dimension),
#' matrix beta (coefficients), 3d-array CI (confidence intervals),
#' matrix Table (annualized vec(F'), nu), vector RiskPremia (annualized),
#' vector nuT (annualized time-varying nu).
#'
#' @export
#'

EstTimeVariantMod <- function(R, Rf, F_mat, Z, Zi, dataset = 3) {

  # First pass (see Section 3.2)
  # Estimation of betas
  TSRegress_CN_V_cond_out <- TSRegress_CN_V_cond_cpp(R, Rf, F_mat, Z, Zi)
  list2env(TSRegress_CN_V_cond_out, environment())

  # Second Pass            (see eq. (13))
  # Cross-sectional regression
  CR_reg_cond_out <- CR_reg_cond_cpp(n, nchi, T_obs, K, p, q, d1, d2, d,
    beta1, beta2, V, tau, Z, Zi, F_mat, dataset)
  list2env(CR_reg_cond_out, environment())


  out <- list(
    n = n
    , nchi = nchi
    , ExcR = ExcR
    , beta = beta
    , beta_sderr = beta_sderr
    , a_timevar = a_timevar
    , b_timevar = b_timevar
    , RiskPremia = RiskPremia
    , Yhat = Yhat
    )

}


