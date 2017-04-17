#' This function estimates time-invariant risk premia.
#'
#' @param R = matrix containing a possibly unbalanced panel of asset returns.
#' @param Rf = vector containing risk free returns.
#' @param F_mat = matrix containing factors (without constant).
#' @param dataset = character vector indicating if ExcR are "portfolio" ExcR or "individual" stocks ExcR.
#'
#' @return out = list of integer n, integer nchi(trimmed cross-sectional dimension),
#' vector a (intercept), matrix b (factor loadings), 3d-array CI (confidence intervals),
#' vector nu (annualized nu), vector RiskPremia (annualized risk premia).
#'
#' @export
#'

EstTimeInvariantMod <- function(R, Rf, F_mat, dataset = "individual") {

  #   First pass (see Section 3.2)
  #   Time-series regression
  TSRegress_CN_V_out <- TSRegress_CN_V_cpp(R, Rf, F_mat)
  list2env(TSRegress_CN_V_out, environment())

  # Second pass (see eq. (14))
  # Cross-sectional regression
  CR_reg_out <- CR_reg_cpp(a, b, V, Vxi, tau, F_mat, n, nchi, T_obs, K, dataset)
  list2env(CR_reg_out, environment())

  out <- list(
    n = n
    , nchi = nchi
    , ExcR = ExcR
    , a = a
    , a_sderr = a_sderr
    , b = b
    , b_sderr = b_sderr
    , RiskPremia = RiskPremia
    , Yhat = Yhat
    )

  return(out)

}
