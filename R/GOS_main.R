#' Reproduce GOS 2016
#'
#' @param dataset integer indicating the dataset to analyse.
#' 1 = FF25, 2 = Indu44, 3 [default] = CRSPCMST_ret, 4 = tidy data. Each dataset contains:
#' R=returns, I=indicator for unbalanced characteristic, BM=book-to-mkt
#' Ti= time-series observations for each asset, T=TS size, n=CS size
#' @param K integer declaring the number of factors (mkt return, small minus big, high minus low, momentum).
#' default is 4.
#' @param Model integer indicating the model. 1 [default] = time variant. 0 = time invariant.
#' @param p number of common instruments Z (1, 2 [default]).
#' @param q number of asset specific instruments Zi (0 [default] or 1).
#'
#' @return out = list. If dataset = 1 or 2, out is a list of integer n, integer nchi(trimmed cross-sectional dimension),
#' vector a (intercept), matrix b (factor loadings), 3d-array CI (confidence intervals),
#' vector nu (annualized nu), vector RiskPremia (annualized risk premia).
#' If dataset = 3, out is a list of integer n, integer nchi(trimmed cross-sectional dimension),
#' matrix beta (coefficients), 3d-array CI (confidence intervals),
#' matrix Table (annualized vec(F'), nu), vector RiskPremia (annualized),
#' vector nuT (annualized time-varying nu).
#'
#' @export
#'
#' @useDynLib armagos2016
#' @importFrom Rcpp sourceCpp
#'

GOS_main <- function(dataset = 1, K = 4, model = 0, p = 1, q = 0) {

  input <- GOS_data(dataset, K, model, p, q)
  list2env(input, environment())

  if (model == 0) {

    #Estimation of no-time varying risk premia, nu, and statistic (see Section 3)
    t0 <- Sys.time()
    twopass <- EstTimeInvariantMod(R, Rf, F_mat, dataset)
    t1 <- Sys.time()
    time_elapsed <- difftime(t1, t0)

  } else {

    #Estimation of time varying risk premia, nu, and statistic (see Section 3)
    t0 <- Sys.time()
    twopass <- EstTimeVariantMod(R, Rf, F_mat, Z, Zi, dataset)
    t1 <- Sys.time()
    time_elapsed <- difftime(t1, t0)
  }

  out <- list(
    twopass
    , time_elapsed
  )
  return(out)

}

