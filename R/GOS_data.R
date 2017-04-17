#' GOS 2016 data preparation
#'
#' @param dataset integer indicating the dataset to analyse.
#' 1 = FF25, 2 = Indu44, 3 = CRSPCMST_ret. Each dataset contains:
#' R=returns, I=indicator for unbalanced characteristic, BM=book-to-mkt
#' Ti= time-series observations for each asset, T=TS size, n=CS size
#' @param K integer declaring the number of factors (mkt return, small minus big, high minus low, momentum).
#' default is 4.
#' @param Model integer indicating the model. 1 (default) = time variant. 0 = time invariant.
#' @param p number of common instruments Z (1, 2 or 3 [default]).
#' @param q number of asset specific instruments Zi (0 or 1 [default]).
#'
#' @export
#'

GOS_data <- function(dataset, K, model, p, q) {

  ## data
  if (dataset == 1) {
    dataobj = FF25
  } else if ( dataset == 2) {
    dataobj = Indu44
  } else if (dataset == 3) {
    dataobj = CRSPCMST_ret
  }

  R <- dataobj$R
  # Factors
  F_mat <- factors
  F_mat <- as.matrix(F_mat[, 1:K] * 10)    # rescale factors
  # RiskFree
  Rf <- RiskFree
  dataset <- ifelse((dataset == 1 || dataset == 2), "portfolio", "individual")

  out <- list(
    R = R
    , F_mat = F_mat
    , Rf = Rf
    , dataset = dataset
  )

  if (model == 1) {

    n <- dataobj$n
    T_obs <- dataobj$T_obs
    Ti <- dataobj$Ti
    I <- dataobj$I

    DefSpread <- Instrum$DefSpread
    TermSpread <- Instrum$TermSpread

    ds <- (DefSpread - mean(DefSpread)) / sd(DefSpread)
    ts <- (TermSpread - mean(TermSpread)) / sd(TermSpread)
    if (p == 2) {
      Z <- matrix(c(rep(1, T_obs), ds, ts), nrow = T_obs, ncol = p + 1)
    } else {
      Z <- matrix(c(rep(1, T_obs), ds), nrow = T_obs, ncol = p + 1)
    }

    if (q == 0) Z <- as.matrix(Z[, -1])

    if (q != 0 ) {

      BM <- dataobj$BM
      # asset specific instruments
      Zi <- array(0, dim = c(T_obs, q, n))
      for(i in 1:n) {
        # Standardization of BM (book-to-mkt)
        m <- sum(BM[, i]) / Ti[i]
        repmat <- do.call(rbind, replicate(T_obs, m, simplify = FALSE))
        s <- sum((BM[, i] - I[, i] * repmat)^2 / Ti[i])
        s <- sqrt(s)
        z1 <- I[, i] * ((BM[, i] - m)/s)
        Zi[, , i] <- z1
      }

    } else {
      Zi <- array(1, dim = c(T_obs, 1, n))
    }

    out[["Zi"]] <- Zi
    out[["Z"]] <- Z

  }

  return(out)

}
