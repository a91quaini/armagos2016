#' This function is project specific and it transfroms specific tidy data into wide data for reproducing GOS 2016.
#'
#' @param asset_dataf tidy dataframe containing assets information.
#' it must have (at least) the following columns: date, value (total returns),
#' hash_key_entity (uniquely identifying the assets).
#' @param Rf NULL for this project. Ideally it would be a dataframe containing risk free information.
#' @param factor_dataf tidy dataframe containig factors information.
#' it must have (at least) the following columns: date, value, hash_key_entity (uniquely identifying the factors).
#' @param instrument_dataf tidy dataframe containing information on the instrumental variables
#' for time varying coefficients. If model = 0, is set to NULL. If model = 1,
#' it must have the following columns: date, value, hash_key_entity (uniquely identifying the instruments).
#'
#' @return wide data ready for estimation.
#'
#' @export
#'

GOS_tidy_input <- function(asset_dataf = lhs_data, Rf = NULL, factor_dataf = rhs_data, model = 0
  , instument_dataf = tv_coeff_data, dataset = "individual") {

  # asset_dataf <- asset_dataf %>%
  #   dplyr::select(
  #     date
  #     , hash_key_entity
  #     , value
  #   ) %>%
  #   reshape2::dcast(date ~ hash_key_entity)
  #
  # R <- asset_dataf[, -1] %>% as.matrix()
  # Rf <- rep(0, nrow(asset_dataf))
  #
  # factor_dataf <- factor_dataf %>%
  #   dplyr::select(
  #     date
  #     , hash_key_entity
  #     , value
  #   ) %>%
  #   reshape2::dcast(date ~ hash_key_entity)
  #
  # F_mat <- factor_dataf[, -1] %>% as.matrix()
  #
  # out <- list(
  #   R = R
  #   , Rf = Rf
  #   , F_mat = F_mat
  #   )
  #
  # if (model == 1) {
  #
  #   instrument_dataf <- instrument_dataf %>%
  #     dplyr::select(
  #       date
  #       , hash_key_entity
  #       , value
  #     ) %>%
  #     reshape2::dcast(date ~ hash_key_entity)
  #
  # }


}

