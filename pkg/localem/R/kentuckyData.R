#' @title Lung cancer cases in Kentucky
#'
#' @description The dataset contains lung cancer cases for the counties of Kentucky.
#'
#' @name kentuckyCounty
#'
#' @docType data
#'
#' @usage data(kentuckyCounty)
#'
#' @format A SpatialPolygonsDataFrame of Kentucky county boundaries, populations and cancer cases.
#'
#' @details The lung cancer case data are obtained from www.cancer-rates.info and were simulated from a log Gaussian Cox process for deliberately unspecified years.
#'
#' @keywords datasets
#'
#' @examples
#' data(kentuckyCounty)
#' kentuckyCounty
#' spplot(kentuckyCounty, 'count')
#'
NULL

#' @title Population data in Kentucky
#'
#' @description The dataset contains population data for the census tracts of Kentucky.
#'
#' @name kentuckyTract
#'
#' @docType data
#'
#' @usage data(kentuckyTract)
#'
#' @format A SpatialPolygonsDataFrame of Kentucky census tract boundaries and populations.
#'
#' @details The lung cancer case data are obtained from www.cancer-rates.info and were simulated from a log Gaussian Cox process for deliberately unspecified years.
#'
#' @keywords datasets
#'
#' @examples
#' data(kentuckyTract)
#' spplot(kentuckyTract, 'logOffset')
#'
NULL
