#' @title Computes the relative risk estimation on the raster of fine polygons
#'
#' @description The \code{lemEst} function first creates the smoothing matrix at the final iteration with the input bandwidth, and then,
#'  computes the estimates of the relative risk on the cells of the fine raster.
#'
#' @param x Spatial polygons of case data
#' @param lemObjects List of arrays for the smoothing matrix and
#'  raster stacks for the partition and smoothed offsets
#' @param bw Numeric value of bandwidth
#' @param ncores Number of cores/threads for parallel processing
#' @param tol tolerance for convergence
#' @param maxIter maximum number of iterations
#' @param verbose verbose output
#' 
#'
#' @details After using the \code{lemEst} function, the raster of risk estimations is done on cells of the raster on the fine polygons.

#' @return The \code{lemEst} function returns a list containing the raster of risk estimations,
#'  array of smoothing matrix at the final iteration, and
#'  input rasters of partitions and smoothed offsets.
#'
#' @examples
#' data(kentuckyCounty)
#'
#' data(kentuckyTract)
#'
#' \dontrun{
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, polyFine = kentuckyTract,
#'                    cellsCoarse = 40, cellsFine = 400,
#'                    bw = c(10, 15, 20, 25) * 1000,
#'                    ncores = 4,
#'                    idFile = 'id.grd', offsetFile = 'offset.grd')
#'
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster,
#'                    ncores = 4)
#'
#' lemCv = lemXv(x = kentuckyCounty,
#'                    lemObjects = lemSmoothMat,
#'                    ncores = 4)
#'
#' lemRisk = lemEst(x = kentuckyCounty,
#'                    lemObjects = lemSmoothMat,
#'                    bw = 15000,
#'                    ncores = 4)
#'
#' plot(lemRisk$risk)
#'}
#'
lemEst = function(x,
  lemObjects,
  bw,
  ncores = 2,
  tol = 1e-6,
  maxIter = 2000,
  verbose = FALSE
) {

  if(length(bw) > 1) {
    stop("Bandwidth must be numeric and length 1")
  }


  if(verbose) {
    cat(date(), "\n")
    cat("generating smoothing matrix at final iteration\n")
  }

  #smoothing matrix for the final iteration
  theFinalMat = smoothingFinalMat(
    lemObjects=lemObjects,
    bw=bw,
    ncores=ncores
  )

  if(verbose) {
    cat(date(), "\n")
    cat("obtaining risk estimation at final iteration\n")
  }

  #risk estimation
  theRisk = riskEst(
    x=x,
    lemObjects=theFinalMat,
    bw=dimnames(theFinalMat$smoothingArray)[[3]],
    tol=tol,
    maxIter=maxIter
  )

  result = c(
      list(risk = theRisk),
      theFinalMat
    )

  return(result)

  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
}
