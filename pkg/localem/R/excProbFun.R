#' @title Computes the exceedance probabilities on the raster of fine polygons
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds, and then, 
#'  computes the exceedance probablities with the same bandwidth as the risk estimation on the cells of the fine raster. 
#' 
#' @param x Spatial polygons of case data
#' @param lemObjects List of arrays for the smoothing matrix and 
#'  raster stacks for the partition and smoothed offsets
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param ncores Number of cores/threads for parallel processing
#' @param tol tolerance for convergence
#' @param maxIter maximum number of iterations
#' @param verbose verbose output
#' 
#' 
#' @details After using the \code{excProb} function, the raster of exceedance probabilities is done on cells of the raster on the fine polygons. 
#' 
#' @return The \code{excProb} function returns a raster stack of the risk estimation and 
#'  exceedance probabilities of specified risk thresholds. 
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
#' lemRisk = riskEst(x = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#'                    bw = 15000, 
#'                    ncores = 4) 
#'                    
#' lemExcProb = excProb(x = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#' 										estimate=lemRisk,
#'                    threshold = c(1, 1.1, 1.25), 
#'                    Nboot = 100, 
#'                    ncores = 4) 
#'                    
#' plot(lemExcProb)
#'}
#'
excProb = function(
		x, 
  	estimate,
    lemObjects,
		bw,
    threshold = 1, 
    Nboot = 100, 
    ncores = 2, 
    tol = 1e-6, 
    maxIter = 2000, 
    verbose = FALSE
){
  
  #observed risk surface
  theLemRisk = estimate
  
  #risk and exceedance probabilities
  result = theLemRisk
  names(result) = paste("risk.", bw, sep = "")
  
  for(inT in 1:length(threshold)) {
    
    if(verbose) {
      cat(date(), "\n")
      cat("obtaining risk estimation of simulated counts for threshold: ", threshold[inT], "\n")
    }
		
    #simulated risk surface
    estRiskBoot = simplify2array(
      	parallel::mclapply(1:Nboot, 
            simLemEst, 
            x = x, 
            lemObjects = lemObjects, 
            threshold = threshold[inT], 
            bw = bw, 
            tol = tol, 
            maxIter = maxIter, 
            mc.cores=ncores
      	))
    
    if(verbose) {
      cat(date(), "\n")
      cat("obtaining exceedance probabilities\n")
    }
    
    #exceedance probabilities
    theExcProb = theLemRisk
    values(theExcProb) = apply(estRiskBoot <= values(theLemRisk), 1, mean)
    names(theExcProb) = paste("threshold.", threshold[inT], sep = "")
    
    result = raster::addLayer(result, theExcProb)
  }

  return(result)
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
}
