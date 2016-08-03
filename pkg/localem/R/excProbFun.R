#' @title Computes the exceedance probabilities on the raster of fine polygons
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds, and then, 
#'  computes the exceedance probablities with the same bandwidth as the risk estimation on the cells of the fine raster. 
#' 
#' @param x Spatial polygons of case data
#' @param estimate Estimated intensity surface
#' @param lemObjects List of arrays for the smoothing matrix and 
#'  raster stacks for the partition and smoothed offsets
#' @param bw bandwidth, specifying which smoothing matrix in \code{lemObjects} to use
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param ncores Number of cores/threads for parallel processing
#' @param tol Tolerance for convergence
#' @param maxIter Maximum number of iterations
#' @param verbose Verbose output
#' @param filename Passed to writeRaster
#' 
#' 
#' @details After using the \code{excProb} function, the raster of exceedance probabilities is done on cells of the raster on the fine polygons. 
#' 
#' @return The \code{excProb} function returns a raster stack of exceedance probabilities of specified risk thresholds. 
#' 
#' @examples 
#' data(kentuckyCounty)
#' 
#' data(kentuckyTract)
#' 
#' \dontrun{
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, polyFine = kentuckyTract,
#'                    cellsCoarse = 6, cellsFine = 100,
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
#' bestBw = lemCv$bw[which.min(lemCv$cv)]
#'
#' lemRisk = riskEst(x = kentuckyCounty,
#'                    lemObjects = lemSmoothMat,
#'                    bw = bestBw)
#'                    
#' lemExcProb = excProb(x = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#' 					estimate=lemRisk,
#'                    threshold = c(1, 1.1, 1.25), 
#'                    Nboot = 100, 
#'                    ncores = 4) 
#'                    
#' plot(lemExcProb)
#'}
#'
#' @export
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
    verbose = FALSE,
	filename=''
){
  
  #observed risk surface
  theLemRisk = estimate
  
  #risk and exceedance probabilities
  result = theLemRisk
  names(result) = paste("risk.", bw, sep = "")
  
	offset = stats::na.omit(as.data.frame(
					stack(
							lemObjects$offset$offset,
							lemObjects$rasterFine$idCoarse
					)
			)
	)
	offset = drop(tapply(
			offset[,'offset'], 
			lemObjects$polyCoarse$id[offset[,'idCoarse']], 
			sum
			)
	)
	offset = offset * prod(res(lemObjects$offset))
	offset = outer(offset, threshold)

	bootCounts = stats::rpois(
			length(offset)* Nboot,
			rep(offset, Nboot)
			)
	bootCounts = matrix(bootCounts, nrow(offset),
			dimnames=list(
					rownames(offset), 
					paste('threshold.',rep(threshold, Nboot),
							'.sim.',rep(1:Nboot, rep(length(threshold), Nboot)), 
									sep='')
			)
	)
	
	if(verbose) {
		cat(date(), "\n")
		cat("obtaining risk estimation of simulated counts\n")
	}
	
	estRiskBoot = riskEst(
			x=bootCounts, lemObjects=lemObjects,
			bw=bw, tol=tol, maxIter=maxIter, ncores=ncores
	)
	
	
	if(verbose) {
		cat(date(), "\n")
		cat("calculating exceedances\n")
	}
			
	values(estRiskBoot)  = values(estRiskBoot) < rep(values(estimate), nlayers(estRiskBoot))		

	tIndex = gsub(
			"^risk.threshold.|.sim.[[:digit:]]+$", "", 
					names(estRiskBoot))
	
	tIndex = factor(tIndex)		
	excProb = stackApply(
			estRiskBoot,
			as.numeric(tIndex),
			mean)		
	
	
	names(excProb) = paste("threshold.", levels(tIndex), sep='')		
	

	if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
	
	excProb
	  
}
