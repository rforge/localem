#' @title Computes the exceedance probabilities on the raster of fine polygons
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds and expected counts of the fine raster, and then, computes the exceedance probablities with the same bandwidth as the risk estimation on the cells of the fine raster. 
#' 
#' @param lemEst Estimated risk intensity surface
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param ncores Number of cores/threads for parallel processing
#' @param tol Tolerance for convergence
#' @param maxIter Maximum number of iterations for convergence
#' @param filename Passed to writeRaster
#' 
#' 
#' @details After using the \code{excProb} function, the exceedance probabilities are computed on the raster cells of the fine polygons. 
#' 
#' @return The \code{excProb} function returns a raster brick of exceedance probabilities of input risk thresholds. 
#' 
#' @examples 
#' \dontrun{ 
#' data(kentuckyCounty)
#' data(kentuckyTract)
#' 
#' ncores = 1 + (.Platform$OS.type == 'unix')
#' 
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, 
#'                            polyFine = kentuckyTract, 
#'                            cellsCoarse = 6, 
#'                            cellsFine = 100, 
#'                            bw = c(10, 12, 15, 17, 20, 25) * 1000, 
#'                            ncores = ncores, 
#'                            idFile = 'id.grd', 
#'                            offsetFile = 'offset.grd', 
#'                            verbose = TRUE)
#'
#'
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'                                ncores = ncores, 
#'                                verbose = TRUE)
#'
#' lemCv = lemXv(x = kentuckyCounty, 
#'              lemObjects = lemSmoothMat, 
#'              Nxv = 5, 
#'              ncores = ncores, 
#'              verbose = TRUE)
#' bestBw = lemCv$bw[which.min(lemCv$cv)]
#'
#' lemRisk = riskEst(x = kentuckyCounty, 
#'                  lemObjects = lemSmoothMat, 
#'                  bw = bestBw, 
#'                  ncores = ncores)
#'                    
#' lemExcProb = excProb(lemEst = lemRisk, 
#'                     lemObjects = lemSmoothMat, 
#'                     threshold = c(1, 1.25, 1.5), 
#'                     Nboot = 200, 
#'                     ncores = ncores)
#'
#' pCol = mapmisc::colourScale(lemExcProb, 
#'                            breaks = c(0, 0.2, 0.8, 0.95, 1), style = 'fixed', 
#'                            col = c('green', 'yellow', 'orange', 'red'))								
#' plot(lemExcProb[[1]], 
#'     main = 'Exceedance Probabilities, t = 1', 
#'     col = pCol$col, breaks = pCol$breaks, 
#'     legend = TRUE)
#'}
#'
#' @export
excProb = function(
  lemEst, 
  lemObjects, 
  threshold = 1, 
  Nboot = 100, 
  ncores = 1, 
  tol = 1e-6, 
  maxIter = 2000, 
  filename = ''
){
  
  # #observed risk surface
  # theLemRisk = estimate
  
 # result = theLemRisk
  # names(result) = paste("risk.", bw, sep = "")

  
  # bootstrap cases
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

	 # bandwidth of the risk estimate
  bw = as.numeric(
    gsub("^risk.", "", names(lemEst))
  )

	#risk and exceedance probabilities
	estRiskBoot = riskEst(
			x=bootCounts, lemObjects=lemObjects,
			bw=bw, tol=tol, maxIter=maxIter, ncores=ncores
	)
	
	values(estRiskBoot)  = values(estRiskBoot) < rep(values(lemEst), nlayers(estRiskBoot))		

	tIndex = gsub(
			"^risk.threshold.|.sim.[[:digit:]]+$", "", 
					names(estRiskBoot))
	
	tIndex = factor(tIndex)		
	theExcProb = stackApply(
			estRiskBoot,
			as.numeric(tIndex),
			mean)		
	
	
	names(theExcProb) = paste("threshold.", levels(tIndex), sep='')		

	result = do.call(brick, 
			c(unstack(theExcProb), list(filename=filename)))		
	
	return(result)
}
