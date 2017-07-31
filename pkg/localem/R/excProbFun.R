#' @title Computes the exceedance probabilities
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds and expected counts from the rasterization of the spatial polygons of population data, and then, computes the exceedance probablities with the same bandwidth as the risk estimation on the raster cells. 
#' 
#' @param cases Spatial polygons, data frame or vector of case data
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param bw Bandwidth specifying which smoothing matrix in \code{lemObjects} to use
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param iterations Convergence tolerance, number of iterations, and use of gpuR package for running local-EM recursions
#' @param verbose Verbose output
#' @param path Folder for storing rasters
#' 
#' 
#' @details After using the \code{excProb} function, the exceedance probabilities are computed on the raster cells of the fine polygons. 
#' 
#' @return The \code{excProb} function returns a raster brick of exceedance probabilities of input risk thresholds. 
#' 
#' @examples 
#' \dontrun{ 
#' data('kentuckyCounty')
#' data('kentuckyTract')
#' 
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, 
#'								polyFine = kentuckyTract, 
#'  	                        cellsCoarse = 6, 
#'                              cellsFine = 100, 
#'                              bw = c(10, 15) * 1000, 
#'                              ncores = 2, 
#'                            	path = 'example', 
#'                              verbose = TRUE)
#'
#'
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'                                 ncores = 2, 
#'                                 verbose = TRUE)
#'
#' lemExcProb = excProb(cases = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#'                    bw = 10000, 
#'					  threshold = c(1, 1.25, 1.5), 
#'					  Nboot = 100, 
#'                    ncores = 2, 
#'                    path = 'example', 
#'                    verbose = TRUE)
#'
#' pCol = mapmisc::colourScale(lemExcProb, 
#'                            breaks = c(0, 0.2, 0.8, 0.95, 1), style = 'fixed', 
#'                            col = c('green', 'yellow', 'orange', 'red'))								
#' plot(lemExcProb$excProb, 
#'     main = 'Exceedance Probabilities, t = 1', 
#'     col = pCol$col, breaks = pCol$breaks, 
#'     legend = TRUE)
#'}
#'
#' @export
excProb = function(
  cases, 
  lemObjects, 
  bw, 
  threshold = 1, 
  Nboot = 100, 
  ncores = 1, 
  iterations = list(tol = 1e-5, maxIter = 1000, gpu=FALSE), 
  verbose=FALSE, 
  path = getwd()

){
  
  dir.create(path, showWarnings=FALSE, recursive=TRUE)

  # warning messages
  if(missing(lemObjects)) {
	stop("smoothing matrix not supplied")
  }
  if(length(bw) > 1) {
    stop("bw must be length 1")
  }
  bwString = paste('bw', bw, sep = '')
  
  # bootstrap cases
	offsetDf = aggregate(x = values(lemObjects$offset$offset) * prod(res(lemObjects$offset)), 
					by = list(ID = values(lemObjects$rasterFine)), 
					FUN = sum)
	colnames(offsetDf) = c('ID','offset')
	
	rasterDf = levels(lemObjects$rasterFine)[[1]]
	offsetDf = merge(offsetDf, rasterDf, by = 'ID')
	
	
	offset = drop(tapply(offsetDf$offset, offsetDf$idCoarse, sum))
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
	estRisk = riskEst(cases = cases, 
					lemObjects = lemObjects, 
					bw = bw,  
					ncores = 1, 
					path = 'example', 
					verbose = verbose)

	#risk and exceedance probabilities
	estRiskBoot = riskEst(cases = bootCounts, 
					lemObjects = lemObjects, 
					bw = bw,  
					ncores = ncores, 
					path = 'example', 
					verbose = verbose)
	
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
