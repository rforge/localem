#' @title Computes the exceedance probabilities
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds and expected counts from the rasterization of the spatial polygons of population data, and then, computes the exceedance probablities with the same bandwidth as the risk estimation on the raster cells. 
#' 
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition, smoothed offsets and risk estimation
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param ncores Number of cores/threads for parallel processing
#' @param iterations Convergence tolerance, number of iterations, and use of gpuR package for running local-EM recursions
#' @param verbose Verbose output
#' @param path Folder for storing rasters
#' @param filename Filename (must have .grd extension) of the exceedance probabilities
#' 
#' @details After using the \code{excProb} function, the exceedance probabilities are computed on a fine resolution based on the rasterization of the spatial polygons of population data.
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
#' lemExcProb = excProb(lemObjects = lemSmoothMat, 
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
  lemObjects, 
  threshold = 1, 
  Nboot = 100, 
  ncores = 1, 
  iterations = list(tol = 1e-5, maxIter = 1000, gpu=FALSE), 
  verbose=FALSE, 
  path = getwd(), 
  filename = file.path(path, "prob.grd")
){
  
	dir.create(path, showWarnings=FALSE, recursive=TRUE)

	# warning messages
	if(missing(lemObjects)) {
		stop("smoothing matrix and rasters not supplied")
	}

	# bandwidth of interest
	bw = as.numeric(gsub('^bw', '', lemObjects$bw))
	
	# risk estimate of interest
	theEstRisk = lemObjects$estimate

	if(verbose) {
		cat("generating bootstrap cases for input thresholds\n")
	}
	
	# offsets of spatial polygons of case data based on population data
 	idCoarse = 1:length(lemObjects$smoothingMatrix$polyCoarse)
		  
	offsetRaster = raster::stack(lemObjects$smoothingMatrix$offset$offset, 
									raster::deratify(lemObjects$smoothingMatrix$rasterFine))
	offsetDf = aggregate(x = values(offsetRaster$offset) * prod(res(offsetRaster)), 
							by = list(idCoarse = values(offsetRaster$idCoarse)), 
							FUN = sum)
	colnames(offsetDf) = c('idCoarse','offset')
	offsetDf = merge(data.frame(idCoarse = idCoarse), offsetDf, by = 'idCoarse', all = TRUE)
	offsetDf$offset[is.na(offsetDf$offset)] = 0
	
	# bootstrap cases
	offsetT = outer(offsetDf$offset, threshold)

	bootCountsDf = matrix(data = stats::rpois(length(offsetT) * Nboot, rep(offsetT, Nboot)), 
							nrow = nrow(offsetDf), 
							ncol = length(threshold) * Nboot)
	dimnames(bootCountsDf) = list(rownames(offsetDf), 
									paste('count', rep(1:Nboot, rep(length(threshold), Nboot)), 
											'_threshold', rep(threshold, Nboot), 
											sep = ''))

	# estimate risk from bootstrap cases
	if(verbose) {
		cat("running local-EM estimation for bootstrap cases with original bw\n")
	}

	bootLemRisk = riskEst(cases = bootCountsDf, 
						lemObjects = lemObjects$smoothingMatrix, 
						bw = bw,  
						ncores = ncores, 
						iterations = iterations, 
						path = path, 
						filename = file.path(path, "riskBoot.grd"), 
						verbose = FALSE)
	bootEstRisk = bootLemRisk$estimate
	
	# exceedance probabilities
	if(verbose) {
		cat("computing exceedance probabilities with input thresholds\n")
	}
	
	indexT = gsub('^bw[[:digit:]]+_count[[:digit:]]+_', '', names(bootEstRisk))
	
	bootEstRisk = raster::overlay(x = bootEstRisk, y = theEstRisk, 
								fun = function(x,y) return(x < y), 
								filename = file.path(path, "probBoot.grd"), 
								overwrite = TRUE)

	theExcProb = raster::stackApply(bootEstRisk, 
								indices = indexT, 
								fun = mean, 
								filename = filename, 
								overwrite = TRUE)
	names(theExcProb) = gsub('index', paste('bw', bw, sep = ''), names(theExcProb))
	
	result = list(
				estimate = theEstRisk, 
				excProb = theExcProb
			)
	
	if(verbose) {
		cat("done\n")
	}	

	return(result)
}
