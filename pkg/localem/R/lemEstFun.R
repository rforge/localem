#' @title Computes the relative risk estimation
#'
#' @description The \code{riskEst} function computes the estimations of the relative risk with high raster resolution. 
#'
#' @param cases Spatial polygons, data frame or vector of case data
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param bw Vector of bandwidths specifying which smoothing matrix in \code{lemObjects} to use
#' @param ncores Number of cores/threads for parallel processing
#' @param iterations List of convergence tolerance, number of iterations, and use of gpuR package for running local-EM recursions
#' @param path Folder for storing rasters
#' @param filename Filename (must have .grd extension) of the risk estimation
#' @param verbose Verbose output
#'
#' @details After using the \code{riskEst} function, the risk estimations are computed on a fine resolution based on the rasterization of the spatial polygons of population data.
#'
#' @return The \code{riskEst} function returns a raster brick of risk estimations for the input bandwidths.
#'
#' @examples 
#' \dontrun{ 
#' # case and population data
#' data('kentuckyCounty')
#' data('kentuckyTract')
#'
#' # parameters
#' ncores = 2
#' cellsCoarse = 8
#' cellsFine = 100
#' bw = c(10, 15, 17.5, 20) * 1000
#' path = 'example'
#' 
#' # rasters of case and population data
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, 
#'								polyFine = kentuckyTract, 
#'								cellsCoarse = cellsCoarse, 
#'								cellsFine = cellsFine, 
#'								bw = bw, 
#'								ncores = ncores, 
#'								path = path, 
#'								idFile = 'lemId.grd', 
#'								offsetFile = 'lemOffsets.grd', 
#'								verbose = TRUE)
#'
#' # smoothing matrix
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'									ncores = ncores, 
#'									path = path, 
#'									filename = 'lemSmoothMat.grd', 
#'									verbose = TRUE)
#'
#' # risk estimation
#' lemRisk = riskEst(cases = kentuckyCounty[,c('id','count')], 
#'						lemObjects = lemSmoothMat, 
#'						bw = bw,  
#'						ncores = ncores, 
#'						path = path, 
#'						filename = 'lemRisk.grd', 
#'						verbose = TRUE)
#'
#' # plot risk
#' rCol = mapmisc::colourScale(lemRisk$riskEst, 
#'							breaks = 5, style = 'quantile', dec = 2)
#'
#' par(mfrow = c(2,2), mar = c(0.5,0.5,3,0.5))
#' for(inBw in 1:length(bw)) {
#'		plot(lemRisk$riskEst[[inBw]], 
#'			main = paste('Risk, bw=', bw[inBw], 'km', sep = ''), 
#'			col = rCol$col, breaks = rCol$breaks, 
#'			axes = FALSE, box = FALSE, legend = FALSE, 
#'			add = FALSE)
#' }
#' mapmisc::legendBreaks('right', rCol)
#'}
#'
#' @export
riskEst = function(
  cases, 
  lemObjects, 
  bw, 
  ncores = 1, 
  iterations = list(tol = 1e-5, maxIter = 1000, gpu = FALSE), 
  path = getwd(), 
  filename,
  verbose = FALSE 
) {
  
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

 	if(missing(filename)) {
		filename = paste(tempfile('lemRisk', path), '.grd', sep = '')
	}
	if(!length(grep('/', filename))) {
		filename = file.path(path, filename)
	}
	if(!(length(grep("\\.gr[id]$", filename)))){
		warning("filename should have .grd extension")
	}
     
  # warning messages
  if(missing(lemObjects)) {
	stop("smoothing matrix not supplied")
  }
  
  if(missing(bw)) {
	cat("using bw from smoothing array\n")
	
	bwString = sort(unique(gsub("xv[[:digit:]]+$", "", lemObjects$bw)))
  } else { 
	  bwString = paste('bw', bw, sep = '')
	  bwMatch = match(bwString, lemObjects$bw)
	  
	  if(any(is.na(bwMatch))) {
	  
		cat("cant match: ", bw[!is.na(bwMatch)], "\n")
		stop("cant match all bw to smoothing array")
	  }
  }

  # checking structure of case data
   if(is.vector(cases)) {
    cases = data.frame(cases=cases)
  }
   if(is.matrix(cases)) {
    cases = as.data.frame(cases)
  }
  if(class(cases) == 'SpatialPolygonsDataFrame') {
    polyCoarse = cases
    cases = cases@data
  } else {
    polyCoarse = NULL
  }
  
  #names of interest for regions in the coarse shapefile
  countcol = grep('^(count|cases)', names(cases), value=TRUE, ignore.case=TRUE)
  # if(!length(countcol)){
    # countcol = grep("^(id|name)", names(cases), invert=TRUE, value=TRUE)[1]
  # }

  if(!length(countcol)) {
	stop("column name of case data must start with 'count' or 'cases'")
  }
  
  
  # parameter specifications (if not supplied)
  if(!length(iterations$tol)) iterations$tol = 1e-6
  if(!length(iterations$maxIter)) iterations$maxIter = 2000
  if(!length(iterations$gpu)) iterations$gpu = FALSE
  
  # initializing parallel features (if specified)
  endCluster = FALSE
  theCluster = NULL
  if(length(grep("cluster", class(ncores))) ) {
    if(verbose) cat("using existing cluster\n")
    theCluster = ncores
	parallel::clusterEvalQ(theCluster, library('raster'))
  } else if(!is.null(ncores)) {
    if(ncores > 1) {
      if(verbose) cat("starting new cluster\n")
      theCluster = parallel::makeCluster(spec=ncores, type='PSOCK', methods=TRUE)
      parallel::setDefaultCluster(theCluster)
	  parallel::clusterEvalQ(theCluster, library('raster'))
      endCluster = TRUE
    }
  }
  
  if(verbose) {
	cat(date(), "\n")
    cat("running local-EM estimation for input bw\n")
  }
  
  # estimate risk (by partition, not continuous) for each bw combination
  forMoreArgs = list(
    x=cases[,countcol, drop=FALSE],
    lemObjects = lemObjects,
    tol = iterations$tol, 
    maxIter = iterations$maxIter,
    gpu = iterations$gpu,
    verbose=verbose)
  
  if(!is.null(theCluster)) {
    estList = parallel::clusterMap(
      theCluster,
      finalLemIter,
      bw = bwString,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE
    )
    
  } else {
    estList = mapply(
      finalLemIter,
      bw = bwString,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE
    )
  }

   if(any(unlist(lapply(estList, class)) == 'try-error') ) {
    warning("errors in local-em estimation")
  }
  
  estListRisk = lapply(estList, function(x) x$risk)
 
  riskDf = as.matrix(do.call(cbind, estListRisk))
  colnames(riskDf) = paste(
    rep(names(estList), unlist(lapply(estListRisk, function(xx) dim(xx)[2]))),
    unlist(lapply(
        estListRisk, colnames
      )), sep = '_')
  rownames(riskDf) = colnames(lemObjects$regionMat)
  
   newDf <- as.data.frame(riskDf[
    as.character(raster::levels(lemObjects$rasterFine)[[1]]$partition),])
	dimnames(newDf) = dimnames(riskDf)
  
  riskRaster = lemObjects$rasterFine
  levels(riskRaster)[[1]] =  as.data.frame(cbind(
      ID = raster::levels(lemObjects$rasterFine)[[1]]$ID,
      newDf))
  
  result = list(
    riskAll = riskRaster,
    smoothingMatrix = lemObjects
  )
  
	# final smoothing step
	result$riskEst = finalSmooth(
						x = result, 
						Slayers = dimnames(newDf)[[2]], 
						filename = filename, 
						ncores = theCluster)
	names(result$riskEst) = dimnames(newDf)[[2]]
   
	result$bw = bwString
 
   # done with the cluster
	if(endCluster) parallel::stopCluster(theCluster)

  if(verbose) {
	cat(date(), "\n")
    cat("done\n")
  }	
    
  return(result)

 
}
