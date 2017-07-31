#' @title Computes the relative risk estimation
#'
#' @description The \code{riskEst} function computes the estimations of the relative risk with fine raster resolution. 
#'
#' @param cases Spatial polygons, data frame or vector of case data
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param bw Vector of bandwidths specifying which smoothing matrix in \code{lemObjects} to use
#' @param iterations Convergence tolerance, number of iterations, and use of gpuR package for running local-EM recursions
#' @param verbose Verbose output
#' @param path Folder for storing rasters
#'
#' @details After using the \code{riskEst} function, the risk estimations are computed on a fine resolution based on the rasterization of the spatial polygons of population data.
#'
#' @return The \code{riskEst} function returns a raster of risk estimations for the input bandwidths.
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
#' lemRisk = riskEst(cases = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#'                    bw = c(10, 15) * 1000),  
#'                    ncores = 2, 
#'                    path = 'example', 
#'                    verbose = TRUE)
#'
#'}
#'
#' @export
riskEst = function(
  cases, 
  lemObjects, 
  bw, 
  ncores = 1, 
  iterations = list(tol = 1e-5, maxIter = 1000, gpu=FALSE), 
  verbose=FALSE, 
  path = getwd()
) {
  
  dir.create(path, showWarnings=FALSE, recursive=TRUE)

  # warning messages
  if(missing(lemObjects)) {
	stop("smoothing matrix not supplied")
  }
  
  if(missing(bw)) {
	cat("using bw from smoothing array\n")
	
	bwString = sort(unique(gsub("xv[[:digit:]]+$", "", lemObjects$bw)))
  } else { 
	  bwString = paste('bw', bw, sep='')
	  bwMatch = match(bwString, lemObjects$bw)
	  
	  if(any(is.na(bwMatch))) {
	  
		cat("cant match: ", bw[!is.na(bwMatch)], "\n")
		stop("cant match all bw to smoothing array")
	  }
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
  } else if(!is.null(ncores)) {
    if(ncores > 1) {
      if(verbose) cat("starting new cluster\n")
      theCluster = parallel::makeCluster(spec=ncores, type='PSOCK', methods=TRUE)
      parallel::setDefaultCluster(theCluster)
      endCluster = TRUE
    }
  }
  
  # checking structure of case data
   if(is.vector(cases)) {
    cases = data.frame(cases=cases)
  }
  #names of interest for regions in the coarse shapefile
  countcol = grep('^(count|cases)[[:digit:]]+?$', 
    names(cases), value=TRUE, ignore.case=TRUE)
  if(!length(countcol)){
    countcol = grep(
      "^(id|name)", names(cases), 
      invert=TRUE, value=TRUE
    )[1]
  }
  if(class(cases) == 'SpatialPolygonsDataFrame') {
    polyCoarse = cases
    cases = cases@data
  } else {
    polyCoarse = NULL
  }
  
  if(verbose) {
    cat("running local-EM estimation for input bandwidths\n")
  }
  
  # estimate risk (by partition, not continuous) for each bw/cv combinantion
  forMoreArgs = list(
    x=cases[,countcol, drop=FALSE],
    lemObjects = lemObjects,
    tol = iterations$tol, 
    maxIter = iterations$maxIter,
    gpu = iterations$gpu,
    verbose=verbose)
  
  if(!is.null(theCluster)) {
    
	parallel::clusterEvalQ(theCluster, library('raster'))

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
      )), sep='_')
  rownames(riskDf) = colnames(lemObjects$regionMat)
  
   newDf <- riskDf[
    as.character(raster::levels(lemObjects$rasterFine)[[1]]$partition),]
  
  riskRaster = lemObjects$rasterFine
  levels(riskRaster)[[1]] =  as.data.frame(cbind(
      ID = raster::levels(lemObjects$rasterFine)[[1]]$ID,
      newDf))
  
  result = list(
    riskAll = riskRaster,
    smoothingMatrix = lemObjects
  )
  
  # final smoothing step
   result$estimate = finalSmooth(
    x = result, 
	counts = rep(countcol, length(bwString)),
	bw = bwString, 
    filename = file.path(path, "risk.grd"),
    ncores = theCluster)
   names(result$estimate) = bwString
 
   # done with the cluster
	if(endCluster) parallel::stopCluster(theCluster)

  if(verbose) {
    cat("done\n")
  }	
    
  return(result)

 
}
