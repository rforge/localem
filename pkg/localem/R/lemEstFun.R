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
#  theFinalMat = smoothingFinalMat(
#    lemObjects=lemObjects,
#    bw=bw,
#    ncores=ncores
#  )

  if(verbose) {
    cat(date(), "\n")
    cat("obtaining risk estimation at final iteration\n")
  }

  #risk estimation
  theRisk = riskEst(
    x=x,
    lemObjects=lemObjects,
    bw=bwdimnames(theFinalMat$smoothingArray)[[3]],
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

# Computes the relative risk estimation on the raster of fine polygons
riskEst = function(x,
    lemObjects,
    bw,
		tol = 1e-6,
  	maxIter = 2000
) {
	
  regionMat = lemObjects$regionMat
  offsetMat = lemObjects$offsetMat
  smoothingMat = lemObjects$smoothingArray[,,
			paste('bw', bw, sep='')]
	
  idCoarse = lemObjects$polyCoarse$id
	
  #fine raster did not include all regions in the coarse shapefile
  if(length(idCoarse) != dim(regionMat)[2]) {
		
    polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)
		
    idMatch = idCoarse[as.numeric(dimnames(regionMat)[[2]])]
    idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
		
    obsCounts = as.matrix(x$count[match(idMatch, x[['id']])])
		
    for(inD in idNotMatch) {
			
      polyNotMatch = lemObjects$polyCoarse[idCoarse == inD,]
      idNeighNotMatch = idCoarse[values(intersect(lemObjects$rasterFine[["idCoarse"]], polyNotMatch))]
      idNeighNotMatch = idNeighNotMatch[!is.na(idNeighNotMatch)]
			
      #if no match found in fine raster, use neighbouring coarse shapefile regions
      if(length(idNeighNotMatch) == 0) {
        idNeighNotMatch = idCoarse[polyNeigh[[which(idCoarse == inD)]]]
        idNeighNotMatch = idMatch[idMatch %in% idNeighNotMatch]
      }
			
      #re-assign counts
      if(length(idNeighNotMatch) == 1) {
				
        obsCounts[idMatch == idNeighNotMatch,] =
          	obsCounts[idMatch == idNeighNotMatch,] + x$count[x$id == inD]
				
      } else if(length(idNeighNotMatch) > 1) {
				
        #if conflict, assign counts to coarse shapefile region whose centroid is closest to the one of interest
        polyNeighNotMatch = lemObjects$polyCoarse[idCoarse %in% idNeighNotMatch,]
        coordsNeighNotMatch = coordinates(rgeos::gCentroid(polyNeighNotMatch, byid = TRUE))
				
        coordsNotMatch = matrix(
          	rep(coordinates(rgeos::gCentroid(polyNotMatch, byid = TRUE)), each = length(polyNeighNotMatch)),
          	nrow = length(polyNeighNotMatch),
          	ncol = 2,
          	dimnames = list(1:length(polyNeighNotMatch), c("x","y"))
        )
				
        distNeighNotMatch = apply((coordsNeighNotMatch - coordsNotMatch)^2, 1, sum)
				
        obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] =
          	obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] + x$count[x$id == inD]
      }
    }
  } else {
		#fine raster does include all regions in the coarse shapefile
    obsCounts = as.matrix(x$count[match(idCoarse, x[['id']])])
  }
	
  #risk estimation for aggregated regions
  oldLambda = offsetMat %*%
    	matrix(1,
          nrow = dim(offsetMat)[1],
          ncol = 1,
          dimnames = list(dimnames(offsetMat)[[1]], bw)
    	)
	
  Diter = 1
  absdiff = 1
	
#	smoothingMat = smoothingMat / prod(res(lemObjects$offset))
	
  while((absdiff > tol) && (Diter < maxIter)) {
		
    Lambda = oneLemIter(
      	Lambda = oldLambda,
      	smoothingMat = smoothingMat,
      	regionMat = regionMat,
      	offsetMat = offsetMat,
      	counts = obsCounts)
		
    absdiff = mean(abs(oldLambda - Lambda))
		
    oldLambda = Lambda
    Diter = Diter + 1
  }
	
	
	lambdaMult = offsetMat %*% Lambda
	
	# convert lambdas to rasters
	
	resultRaster = lemObjects$rasterFine$idFine + 
			(maxValue(lemObjects$rasterFine$idFine)+10)*lemObjects$rasterFine$idCoarse
	
	resultRaster = resultRaster + (
				maxValue(resultRaster)+10) * lemObjects$rasterFine$cellCoarse
	
	resultRaster = ratify(resultRaster, count=TRUE)
	
	uniqueCells = levels(resultRaster)[[1]]$ID
	uniqueCells = match(uniqueCells, values(resultRaster))
	
	for(Dname in names(lemObjects$rasterFine))
		levels(resultRaster)[[1]][[Dname]] =
				values(lemObjects$rasterFine[[Dname]])[uniqueCells]
	
	levels(resultRaster)[[1]]$expected = 
			values(lemObjects$offset$offset)[uniqueCells] * 
			prod(res(lemObjects$offset))
	
	levels(resultRaster)[[1]]$partition = paste(
			'c', levels(resultRaster)[[1]]$cellCoarse,
			'p', levels(resultRaster)[[1]]$idCoarse,
			'.', levels(resultRaster)[[1]]$idFine,
			sep=''
	)
	
	levels(resultRaster)[[1]]$em = attributes(Lambda)$em[levels(resultRaster)[[1]]$partition,]		
	levels(resultRaster)[[1]]$bigLambda = lambdaMult[levels(resultRaster)[[1]]$partition,]		
	levels(resultRaster)[[1]]$lambda = levels(resultRaster)[[1]]$bigLambda / 
			(levels(resultRaster)[[1]]$COUNT * levels(resultRaster)[[1]]$expected)
	
	levels(resultRaster)[[1]]$emScale = levels(resultRaster)[[1]]$em / (
				levels(resultRaster)[[1]]$COUNT * prod(res(resultRaster)))
	
	emScale = deratify(resultRaster, 'emScale')
	
	emSmooth = focal(
			x=emScale, 
			w=lemObjects$focal$focal[[paste('bw',bw, sep='')]],
			na.rm=TRUE,pad=TRUE
	)
	
	result = emSmooth/lemObjects$offset[[paste('offset.bw',bw, sep='')]]
	
	
  return(result)
}
