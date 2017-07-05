#' @title Computes the relative risk estimation
#'
#' @description The \code{riskEst} function computes the estimations of the relative risk with fine raster resolution. 
#'
#' @param x Spatial polygons of case data
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param bw Bandwidth specifying which smoothing matrix in \code{lemObjects} to use
#' @param ncores Number of cores/threads for parallel processing
#' @param tol Tolerance for convergence
#' @param maxIter Maximum number of iterations for convergence
#'
#' @details After using the \code{riskEst} function, the risk estimations are computed on a fine resolution based on the rasterization of the spatial polygons of population data.

#' @return The \code{riskEst} function returns a raster of risk estimations for the input bandwidth.
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
#'                              path = tempdir(), 
#'                              verbose = TRUE)
#'
#'
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'                                 ncores = 2, 
#'                                 verbose = TRUE)
#'
#' lemRisk = riskEst(x = kentuckyCounty, 
#'                    lemObjects = lemSmoothMat, 
#'                    bw = 15 * 1000) 
#'
#'}
#'
#' @export
riskEst = function(
    x, 
    lemObjects, 
    bw, 
    tol = 1e-6, 
    maxIter = 2000
) {
  
  if(length(bw) > 1) {
    stop("Bandwidth must be length 1")
  }
  
  if(is.character(bw)) {
    bwString = bw
    bw = match(bw, names(lemObjects$smoothingArray))
    if(is.na(bw))
      bw = match(bw, lemObjects$bw)
    if(is.na(bw))
      warning("cant match bw to smoothing array")
  } else {
    bwOrig = bw
    bw = paste('bw', bw, sep='')
    bw = match(bw, names(lemObjects$smoothingArray))
    bwString = names(lemObjects$smoothingArray)[bw]
    if(is.na(bw)) {
      bw = match(bw, lemObjects$bw)
      bwString = lemObjects$bw[bw]
    }
    if(is.na(bw)) {
      bwString = lemObjects$bw[bwOrig] 
      bw = bwOrig
    }
  }
  
  regionMat = lemObjects$regionMat
  smoothingMat = lemObjects$smoothingArray[[bw]]

  if(class(smoothingMat) == 'RasterLayer') {
  
  smoothingMat = Matrix::Matrix(raster::values(lemObjects$smoothingArray[[bw]]),
      nrow = nrow(lemObjects$smoothingArray),
      ncol = ncol(lemObjects$smoothingArray),
      byrow = FALSE,
      dimnames = list(
          lemObjects$partitions, lemObjects$partitions
      ))
  }
  
  if(length(grep("xv[[:digit:]]+$", bwString))) {
    offsetMat = lemObjects$offsetMat[[
        paste('xvOffset', 
            gsub('^bw[[:digit:]]+xv', '', bwString), sep='')]]
  } else {
    offsetMat = lemObjects$offsetMat[['offset']]
  }
  
  
# multiple observations
  if(is.data.frame(x)) x = as.matrix(x)
  if(is.matrix(x)) {
    
    idCoarse = 1:nrow(x)
    
    if(length(idCoarse) != dim(regionMat)[1]) {		
      
      #fine raster did not include all regions in the coarse shapefile
      idMatch = idCoarse[as.numeric(dimnames(regionMat)[[1]])]
      
      obsCounts = x[idMatch,, drop=FALSE]
      
    } else {
      obsCounts = x[idCoarse,,drop=FALSE]
    }
    
  } else if(class(x) == 'SpatialPolygonsDataFrame') {
    polyCoarse = x
    idCoarseCol = names(lemObjects$polyCoarse)[1]
    idCoarse = lemObjects$polyCoarse@data[[idCoarseCol]]
    
    x = data.frame(x)
    
    countcol = grep('^(count|cases)$', names(x), value=TRUE, ignore.case=TRUE)
    if(length(countcol)){
      countcol = countcol[1]
    } else {
      countcol = grep("^(id|name)", names(x), invert=TRUE, value=TRUE)[1]
    }
    
    idColX = grep("^id", names(x), value=TRUE)
    if(length(idColX)) {
      idColX = idColX[1]
    } else {
      idColX = names(x)[1]
    }
    
    #fine raster did not include all regions in the coarse shapefile
    if(length(idCoarse) != dim(regionMat)[2]) {
      
      polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)
      
      idMatch = idCoarse[as.numeric(dimnames(regionMat)[[2]])]
      idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
      
      obsCounts = as.matrix(x[match(idMatch, x[[idColX]]),countcol])
      dimnames(obsCounts) = list(idMatch, bw)
      
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
              obsCounts[idMatch == idNeighNotMatch,] + x[x[[idColX]] == inD,countcol]
          
        } else if(length(idNeighNotMatch) > 1) {
          
          #if conflict, assign counts to coarse shapefile region whose centroid is closest to the one of interest
          polyNeighNotMatch = lemObjects$polyCoarse[idCoarse %in% idNeighNotMatch,]
          coordsNeighNotMatch = sp::coordinates(rgeos::gCentroid(polyNeighNotMatch, byid = TRUE))
          
          coordsNotMatch = matrix(
              rep(coordinates(rgeos::gCentroid(polyNotMatch, byid = TRUE)), each = length(polyNeighNotMatch)),
              nrow = length(polyNeighNotMatch),
              ncol = 2,
              dimnames = list(1:length(polyNeighNotMatch), c("x","y"))
          )
          
          distNeighNotMatch = apply((coordsNeighNotMatch - coordsNotMatch)^2, 1, sum)
          
          obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] =
              obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] + x[x[[idColX]] == inD,countcol]
        }
      }
    } else {
      
      obsCounts = as.matrix(x[match(idCoarse, x[[idColX]]),countcol])
      dimnames(obsCounts) = list(idCoarse, countcol)
    }
  } else { # x is a vector
    idCoarse = 1:length(x)
    obsCounts = as.matrix(x, ncol=1)
    colnames(obsCounts) = 'y'
  }
  xvSet = gsub("^bw[[:digit:]]+(xv)?", "", bwString)
  if(nchar(xvSet)) {
    obsCounts = obsCounts * (!lemObjects$xv[,xvSet])
  }
  
  # Lambda = \int_{Q_\ell} O(u) du, = O_{\ell} (Lambda in Biometrics paper)
  
  partitionAreasMat = Matrix::Diagonal(
      length(lemObjects$partitionAreas),
      lemObjects$partitionAreas)
  dimnames(partitionAreasMat) = list(
      names(lemObjects$partitionAreas), 
      names(lemObjects$partitionAreas)
  )
  
  # make sure everything's ordered correctly
  
  offsetMat = offsetMat[colnames(regionMat), colnames(regionMat)]
  partitionAreasMat = partitionAreasMat[colnames(regionMat), colnames(regionMat)]
 
  smoothingMat = smoothingMat[colnames(regionMat), colnames(regionMat)]
  if(requireNamespace("gpuR", quietly=TRUE)) {
#    smoothingMat = gpuR::gpuMatrix(as.matrix(smoothingMat))
    smoothingMat = gpuR::vclMatrix(as.matrix(smoothingMat))
  }

  
  
  #risk estimation for aggregated regions
  oldLambda = partitionAreasMat %*%
      matrix(1,
          nrow = dim(offsetMat)[1],
          ncol = ncol(obsCounts),
          dimnames = list(
              dimnames(offsetMat)[[1]], 
              colnames(obsCounts))
      )
  
  Diter = 1
  absdiff = Inf
  
  
  regionOffset = regionMat %*% offsetMat
  
  
  while((absdiff > tol) && (Diter < maxIter)) {
    
    Lambda = oneLemIter(
        Lambda = oldLambda,
        smoothingMat = smoothingMat,
        regionOffset = regionOffset,
        counts = obsCounts)
    
    absdiff = mean(abs(oldLambda - Lambda))
    if(all(is.na(absdiff))) {
      warning("missing values in Lambda")
      absdiff = -Inf
    }
    
    oldLambda = Lambda
    Diter = Diter + 1
  }
  colnames(Lambda) = colnames(obsCounts)
  
  littleLambda = solve(partitionAreasMat) %*% Lambda
  
    # expected count using full offsets, not xv offests
    expectedCoarseRegions = 
        (regionMat %*% lemObjects$offsetMat[['offset']]) %*% 
        littleLambda
    
    return(list(
            expected = as.matrix(expectedCoarseRegions),
            risk = as.matrix(littleLambda)))
}
