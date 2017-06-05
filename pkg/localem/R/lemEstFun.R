#' @title Computes the relative risk estimation on the raster of fine polygons
#'
#' @description The \code{riskEst} function computes the estimations of the relative risk on the cells of the fine raster. 
#'
#' @param x Spatial polygons of case data
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition and smoothed offsets
#' @param bw Bandwidth specifying which smoothing matrix in \code{lemObjects} to use
#' @param ncores Number of cores/threads for parallel processing
#' @param tol Tolerance for convergence
#' @param maxIter Maximum number of iterations for convergence
#' @param filename Passed to writeRaster
#'
#' @details After using the \code{riskEst} function, the risk estimations are computed on the raster cells of the fine polygons.

#' @return The \code{riskEst} function returns a raster of risk estimations for the input bandwidth.
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
#' rCol = mapmisc::colourScale(lemRisk, 
#'                            breaks = 10, style = 'equal', dec = 1)
#' plot(lemRisk, 
#'     col = rCol$col, breaks = rCol$breaks, 
#'     legend = TRUE)
#'}
#'
#' @export
riskEst = function(
    x, 
    lemObjects, 
    bw, 
    tol = 1e-6, 
    maxIter = 2000,
    ncores = 1,
    type = c('final','unsmoothed','expected'), 
    filename = ''
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
  smoothingMat = t(as.matrix(lemObjects$smoothingArray[[bw]]))
  dimnames(smoothingMat) = list(
      lemObjects$partitions, lemObjects$partitions
  )
  
  if(length(grep("xv[[:digit:]]+$", bwString))) {
    offsetMat = lemObjects$offsetMat[[
        paste('xvOffset', 
            gsub('^bw[[:digit:]]+xv', '', bwString), sep='')]]
  } else {
    offsetMat = lemObjects$offsetMat[['offset']]
  }
  
  
  # simulated bootstrap data 
  if(is.matrix(x)) {
    
    idCoarse = 1:nrow(x)
    
    if(length(idCoarse) != dim(regionMat)[2]) {		
      
      #fine raster did not include all regions in the coarse shapefile
      idMatch = idCoarse[as.numeric(dimnames(regionMat)[[2]])]
      
      obsCounts = x[idMatch,, drop=FALSE]
      
    } else {
      
      obsCounts = x[idCoarse,, drop=FALSE]
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
          coordsNeighNotMatch = coordinates(rgeos::gCentroid(polyNeighNotMatch, byid = TRUE))
          
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
      dimnames(obsCounts) = list(idCoarse, bw)
    }
  } else { # x is a vector
    idCoarse = 1:length(x)
    obsCounts = as.matrix(x, ncol=1)  
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
  offsetMat = offsetMat[colnames(regionMat), colnames(regionMat)]
  smoothingMat = smoothingMat[colnames(regionMat), colnames(regionMat)]
  partitionAreasMat = partitionAreasMat[rownames(regionMat), rownames(regionMat)]
  
  
  regionOffset = Matrix::crossprod(regionMat, offsetMat) # 
  
  # make sure everything's ordered correctly
  
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
  absdiff = 1
  
#	smoothingMat = smoothingMat / prod(res(lemObjects$offset))
  
  
  while((absdiff > tol) && (Diter < maxIter)) {
    
    Lambda = oneLemIter(
        Lambda = oldLambda,
        smoothingMat = smoothingMat,
        regionOffset = regionOffset,
        counts = obsCounts)
    
    absdiff = mean(abs(oldLambda - Lambda))
    
    oldLambda = Lambda
    Diter = Diter + 1
  }
  
  littleLambda = solve(partitionAreasMat) %*% Lambda
  
  if(type[1] == 'expected') {
    # expected count using full offsets, not xv offests
    expectedCoarseRegions = 
        Matrix::crossprod(regionMat, lemObjects$offsetMat[['offset']]) %*% 
        littleLambda
    return(list(
            expected = expectedCoarseRegions,
            risk = littleLambda))
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
  
  levelsEm = as.matrix(
      attributes(Lambda)$em
  )[levels(resultRaster)[[1]]$partition,,drop=FALSE]
  
  emScale = levelsEm /(
        levels(resultRaster)[[1]]$COUNT * prod(res(resultRaster))
        )
  
  colnames(levelsEm) = 
      paste('em.', colnames(levelsEm), sep='')
  colnames(emScale) =
      paste('emScale.', colnames(emScale), sep='')
  
  
  levels(resultRaster)[[1]] = cbind(
      levels(resultRaster)[[1]],
      emScale#, levelsEm, bigLambda, littleLambda
  )
  
  
#	bigLambda = as.matrix(lambdaMult)[levels(resultRaster)[[1]]$partition,]
  
#	littleLambda = bigLambda / (levels(resultRaster)[[1]]$COUNT * levels(resultRaster)[[1]]$expected)
  
#	colnames(bigLambda) = paste('bigLambda.', colnames(bigLambda), sep='')
#	colnames(littleLambda) = paste('lambda.', colnames(littleLambda), sep='')
  
  emScale = deratify(resultRaster, 
      grep('^emScale', 
          colnames(levels(resultRaster)[[1]]), 
          value=TRUE)
  )
  
  
  wMat=lemObjects$focal$focal[[paste('bw',bw, sep='')]]
  offsetSmooth = lemObjects$offset[[paste('offset.bw',bw, sep='')]]
  
#	stuff = oneLastStepSmooth(
#			Dlayer = names(emScale)[100],
#			emScale=emScale,
#			w=wMat,
#			offsetSmooth=offsetSmooth			
#			)
  
  emSmooth = parallel::mcmapply(
      oneLastStepSmooth,
      Dlayer = names(emScale),
      MoreArgs = list(
          emScale=emScale,
          w=wMat,
          offsetSmooth=offsetSmooth
      ),
      SIMPLIFY=FALSE,
#			USE.NAMES=FALSE,
      # for some reason ncores=4 fails
#			mc.cores=pmin(2,ncores)
      mc.cores=ncores
  )
  theNames = names(emSmooth)
  names(emSmooth) = NULL
  result = do.call(brick, 
      c(emSmooth, list(filename=filename)))		
  names(result) = gsub("^emScale", "risk", theNames)
  
  return(result)
}
