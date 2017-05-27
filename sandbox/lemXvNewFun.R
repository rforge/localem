#' Additional options and packages
options(scipen = 22)

library('geostatsp')
library('Matrix')
library('rgeos')


#' Computes the likelihood cross-validation by removing spatial regions and estimating risk for those regions
#' This requires computations of new smoothing matrix for each training set
#' Sampling Ratio: 75/25% of data for training and testing sets
lemXvNew = function(
  polyCoarse, 
  polyFine, 
  cellsCoarse, 
  cellsFine, 
  bw, 
  Nxv = 1, 
  ncores = 1, 
  tol = 1e-6, 
  maxIter = 2000, 
  randomSeed = NULL, 
  verbose = FALSE
){
  
  #names of interest for regions in the coarse shapefile
  countcol = grep('^(count|cases)$', names(polyCoarse), value=TRUE, ignore.case=TRUE)
  if(length(countcol)){
    countcol = countcol[1]
  } else {
    countcol = grep(
      "^(id|name)", names(polyCoarse), 
      invert=TRUE, value=TRUE
    )[1]
  }
  
  idColX = grep("^id", names(polyCoarse), value=TRUE)
  if(length(idColX)) {
    idColX = idColX[1]
  } else {
    idColX = names(polyCoarse)[1]
  }

  #cross-validation scores
  if(!is.null(randomSeed)) {
    set.seed(randomSeed)
  }
  
  xvPct = 0.25
  
  xvCv = array(NA, 
               dim = c(length(bw) + 1, Nxv), 
               dimnames = list(c(paste("bw", bw, sep = ""),"bwInf"), 1:Nxv)
  )

  for(inXv in 1:Nxv) {
    
    #IDs for regions in the coarse shapefile
    idCoarse = polyCoarse[[idColX]]

        
    #test set
    ##random sample (25% used as testing data)
    xvId = sample(idCoarse, round(xvPct * length(idCoarse)), replace = FALSE)

    ##raster partition
    xvLemRaster = rasterPartition(polyCoarse = polyCoarse, 
                                  polyFine = polyFine, 
                                  cellsCoarse = cellsCoarse, 
                                  cellsFine = cellsFine, 
                                  bw = min(bw[1], 10000), 
                                  ncores = ncores, 
                                  idFile = 'idXv.grd', 
                                  offsetFile = 'offsetXv.grd', 
                                  verbose = FALSE)
 
    xvRasterFine = xvLemRaster$rasterFine
    
    Scells = ifelse(is.na(values(xvRasterFine[["cellCoarse"]])) | is.na(values(xvRasterFine[["idCoarse"]])), 
                    NA, 
                    paste("c", values(xvRasterFine[["cellCoarse"]]), 
                          "p", values(xvRasterFine[["idCoarse"]]), 
                          ".", values(xvRasterFine[["idFine"]]), 
                          sep = "")
    )

        
    #linkage for regions in the coarse and fine shapefiles (based on raster)
    polyFineCentre = SpatialPoints(polyFine, proj4string = polyFine@proj4string)
    polyFineCell = cellFromXY(xvRasterFine, polyFineCentre@coords)
    
    
    #training set
    ##offset and smoothing matrix
    trainPolyFine = polyFine
    
    trainRasterIdCoarse = which(values(xvRasterFine[['idCoarse']]) %in% match(xvId, idCoarse))
    trainPolyFine$expected[polyFineCell %in% trainRasterIdCoarse] = 0
    
    trainLemRaster = rasterPartition(polyCoarse = polyCoarse, 
                                     polyFine = trainPolyFine, 
                                     cellsCoarse = cellsCoarse, 
                                     cellsFine = cellsFine, 
                                     bw = bw, 
                                     ncores = ncores, 
                                     idFile = 'idTrain.grd', 
                                     offsetFile = 'offsetTrain.grd', 
                                     verbose = FALSE)
    trainLemRaster$rasterFine = xvRasterFine
    
    trainLemSmoothMat = smoothingMatrix(rasterObjects = trainLemRaster, 
                                        ncores = ncores, 
                                        verbose = FALSE)
  
    trainOffsetMat = trainLemSmoothMat$offsetMat
    
    trainSmoothingMat = trainLemSmoothMat$smoothingArray
    Spartitions = dimnames(trainSmoothingMat)[[1]]
    
    #test
    ##offsets
    xvRasterOffsets = xvLemRaster$offset[['offset']]
    values(xvRasterOffsets)[!(values(xvRasterFine[['idCoarse']]) %in% match(xvId, idCoarse))] = 0
    
    xvMeanOffsets = tapply(values(xvRasterOffsets), 
                           list(Scells), 
                           mean, na.rm=TRUE)
    
    xvMeanOffsets = xvMeanOffsets[match(Spartitions, names(xvMeanOffsets))]
    xvOffsetMat = Diagonal(length(Spartitions), xvMeanOffsets)
    
    
    #counts
    x = polyCoarse
    if(class(x) == 'SpatialPolygonsDataFrame') {
      x = data.frame(x)
    }
    
    ##fine raster did not include all regions in the coarse shapefile
    regionMat = trainLemSmoothMat$regionMat
        
    if(length(idCoarse) != dim(regionMat)[[2]]) {
      
      polyNeigh = spdep::poly2nb(trainLemSmoothMat$polyCoarse, row.names = idCoarse)
      
      idMatch = idCoarse[as.numeric(dimnames(regionMat)[[2]])]
      idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
      
      countCoarse = x[match(idMatch, x[[idColX]]),countcol]
      
      for(inD in idNotMatch) {
        
        polyNotMatch = trainLemSmoothMat$polyCoarse[idCoarse == inD,]
        idNeighNotMatch = idCoarse[values(intersect(trainLemSmoothMat$rasterFine[["idCoarse"]], polyNotMatch))]
        idNeighNotMatch = idNeighNotMatch[!is.na(idNeighNotMatch)]
        
        #if no match found in fine raster, use neighbouring coarse shapefile regions 
        if(length(idNeighNotMatch) == 0) {
          idNeighNotMatch = idCoarse[polyNeigh[[which(idCoarse == inD)]]]
          idNeighNotMatch = idMatch[idMatch %in% idNeighNotMatch]
        }
        
        #re-assign counts
        if(length(idNeighNotMatch) == 1) {
          
          countCoarse[idMatch == idNeighNotMatch] = 
            countCoarse[idMatch == idNeighNotMatch] + x[x[[idColX]] == inD,countcol]
          
        } else if(length(idNeighNotMatch) > 1) {
          
          #if conflict, assign counts to coarse shapefile region whose centroid is closest to the one of interest
          polyNeighNotMatch = trainLemSmoothMat$polyCoarse[idCoarse %in% idNeighNotMatch,]
          coordsNeighNotMatch = coordinates(rgeos::gCentroid(polyNeighNotMatch, byid = TRUE))
          
          coordsNotMatch = matrix(
            rep(coordinates(rgeos::gCentroid(polyNotMatch, byid = TRUE)), each = length(polyNeighNotMatch)), 
            nrow = length(polyNeighNotMatch), 
            ncol = 2, 
            dimnames = list(1:length(polyNeighNotMatch), c("x","y"))
          )
          
          distNeighNotMatch = apply((coordsNeighNotMatch - coordsNotMatch)^2, 1, sum)
          
          countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] = 
            countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] + x[x[[idColX]] == inD,countcol]
          
        }
      }
      
      idCoarse = idMatch
      
    } else {
      countCoarse = x[match(idCoarse, x[[idColX]]),countcol]
    }
  

    #risk estimation
    ##finite bandwidths
    if(verbose) {
      cat(date(), "\n")
      cat("obtaining CV for finite bandwidths\n")
    }
    
    ##training counts of interest
    trainCounts = as.matrix(countCoarse)
    trainCounts[idCoarse %in% xvId,] = 0
    
    ##starting values
    trainRegionOffset = crossprod(regionMat, trainOffsetMat)
    xvRegionOffset = crossprod(regionMat, xvOffsetMat) 
  
    startingValue = matrix(apply(trainRegionOffset, 2, sum), 
                           ncol(trainRegionOffset), 1)

    xvList = parallel::mclapply(dimnames(trainSmoothingMat)[[3]], 
                                xvLemEstOneBw, 
                                trainId = 1, 
                                trainCounts = trainCounts, 
                                regionOffset = trainRegionOffset, 
                                smoothingMat = trainSmoothingMat, 
                                regionXvOffset = xvRegionOffset, 
                                startingValue = startingValue, 
                                tol = tol, 
                                maxIter = maxIter, 
                                mc.cores = ncores
    )

		names(xvList) = dimnames(trainSmoothingMat)[[3]]
		xvMeans = simplify2array(xvList)[,1,]
		
		#cross-validation scores
		##finite bandwidths
		
		##test counts of interest
		xvCounts = matrix(countCoarse, length(countCoarse), dim(trainSmoothingMat)[3])
		xvCounts[!(idCoarse %in% xvId),] = 0
		
		xvRes = dpois(xvCounts, xvMeans, log = TRUE)
		dimnames(xvRes) = list(1:length(countCoarse), dimnames(trainSmoothingMat)[[3]])
		xvCv[dimnames(xvRes)[[2]],inXv] = -apply(xvRes, 2, sum)

		
		#risk estimation
    ##infinite bandwidth		
		if(verbose) {
      cat(date(), "\n")
      cat("obtaining CV for infinite bandwidth\n")
    }
  
    Scw2 = prod(res(xvRasterFine))
    SNcells = tapply(Scells, list(Scells), length)
    SNcells = SNcells[match(dimnames(trainSmoothingMat)[[1]], names(SNcells))]
    Sarea = SNcells * Scw2 
  
    Sexpected = trainOffsetMat@x * Sarea
  
    xvLambdaInf = as.matrix(Sarea * sum(trainCounts) / sum(Sexpected))
  
    #cross-validation scores
    ##infinite bandwidths
    xvMeansInf = xvRegionOffset %*% xvLambdaInf
    
    ##test counts of interest
    xvCountsInf = matrix(countCoarse, length(countCoarse), 1)
    xvCountsInf[!(idCoarse %in% xvId),] = 0
    
		xvResInf = dpois(as.vector(xvCountsInf), as.vector(xvMeansInf), log = TRUE)
    xvCv["bwInf",inXv] = -sum(xvResInf)
  }
  
  theXv = apply(xvCv, 1, sum)
  
  result = data.frame(
    bw = as.numeric(gsub("^bw", "", names(theXv))), 
    cv = theXv
  )
  
	if(verbose) { 
   cat(date(), "\n")
   cat("done\n")
  }
		
  return(result)
  
}




#' Additional functions from the 'localEM' package that are required but were not globally loaded
## Computes the risk estimation aggregated to the partitions for one iteration
oneLemIter = function(
  Lambda, smoothingMat, regionMat, offsetMat, counts,
  regionOffset = crossprod(regionMat, offsetMat)
){
  
  denom = regionOffset %*% Lambda
  
  em = crossprod(regionOffset, counts/denom) * Lambda
  em[as.vector(!is.finite(em))] = 0
  
  result = crossprod(smoothingMat, em)
  attributes(result)$em = em
  
  return(result)
}

## Computes the expected counts for the test set based on aggregated risk estimation and smoothing matrix of training set
xvLemEstOneBw = function(
  Dbw,
  trainId,
  trainCounts,
  regionOffset, # = crossprod(regionMat, offsetMat)
  smoothingMat,
  regionXvOffset, # = crossprod(regionMat, xvOffsetMat) 
  startingValue,
  tol,
  maxIter) {
  
  xvLambda = xvLemEst(
    trainId=trainId,
    trainCounts=trainCounts,
    regionOffset =regionOffset,
    smoothingMat = Matrix(drop(smoothingMat[,,Dbw])),
    startingValue = startingValue,
    tol=tol,
    maxIter = maxIter)
  
  as.matrix(regionXvOffset %*% xvLambda)
  
}

## Computes the aggregated risk estimation of the training set
xvLemEst = function(
  trainId,
  trainCounts,
  regionOffset,
  smoothingMat,
  offsetMat,
  startingValue = matrix(
    apply(regionOffset, 2, sum),
    ncol(regionOffset), length(trainId)
  ),
  tol,
  maxIter) {
  
  if(missing(trainId))
    trainId = 1:ncol(trainCounts)
  obsCounts = as.matrix(trainCounts[,trainId, drop=FALSE])
  
  oldRes = startingValue
  
  Diter = 1
  absdiff = 1
  
  while((absdiff > tol) && (Diter < maxIter)) {
    
    newRes = oneLemIter(
      Lambda = oldRes,
      smoothingMat = smoothingMat,
      regionOffset = regionOffset,
      counts = obsCounts
    )
    
    absdiff = mean(abs(oldRes - newRes))
    
    oldRes = newRes
    Diter = Diter + 1
  }
  
  result = as.matrix(newRes)
  
  return(result)
}
