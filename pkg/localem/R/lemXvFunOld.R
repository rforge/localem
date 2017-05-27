
lemXvOld = function(
  x, 
  lemObjects, 
  Nxv = 5, 
  ncores = 1, 
  tol = 1e-6, 
  maxIter = 2000, 
  randomSeed = NULL, 
  verbose = FALSE
){
  
  regionMat = lemObjects$regionMat
  
  trainOffsetMat = lemObjects$offsetMat * (Nxv - 1) / Nxv
  xvOffsetMat = lemObjects$offsetMat / Nxv
  
  trainSmoothingMat = lemObjects$smoothingArray / ((Nxv - 1) / Nxv)
  
  idCoarseCol = names(lemObjects$polyCoarse)[1]
  idCoarse = lemObjects$polyCoarse[[idCoarseCol]]
  
  countcol = grep('^(count|cases)$', names(x), value=TRUE, ignore.case=TRUE)
  if(length(countcol)){
    countcol = countcol[1]
  } else {
    countcol = grep(
      "^(id|name)", names(x), 
      invert=TRUE, value=TRUE
    )[1]
  }
  
  idColX = grep("^id", names(x), value=TRUE)
  if(length(idColX)) {
    idColX = idColX[1]
  } else {
    idColX = names(x)[1]
  }
  
  if(class(x) == 'SpatialPolygonsDataFrame') {
    x = data.frame(x)
  }
  
  #fine raster did not include all regions in the coarse shapefile
  if(length(idCoarse) != dim(regionMat)[[2]]) {
    
    polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)
    
    idMatch = idCoarse[as.numeric(dimnames(regionMat)[[2]])]
    idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
    
#   countCoarse = x$count[match(idMatch, x[['id']])]
    
    countCoarse = x[match(idMatch, x[[idColX]]),countcol]
    
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
        
#        countCoarse[idMatch == idNeighNotMatch] = 
#          countCoarse[idMatch == idNeighNotMatch] + x$count[x$id == inD]
        
        countCoarse[idMatch == idNeighNotMatch] = 
          countCoarse[idMatch == idNeighNotMatch] + x[x[[idColX]] == inD,countcol]
        
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
        
#        countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] = 
#          countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] + x$count[x$id == inD]
        
        countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] = 
          countCoarse[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)]] + x[x[[idColX]] == inD,countcol]
        
      }
    }
    
    idCoarse = idMatch
    
  } else {
    countCoarse = x[match(idCoarse, x[[idColX]]),countcol]
  }
  
  #partition area
  rasterFine = lemObjects$rasterFine
  
  # Scells = ifelse(is.na(values(rasterFine[["cellCoarse"]])) | is.na(values(rasterFine[["idCoarse"]])), 
  # NA, 
  # paste("c", values(rasterFine[["cellCoarse"]]), "p", values(rasterFine[["idCoarse"]]), sep = "")
  # )
  
  Scw2 = prod(res(rasterFine))
  Scells = ifelse(is.na(values(rasterFine[["cellCoarse"]])) | is.na(values(rasterFine[["idCoarse"]])), 
    NA, 
    paste("c", values(rasterFine[["cellCoarse"]]), 
      "p", values(rasterFine[["idCoarse"]]), 
      ".", values(rasterFine[["idFine"]]), 
      sep = "")
  )
  
  SNcells = tapply(Scells, list(Scells), length)
  SNcells = SNcells[match(dimnames(trainSmoothingMat)[[1]], names(SNcells))]
  Sarea = SNcells * Scw2 
  
  #test set
  if(!is.null(randomSeed)) {
    set.seed(randomSeed)
  }
  xvCounts = t(mapply(rmultinom, size = countCoarse, MoreArgs = list(n = 1, prob = rep(1/Nxv, Nxv))))
  dimnames(xvCounts) = list(idCoarse, 1:Nxv)
  
  #training set
  trainCounts = matrix(0,
    nrow = length(idCoarse), 
    ncol = Nxv,
    dimnames = list(idCoarse, 1:Nxv)
  )
  trainCounts[as.character(idCoarse),] = countCoarse
  trainCounts[dimnames(xvCounts)[[1]],] = trainCounts[dimnames(xvCounts)[[1]],] - xvCounts
  
  #risk estimation  
  xvLambda = array(NA, 
    dim = c(dim(trainSmoothingMat)[1], Nxv, dim(trainSmoothingMat)[3] + 1), 
    dimnames = list(dimnames(trainSmoothingMat)[[1]], dimnames(trainCounts)[[2]], c(dimnames(trainSmoothingMat)[[3]],"bwInf"))
  )
  
  xvRes = array(NA, 
    dim = c(dim(xvCounts)[1], Nxv, dim(xvLambda)[3]), 
    dimnames = list(dimnames(xvCounts)[[1]], dimnames(xvCounts)[[2]], dimnames(xvLambda)[[3]])
  )
  
#		stuff= xvLemEstOneBw(
#				Dbw=dimnames(trainSmoothingMat)[[3]][1],
#				trainId = 1:Nxv, 
#  		trainCounts = trainCounts, 
#  		regionMat = regionMat, 
#  		offsetMat = trainOffsetMat, 
#  		smoothingMat = trainSmoothingMat, 
#  		tol = tol, 
#  		maxIter = maxIter
#		)
  
  
  trainRegionOffset = crossprod(regionMat, trainOffsetMat)
  regionXvOffset = crossprod(regionMat, xvOffsetMat) 
  
  startingValue = matrix(
    apply(trainRegionOffset, 2, sum),
    ncol(trainRegionOffset), Nxv
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("obtaining CV for finite bandwidths\n")
  }
  
  xvList = parallel::mclapply(
    dimnames(trainSmoothingMat)[[3]],
    xvLemEstOneBw,
    trainId = 1:Nxv, 
    trainCounts = trainCounts, 
    regionOffset = trainRegionOffset, 
    smoothingMat = trainSmoothingMat, 
    regionXvOffset = regionXvOffset,
    startingValue=startingValue,
    tol = tol, 
    maxIter = maxIter,
    mc.cores=ncores
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("obtaining CV for infinite bandwidth\n")
  }
  
  names(xvList) = dimnames(trainSmoothingMat)[[3]]
  xvLambda = simplify2array(xvList)
  
  xvRes = xvLambda
  xvRes[,,] = dpois(as.vector(xvCounts), as.vector(xvRes), log = TRUE)
  
  #cross-validation scores
  xvCV = -apply(xvRes, 3, sum)
  
#  for(Dbw in dimnames(trainSmoothingMat)[[3]]) {
  
  
  #risk estimation of training set for finite bandwidths
#  xvLambda[,,Dbw] = simplify2array(
#    parallel::mclapply(1:Nxv, 
#      xvLemEst, 
#      trainCounts = trainCounts, 
#      regionMat = regionMat, 
#      offsetMat = trainOffsetMat, 
#      smoothingMat = trainSmoothingMat[,,Dbw], 
#      tol = tol, 
#      maxIter = maxIter, 
#      mc.cores=ncores
#    ))
#			xvLambda[,,Dbw] = xvLemEst(
#					trainId = 1:Nxv, 
#  			trainCounts = trainCounts, 
#  			regionMat = regionMat, 
#  			offsetMat = trainOffsetMat, 
#  			smoothingMat = trainSmoothingMat[,,Dbw], 
#  			tol = tol, 
#  			maxIter = maxIter)
  
  #likelihood cross-validation scores of test set for finite bandwidths
  #xvMeansOld = t(regionMat) %*% xvOffsetMat %*% xvLambda[,,Dbw]
#   xvMeans = crossprod(regionMat, xvOffsetMat) %*% xvLambda[,,Dbw]
#			
#			xvRes[,,Dbw] = dpois(as.vector(xvCounts), as.vector(xvMeans), log = TRUE)
#  }
  
  #risk estimation of training set for infinity bandwidth
  # Sexpected = diag(trainOffsetMat) * Sarea
  Sexpected = trainOffsetMat@x * Sarea
  
  xvLambdaInf = outer(Sarea, apply(trainCounts, 2, sum) / sum(Sexpected), "*")
  
  #likelihood cross-validation scores of test set for infinity bandwidth
  xvMeans = regionXvOffset %*% xvLambdaInf
  xvResInf = xvMeans
  xvResInf[,] = dpois(as.vector(xvCounts), as.vector(xvMeans), log = TRUE)
  
  
  #cross-validation scores
  xvCVInf = -sum(xvResInf)
  
  
  result = data.frame(
    bw = c(
      as.numeric(gsub("^bw", "", dimnames(xvLambda)[[3]])), Inf), 
    cv = c(xvCV, xvCVInf)
  )
  
  if(verbose) { 
    cat(date(), "\n")
    cat("done\n")
  }
  
  return(result)
  
}