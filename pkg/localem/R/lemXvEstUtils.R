# Generates the matrix of cross-validation sets
getXvMat = function(coarse, Nxv) {
  if(length(coarse)>1) {
    Ncoarse = length(coarse)
  } else{
    Ncoarse = coarse
    coarse = as.character(1:Ncoarse)
  }
  
  Matrix::sparseMatrix(
      i = 1:Ncoarse,
      j=sample(1:Nxv, Ncoarse, replace=TRUE),
      dimnames = list(coarse, as.character(1:Nxv))
  )
}


# Computes the risk estimation aggregated to the partitions for one iteration
oneLemIter = function(
  Lambda, smoothingMat, regionMat, offsetMat, counts,
		regionOffset = regionMat %*% offsetMat
){
# to do: write this all in C, with option to link to gpuR 
	
#  denomOld = t(regionMat) %*% offsetMat %*% Lambda

	 denom = regionOffset %*% Lambda
	
#  emOld = t(t(counts/denom) %*% t(regionMat) %*% offsetMat) * Lambda
 
  # the script M(Lambda) object
	em = crossprod(regionOffset, counts/denom) * Lambda
	
  em[as.vector(!is.finite(em))] = 0

#  resultOld = as.matrix(t(smoothingMat) %*% em)

result = crossprod(smoothingMat, em)

#  if(length(grep('vclMatrix',class(smoothingMat)))) {
#    em = gpuR::vclMatrix(as.matrix(em))
#    result = as.matrix(gpuR::crossprod(smoothingMat, em))
#  } else {
#    result = as.matrix(Matrix::crossprod(smoothingMat, em))    
#  }
	
#	attributes(result)$em = em
	
  return(result)
}


# Computes the expected counts for the test set based on aggregated risk estimation and smoothing matrix of training set
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
   smoothingMat = Matrix::Matrix(drop(smoothingMat[,,Dbw])),
			startingValue = startingValue,
   tol=tol,
   maxIter = maxIter)
	
	as.matrix(regionXvOffset %*% xvLambda)
	
}


# Computes the aggregated risk estimation of the training set
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


# Computes the risk estimation to the fine raster at the final iteration
oneLastStepSmooth = function(Dlayer, emScale, w, offsetSmooth) {
	result = focal(
			x=emScale[[Dlayer]],
			w=w, na.rm=TRUE, pad=TRUE
	)/offsetSmooth
	names(result) = as.character(Dlayer)
	result
}


# Computes the risk estimation to the partitions at the final iteration
finalLemIter = function(
  x, 
  lemObjects, 
  bw, 
  tol = 1e-5,
  maxIter = 1000,
  gpu=FALSE, 
  verbose=FALSE
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
    if(verbose) cat("loading smoothing matrix..")
    smoothingMat = matrix(raster::values(lemObjects$smoothingArray[[bw]]),
      nrow = nrow(lemObjects$smoothingArray),
      ncol = ncol(lemObjects$smoothingArray),
      byrow = FALSE,
      dimnames = list(
        lemObjects$partitions, lemObjects$partitions
      ))
    if(verbose) cat("done.\n")
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
    if(length(idCoarse) != dim(regionMat)[1]) {
      
      polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)
      
      idMatch = idCoarse[as.numeric(dimnames(regionMat)[[1]])]
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
  
  if(!all(rownames(smoothingMat) == colnames(regionMat))) {
    smoothingMat = smoothingMat[colnames(regionMat), colnames(regionMat)]
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
  
  regionOffset = regionMat %*% offsetMat
  
  crossprodFun = Matrix::crossprod
  

# a fix for zero counts zero offset give ratio of zero
  theZerosDenom = which(apply(regionOffset, 1, sum)<=0)
  zerosToAdd = matrix(0, nrow(obsCounts), ncol(obsCounts))		
  zerosToAdd[theZerosDenom,] = -1
  
  
  theZerosCount = which(as.matrix(obsCounts) <=0 )
  theZerosDenomX = sort(outer(
      theZerosDenom, 
      seq(from=0, by=nrow(obsCounts), len = ncol(obsCounts)), 
      '+'))
  
  if(!all(theZerosDenomX %in% theZerosCount)) {
    warning("zero denominator but non-zero count")
  }
  
  
  if(requireNamespace("gpuR", quietly=TRUE) & gpu) {
    if(verbose) cat("using GPU, ")
    smoothingMat = try(gpuR::vclMatrix(as.matrix(smoothingMat), type='double'), silent=TRUE)
    regionOffset = try(gpuR::vclMatrix(as.matrix(regionOffset), type='double'), silent=TRUE)
    oldLambda = try(gpuR::vclMatrix(as.matrix(oldLambda), type='double'), silent=TRUE)
    obsCounts = gpuR::vclMatrix(as.matrix(obsCounts), type='double')
    zerosToAdd = gpuR::vclMatrix(zerosToAdd, type='double')
    crossprodFun = gpuR::crossprod
    
    if(class(smoothingMat) == 'try-error'){
      warning('putting smoothing matrix in gpu failed, probably out of memory\n')
    }
  } else {
    smoothingMat = Matrix::Matrix(smoothingMat)
  }

  Diter = 1
  absdiff = Inf
  
  if(verbose) cat("starting lem, bandwidth", bwString, '\n')
  

  while((absdiff > tol) && (Diter < maxIter) ) {
    # to do: use gpuR's cpp_gpuMatrix_custom_igemm to reuse memory
    
    denom = regionOffset %*% oldLambda
    denom = denom + zerosToAdd
    
    countsDivDenom = obsCounts/denom
    
    roCprod = crossprodFun(regionOffset, countsDivDenom)
    
    em = roCprod * oldLambda
    
    Lambda = crossprodFun(smoothingMat, em)
    
    lDiff = oldLambda - Lambda    
    
    absdiff = max(abs(c(max(lDiff), min(lDiff))))
    if(all(is.na(absdiff))) {
      warning("missing values in Lambda")
      absdiff = -Inf
    }
    
    lDiff = oldLambda
    oldLambda = Lambda
    Lambda = lDiff
    Diter = Diter + 1
    if(verbose) cat(".")
  }
  if(verbose) cat("\ndone lem,", Diter, 'iterations\n')
  colnames(Lambda) = colnames(obsCounts)

  littleLambda = solve(partitionAreasMat) %*% as.matrix(Lambda)

  # expected count using full offsets, not xv offests
  expectedCoarseRegions = 
    (regionMat %*% lemObjects$offsetMat[['offset']]) %*% 
    littleLambda

  return(list(
      expected = as.matrix(expectedCoarseRegions),
      risk = as.matrix(littleLambda)))
}


# Computes the risk estimation to the fine raster by smoothing the risk of the partitions at the final iteration
finalSmooth = function(
  x, 
  Slayers, 
  ncores = 1, 
  filename) {
  
  toSmooth = x$riskAll
  levels(toSmooth)[[1]] = levels(toSmooth)[[1]][, c("ID", Slayers)]
  
  deratifyFile = tempfile("deratify", tempdir(), ".grd")
  
  toSmooth = deratify(toSmooth, filename = deratifyFile)
  
  endCluster = FALSE
  theCluster = NULL
  if(length(grep("cluster", class(ncores))) ) {
    theCluster = ncores
  } else if(!is.null(ncores)) {
    if(ncores > 1) {
      theCluster = parallel::makeCluster(spec=ncores, type='PSOCK', methods=TRUE)
      parallel::setDefaultCluster(theCluster)
      endCluster = TRUE
    }
  }
  
  xOrig = x
  result = focalMult(
    x=toSmooth, 
    w=xOrig$smoothingMatrix$focal$focal, 
    edgeCorrect = TRUE,
    filename = filename,
    cl = theCluster
  )
  
  if(endCluster)  
    parallel::stopCluster(theCluster)
    
  return(result)
  
}
