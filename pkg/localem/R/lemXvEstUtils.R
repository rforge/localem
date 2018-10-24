###
# SINGLE MAP
###

# Generates the matrix of cross-validation sets
getXvMat = function(coarse, coarseRaster, offsetRaster, Nxv) {

  if(length(coarse)>1) {
    Ncoarse = length(coarse)
  } else {
    Ncoarse = coarse
    coarse = as.character(1:Ncoarse)
  }

  expectedCoarse = aggregate(values(offsetRaster[['offset']]) * prod(res(offsetRaster)),
                             by = list(values(coarseRaster[['idCoarse']])),
                             FUN = 'sum')
  colnames(expectedCoarse) = c('idCoarse','expected')
  expectedCoarse = merge(data.frame(idCoarse = 1:Ncoarse), expectedCoarse,
                         by = 'idCoarse',
                         all = TRUE)
  expectedCoarse$expected[is.na(expectedCoarse$expected)] = 0

  xvExpected = expectedCoarse[order(expectedCoarse$expected),]
  xvExpected$idXv = (1:Ncoarse %% Nxv) + 1

  Matrix::sparseMatrix(i = xvExpected$idCoarse,
                       j = xvExpected$idXv,
                       dimnames = list(coarse, as.character(1:Nxv))
  )

}

# getXvMat = function(coarse, Nxv) {
#   if(length(coarse)>1) {
#     Ncoarse = length(coarse)
#   } else{
#     Ncoarse = coarse
#     coarse = as.character(1:Ncoarse)
#   }
#
#   Matrix::sparseMatrix(
#       i = 1:Ncoarse,
#       j=sample(1:Nxv, Ncoarse, replace=TRUE),
#       dimnames = list(coarse, as.character(1:Nxv))
#   )
# }

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


###
# MULTIPLE MAPS
###

# Generates the matrix of cross-validation sets (based on regions with similar expected counts) for each case map
getXvMatOneMap = function(coarse, coarseRaster, offsetRaster, Nxv) {
  
  if(length(coarse) > 1) {
    Ncoarse = length(coarse)
  } else {
    Ncoarse = coarse
    coarse = as.character(1:Ncoarse)
  }
  
  expectedCoarse = aggregate(values(offsetRaster[['offset']]) * prod(res(offsetRaster)),
                             by = list(values(coarseRaster[['idCoarse']])),
                             FUN = 'sum')
  colnames(expectedCoarse) = c('idCoarse','expected')
  expectedCoarse = merge(data.frame(idCoarse = 1:Ncoarse), expectedCoarse,
                         by = 'idCoarse',
                         all = TRUE)
  expectedCoarse$expected[is.na(expectedCoarse$expected)] = 0
  
  xvExpected = expectedCoarse[order(expectedCoarse$expected),]
#  xvExpected$idXv = (1:Ncoarse %/% floor(Ncoarse / (Nxv + 1)) + 1
#  xvExpected$idXv = (1:Ncoarse %% Nxv) + 1
  xvExpected$idXv = (1:Ncoarse %/% (Ncoarse / Nxv)) + 1
  xvExpected$idXv[xvExpected$idXv > Nxv] = Nxv
  
  Matrix::sparseMatrix(i = xvExpected$idCoarse,
                       j = xvExpected$idXv,
                       dimnames = list(coarse, as.character(1:Nxv))
  )
}

# Re-calibrate the matrix of cross-validation sets for all maps
## For a cross-validation set for map of interest, exclude overlaying case regions of other maps for risk estimation
getXvMatUpdate = function(polyCoarse, xvMat, Nxv) {
  
  resList = xvMat
  for(inM in 1:length(polyCoarse)) {
    
    polyCoarseInt = polyCoarse[[inM]]
	
    xvMatInt = xvMat[[inM]]
    
    idMatInt = as.numeric(dimnames(xvMatInt)[[2]]) - 1
    # idMatInt = (idMatInt %/% Nxv) + 1
    idMatInt = (idMatInt %% length(polyCoarse)) + 1
    idMatInt = which(idMatInt == inM)
    
    for(inN in 1:length(polyCoarse)) {
      
      if(inN != inM) {
        
        polyCoarseOth = polyCoarse[[inN]]
		polyCoarseOth$idCoarse = 1:length(polyCoarseOth)
		
        xvMatOth = xvMat[[inN]]
    
        for(inX in idMatInt) {
        
          # idMatOth = sp::over(
          #   sp::spsample(polyCoarseInt[xvMatInt[,inX],], n = 250, type = 'random'), 
          #   polyCoarseOth)$id
		  
		  idMatNpts = 1000
		  idMatCutoff = 0.25 * idMatNpts
		  spInt = sapply(which(xvMatInt[,inX]), 
					function(x) sp::spsample(polyCoarseInt[x,], n = idMatNpts, type = 'regular'))
		  idMatOth = unique(unlist(
						lapply(spInt, 
							function(x, cutoff) {
								result = table(sp::over(x, polyCoarseOth)$idCoarse)
								result = result[result > cutoff]
								return(as.numeric(names(result)))
							}, 
						cutoff = idMatCutoff)))
	          
          xvMatOth[,inX] = FALSE
          xvMatOth[idMatOth,inX] = TRUE
          
          resList[[inN]][,inX] = xvMatOth[,inX]
        }
      }
    }
  }

  return(resList)
}

# Generate the observed counts for regions of interest of the case map
getObsCounts = function(
  x, 
  polyCoarse, 
  regionMat
) {
  
  idCoarse = 1:nrow(polyCoarse)
  idMatch = as.numeric(which(apply(regionMat, 1, sum) > 0))
  
  #fine raster did not include all regions in the coarse shapefile
  result = x
  if(length(idCoarse) != length(idMatch)) {
    
    polyNeigh = spdep::poly2nb(polyCoarse, row.names = idCoarse)
    
    idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
    result[idNotMatch,] = 0
    
    for(inM in idNotMatch) {
      
      idNeighNotMatch = polyNeigh[[inM]]
      idNeighNotMatch = idNeighNotMatch[!(idNeighNotMatch %in% idNotMatch)]
      
      nNeighNotMatch = length(idNeighNotMatch)
      
      #re-assign counts
      if(nNeighNotMatch > 0) {
        
        for(inN in 1:ncol(result)) {
          
          obsNeighNotMatch = stats::rmultinom(n = 1,
                                              size = x[inM,inN],
                                              prob = rep(1/nNeighNotMatch, nNeighNotMatch))
          result[idMatch %in% idNeighNotMatch,inN] =
            result[idMatch %in% idNeighNotMatch,inN] + obsNeighNotMatch
        }
      }
    }
  }
  
  return(result)
}
