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

  #single data
  if(is.vector(x)) {
    x = matrix(x, ncol = 1)
    colnames(x) = 'count'
  }
  # multiple data
  if(is.data.frame(x)) {
    x = as.matrix(x)
  }

  idCoarse = 1:nrow(x)
  idMatch = as.numeric(which(apply(regionMat, 1, sum) > 0))

  obsCounts = x[idCoarse,,drop=FALSE]

  #fine raster did not include all regions in the coarse shapefile
  if(length(idCoarse) != length(idMatch)) {

    polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)

    idNotMatch = idCoarse[!(idCoarse %in% idMatch)]
    obsCounts[idNotMatch,] = 0

    for(inM in idNotMatch) {

      idNeighNotMatch = polyNeigh[[inM]]
      idNeighNotMatch = idNeighNotMatch[!(idNeighNotMatch %in% idNotMatch)]

      nNeighNotMatch = length(idNeighNotMatch)

      #re-assign counts
      if(nNeighNotMatch > 0) {

        for(inN in 1:ncol(obsCounts)) {

          obsNeighNotMatch = stats::rmultinom(n = 1,
                                              size = x[inM,inN],
                                              prob = rep(1/nNeighNotMatch, nNeighNotMatch))
          obsCounts[idMatch %in% idNeighNotMatch,inN] =
            obsCounts[idMatch %in% idNeighNotMatch,inN] + obsNeighNotMatch
        }
      }
    }
  }

# # multiple observations
  # if(is.data.frame(x)) x = as.matrix(x)
  # if(is.matrix(x)) {

    # idCoarse = 1:nrow(x)

    # if(length(idCoarse) != dim(regionMat)[1]) {

      # #fine raster did not include all regions in the coarse shapefile
      # idMatch = idCoarse[as.numeric(dimnames(regionMat)[[1]])]

      # obsCounts = x[idMatch,, drop=FALSE]

    # } else {
      # obsCounts = x[idCoarse,,drop=FALSE]
    # }

  # } else if(class(x) == 'SpatialPolygonsDataFrame') {
    # polyCoarse = x
    # idCoarseCol = names(lemObjects$polyCoarse)[1]
    # idCoarse = lemObjects$polyCoarse@data[[idCoarseCol]]

    # x = data.frame(x)

    # countcol = grep('^(count|cases)$', names(x), value=TRUE, ignore.case=TRUE)
    # if(length(countcol)){
      # countcol = countcol[1]
    # } else {
      # countcol = grep("^(id|name)", names(x), invert=TRUE, value=TRUE)[1]
    # }

    # idColX = grep("^id", names(x), value=TRUE)
    # if(length(idColX)) {
      # idColX = idColX[1]
    # } else {
      # idColX = names(x)[1]
    # }

    # #fine raster did not include all regions in the coarse shapefile
    # if(length(idCoarse) != dim(regionMat)[1]) {

      # polyNeigh = spdep::poly2nb(lemObjects$polyCoarse, row.names = idCoarse)

      # idMatch = idCoarse[as.numeric(dimnames(regionMat)[[1]])]
      # idNotMatch = idCoarse[!(idCoarse %in% idMatch)]

      # obsCounts = as.matrix(x[match(idMatch, x[[idColX]]),countcol])
      # dimnames(obsCounts) = list(idMatch, bw)

      # for(inD in idNotMatch) {

        # polyNotMatch = lemObjects$polyCoarse[idCoarse == inD,]
        # idNeighNotMatch = idCoarse[values(intersect(lemObjects$rasterFine[["idCoarse"]], polyNotMatch))]
        # idNeighNotMatch = idNeighNotMatch[!is.na(idNeighNotMatch)]

        # #if no match found in fine raster, use neighbouring coarse shapefile regions
        # if(length(idNeighNotMatch) == 0) {
          # idNeighNotMatch = idCoarse[polyNeigh[[which(idCoarse == inD)]]]
          # idNeighNotMatch = idMatch[idMatch %in% idNeighNotMatch]
        # }

        # #re-assign counts
        # if(length(idNeighNotMatch) == 1) {

          # obsCounts[idMatch == idNeighNotMatch,] =
              # obsCounts[idMatch == idNeighNotMatch,] + x[x[[idColX]] == inD,countcol]

        # } else if(length(idNeighNotMatch) > 1) {

          # #if conflict, assign counts to coarse shapefile region whose centroid is closest to the one of interest
          # polyNeighNotMatch = lemObjects$polyCoarse[idCoarse %in% idNeighNotMatch,]
          # coordsNeighNotMatch = sp::coordinates(rgeos::gCentroid(polyNeighNotMatch, byid = TRUE))

          # coordsNotMatch = matrix(
              # rep(coordinates(rgeos::gCentroid(polyNotMatch, byid = TRUE)), each = length(polyNeighNotMatch)),
              # nrow = length(polyNeighNotMatch),
              # ncol = 2,
              # dimnames = list(1:length(polyNeighNotMatch), c("x","y"))
          # )

          # distNeighNotMatch = apply((coordsNeighNotMatch - coordsNotMatch)^2, 1, sum)

          # obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] =
              # obsCounts[idMatch == idNeighNotMatch[which.min(distNeighNotMatch)],] + x[x[[idColX]] == inD,countcol]
        # }
      # }
    # } else {

      # obsCounts = as.matrix(x[match(idCoarse, x[[idColX]]),countcol])
      # dimnames(obsCounts) = list(idCoarse, countcol)
    # }
  # } else { # x is a vector
    # idCoarse = 1:length(x)
    # obsCounts = as.matrix(x, ncol=1)
    # colnames(obsCounts) = 'y'
  # }

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

    oldLambda = Lambda
    Diter = Diter + 1

#    if(verbose) cat(".")
  }
  if(verbose) cat("\ndone lem,", Diter, 'iterations\n')
  colnames(Lambda) = colnames(obsCounts)

  littleLambda = solve(partitionAreasMat) %*% as.matrix(Lambda)

  # expected count using full offsets, not xv offests
  offsetMatFull = lemObjects$offsetMat[['offset']]
  offsetMatFull = offsetMatFull[colnames(regionMat), colnames(regionMat)]
  expectedCoarseRegions = (regionMat %*% offsetMatFull) %*% as.matrix(Lambda)

  return(list(
          expected = as.matrix(expectedCoarseRegions),
          risk = as.matrix(littleLambda)))
}


