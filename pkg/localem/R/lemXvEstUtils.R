# Computes the risk estimation aggregated to the partitions for one iteration
oneLemIter = function(
  Lambda, smoothingMat, regionMat, offsetMat, counts
){

  denom = t(regionMat) %*% offsetMat %*% Lambda

  em = t(t(counts/denom) %*% t(regionMat) %*% offsetMat) * Lambda
  em[as.vector(!is.finite(em))] = 0

  result = as.matrix(
    t(smoothingMat) %*% em
  )

  return(result)
}

# Computes the risk estimation of the training set
xvLemEst = function(trainId,
                    trainCounts,
                    regionMat,
                    offsetMat,
                    smoothingMat,
                    tol,
                    maxIter) {

  obsCounts = as.matrix(trainCounts[,trainId])

  oldRes = offsetMat %*%
    matrix(1,
           nrow = dim(offsetMat)[1],
           ncol = 1
    )

  Diter = 1
  absdiff = 1

  while((absdiff > tol) && (Diter < maxIter)) {

    newRes = oneLemIter(
      Lambda = oldRes,
      smoothingMat = smoothingMat,
      regionMat = regionMat,
      offsetMat = offsetMat,
      counts = obsCounts
    )

    absdiff = mean(abs(oldRes - newRes))

    oldRes = newRes
    Diter = Diter + 1
  }

  result = as.vector(newRes)

  return(result)
}

# Computes the entries of the smoothing matrix for all cells in the fine raster at the final iteration
# Specifically, the computations are done with the individual cells of the fine raster found in partitions
smoothingFinalStepEntries = function(
  cell1,
  cell2=cell1,
  offsetRaster,
  fine,
  coarse,
  focalList
) {
  if(is.list(focalList)) {
    # focal is not an array, it's a list of focal objects
    # compute the kernel array
    kernelArray = kernMat(
      cell2, cell1=cell1,
      focalList=focalList,
      coarse = coarse,
      fineRelativeRes = fine,reorder=FALSE)
  } else {
    # a kernel array was supplied
    kernelArray = focalList
  }

  if(is.function(kernelArray)) {
    Sbw = dimnames(kernelArray('q1x'))[[3]]
  } else {
    Sbw = dimnames(kernelArray)[[3]]
  }
  offsetBwNames = paste("offset.", Sbw, sep="")

  # extract data for extent of coarse cells from fine rasters
  cellInfo1=getCellInfo(
    cell1, coarse,
    fine, offsetRaster)

  if(!length(cellInfo1$idCoarse)){
    return(NULL)
  }

  cellsAreEqual = all(cell2==cell1)

  res = list()

  Scw2 = prod(res(fine))

  Spart1 = cellInfo1$idCoarse
  SNcells1 = tapply(cellInfo1$id[,"idCoarse"], list(cellInfo1$id[,"idCoarse"]), length)
  Sarea1 = SNcells1 * Scw2

  for(Dcell2 in cell2) {
    if(cellsAreEqual){
      cellInfo2=cellInfo1
      thisKernelArray = kernelArray
    } else {
      cellInfo2=getCellInfo(
        Dcell2,
        coarse,
        fine, offsetRaster)

      # deal with the reordering stuff
      if(is.function(kernelArray)){
        thisDist = xyFromCell(coarse, Dcell2) -
          xyFromCell(coarse, cell1)
        xoryaxis = c('y','x')[
          1+(abs(thisDist[1]) >= abs(thisDist[2]))
          ]
        quadrant = 1+
          (thisDist[2]<0) +
          (thisDist[1]<0) +
          2*(thisDist[1]<0 &thisDist[2]>=0)
        thisKernelArray = kernelArray(
          paste('q', quadrant, xoryaxis,sep='')
        )
      } else {
        thisKernelArray = kernelArray
      }
    } # end cells not equal

    Spart2 = cellInfo2$idCoarse
    SNcells2 = tapply(cellInfo2$id[,"idCoarse"], list(cellInfo2$id[,"idCoarse"]), length)
    Sarea2 = SNcells2 * Scw2

    if(!length(Spart2)){
      res[[as.character(Dcell2)]]=NULL
    }

    smoothedOffset2 = cellInfo2$offset[, offsetBwNames, drop=FALSE]
    smoothedOffset2 = array(smoothedOffset2,
                            c(dim(smoothedOffset2),cellInfo1$Ncells))
    smoothedOffset2 = aperm(smoothedOffset2, c(3,1,2))

    # the final smoothing matrix (or part thereof)
    Scell1 = ifelse(is.na(cellInfo2$id[,"idCoarse"]),
                    NA,
                    paste("c", Dcell2, "p", cellInfo2$id[,"idCoarse"], sep='')
    )
    Scell2 = ifelse(is.na(cellInfo1$id[,"idCoarse"]),
                    NA,
                    paste("c", cell1, "p", cellInfo1$id[,"idCoarse"], sep='')
    )


    theDim1 = c(
      cellInfo2$Ncells,
      length(Spart1),
      length(Sbw)
    )
    theDimnames1 = list(
      Scell1,
      paste("c", cell1, "p", Spart1, sep=''),
      Sbw)

    theDim2 = c(
      cellInfo1$Ncells,
      length(Spart2),
      length(Sbw)
    )
    theDimnames2 = list(
      Scell2,
      paste("c", Dcell2, "p", Spart2, sep=''),
      Sbw)


    pointsStepOne = complex(real = cellInfo2$points[,1], imaginary = cellInfo2$points[,2])
    names(pointsStepOne) = Scell1
    pointsStepRev = complex(real = cellInfo1$points[,1], imaginary = cellInfo1$points[,2])
    names(pointsStepRev) = Scell2


    smoothingStepOne = array(NA, dim=theDim1, dimnames=theDimnames1)
    smoothingStepRev = array(NA, dim=theDim2, dimnames=theDimnames2)

    for(D1 in 1:length(Spart1)){
      inD1 = which(cellInfo1$id[,'idCoarse']==Spart1[D1])
      smoothedOffset1 = cellInfo1$offset[, offsetBwNames, drop=FALSE]
      smoothedOffset1 = array(smoothedOffset1,
                              c(dim(smoothedOffset1),cellInfo2$Ncells))
      smoothedOffset1 = aperm(smoothedOffset1, c(1,3,2))

      for(D2 in 1:length(Spart2)){
        inD2 = which(cellInfo2$id[,'idCoarse']==Spart2[D2])
        numerator = thisKernelArray[inD1,inD2,,drop=FALSE]
        smoothingStepOne[inD2,D1,] =
          apply(numerator/smoothedOffset2[1:length(inD1),inD2,,drop=FALSE], c(2,3), sum) / Sarea1[D1]
        if(!cellsAreEqual) {
          smoothingStepRev[inD1,D2,] =
            apply(numerator/smoothedOffset1[inD1,1:length(inD2),,drop=FALSE], c(1,3), sum) / Sarea2[D2]
        } #cells Equal
      } #loop Spart22
    } # loop Spart1

    if(cellsAreEqual) {
      pointsStep = pointsStepOne
      smoothingStep = smoothingStepOne
    } else {
      pointsStep = list(pointsStepOne, pointsStepRev)
      names(pointsStep) = c('straightup', 'transpose')

      smoothingStep = list(smoothingStepOne, smoothingStepRev)
      names(smoothingStep) = c('straightup', 'transpose')
    }

    res[[as.character(Dcell2)]] = list(pointsStep, smoothingStep)
  } # loop cell2
  if(cellsAreEqual) {
    res = res[[1]]
  }

  res
}

# Computes all of the diagonal entries of the smoothing matrix for all cells in the fine raster at the final iteration
# Specifically, the computations are done with the cells in the same partitions
smoothingFinalStepDiag = function(
  rasterCoarse,rasterFine,
  focalList,offsetRaster, ncores
) {
  kernelArrayD = kernMat(
    cellDist=c(0,0),
    focalList=focalList,
    fineRelativeRes=rasterFine,
    cell1=c(0,0),
    coarse=rasterCoarse,
    reorder=FALSE
  )

  diagBlocks = parallel::mcmapply(
    smoothingFinalStepEntries,
    cell1 = 1:ncell(rasterCoarse),
    MoreArgs=list(
      focalList=kernelArrayD,
      coarse=rasterCoarse,
      fine=rasterFine,
      offsetRaster=offsetRaster
    ),
    mc.cores=ncores, SIMPLIFY=FALSE
  )

  cellsWithData = which(!unlist(lapply(diagBlocks, is.null)))
  diagBlocks = diagBlocks[cellsWithData]
  names(diagBlocks) = cellsWithData
  Spartitions = unlist(lapply(diagBlocks, function(qq) colnames(qq[[2]])))
  Sbw = unique(unlist(lapply(diagBlocks, function(qq) dimnames(qq[[2]])[[3]])))

  pointsStepArray = unlist(lapply(diagBlocks, function(qq) qq[[1]]))
  names(pointsStepArray) = gsub("[[:digit:]]+\\.", "", names(pointsStepArray))

  smoothingStepArray = array(NA,
                             c(length(pointsStepArray), length(Spartitions), length(Sbw)),
                             dimnames=list(names(pointsStepArray), Spartitions, Sbw)
  )

  for(D in 1:length(diagBlocks)) {
    inD = dimnames(diagBlocks[[D]][[2]])
    inP = match(diagBlocks[[D]][[1]], pointsStepArray)

    smoothingStepArray[inP,inD[[2]],] = diagBlocks[[D]][[2]]
  }

  # matrix of cell distances
  allCells = expand.grid(cellsWithData, cellsWithData)
  colnames(allCells) = gsub("^Var", "cell", colnames(allCells))
  allCells= cbind(allCells,
                  minCell = pmin(allCells[,'cell2'],allCells[,'cell1']),
                  maxCell = pmax(allCells[,'cell2'],allCells[,'cell1'])
  )
  allCells = cbind(allCells,
                   cellPair=1i*allCells[,'minCell'] + allCells[,'maxCell'])

  # row, col of cell1
  toAdd = rowColFromCell(rasterCoarse,allCells[,'cell1'])
  colnames(toAdd) = paste(colnames(toAdd), "1", sep="")
  allCells = cbind(allCells, toAdd)
  # row, col of cell2
  toAdd = rowColFromCell(rasterCoarse,allCells[,'cell2'])
  colnames(toAdd) = paste(colnames(toAdd), "2", sep="")
  allCells = cbind(allCells, toAdd)
  # row, col of min
  toAdd = rowColFromCell(rasterCoarse,allCells[,'minCell'])
  colnames(toAdd) = paste(colnames(toAdd), "min", sep="")
  allCells = cbind(allCells, toAdd)
  # row, col of max
  toAdd = rowColFromCell(rasterCoarse,allCells[,'maxCell'])
  colnames(toAdd) = paste(colnames(toAdd), "max", sep="")
  allCells = cbind(allCells, toAdd)

  # col dist
  allCells = cbind(allCells,
                   colDiff = allCells[,'col2'] -allCells[,'col1'],
                   rowDiff = allCells[,'row1'] -allCells[,'row2'],
                   colDiffUnique = allCells[,'colmax'] -allCells[,'colmin'],
                   rowDiffUnique = allCells[,'rowmin'] -allCells[,'rowmax']
  )

  # dist pair
  allCells = cbind(allCells,
                   dist=1i*allCells[,'rowDiff'] + allCells[,'colDiff'],
                   distUnique = 1i*allCells[,'rowDiffUnique'] + allCells[,'colDiffUnique']
  )

  # get rid of diagonals
  allCells = allCells[allCells[,'dist']!= 0,]
  allCells = allCells[allCells[,'dist']==allCells[,'distUnique'],]

  # unique distance in first quadrant, x dist bigger than y dists
  allCells = cbind(allCells,
                   distOf8 =
                     pmax(abs(Re(allCells[,'dist'])), abs(Re(allCells[,'dist']))) +
                     1i * pmin(abs(Im(allCells[,'dist'])), abs(Im(allCells[,'dist']))),
                   posY = (Im(allCells[,'dist'])>0)+1,
                   YgtX = (abs(Im(allCells[,'dist'])) > abs(Re(allCells[,'dist'])) )+1
  )

  # find the extent of the kernel with the biggest bandwidth
  maxDist = max(unlist(lapply(focalList$focal, dim)))
  # how many coarse cells it is
  maxCoarseCells = maxDist /
    max(ceiling(res(rasterCoarse)/res(rasterFine)))
  # find pairs of coarse cells separated by a distance bigger than this
  cellsAllZeros = pmax(abs(Re(allCells[,'dist'])),
                       abs(Im(allCells[,'dist']))) > maxCoarseCells

  allCells = allCells[!cellsAllZeros,]
  allCells = allCells[order(allCells[,'colDiff'], allCells[,'rowDiff']),]
  allCells = allCells[, c('cell1','cell2','dist','distOf8', 'posY','YgtX')]

  list(
    pointsStepArray=pointsStepArray,
    smoothingStepArray=smoothingStepArray,
    cells=allCells,
    uniqueDist=unique(allCells$dist)
  )
}

# Computes all of the off-diagonal entries of the smoothing matrix for all cells in the fine raster at the final iteration
# Specifically, the computations are done with the cells in neighbouring partitions
smoothingFinalStepOneDist = function(
  x, allCells,
  focal, coarse,
  fine, offsetRaster)
{
  kernelArray = kernMat(
    cellDist=c(Re(x), Im(x)),
    focalList=focal,
    coarse=coarse,
    fineRelativeRes=fine,
    reorder=TRUE)

  whichCells = allCells[which(allCells[,'dist']==x),]

  whichCellsRev = whichCells
  whichCellsRev$cellA = whichCells$cellB = whichCells$cell2
  whichCellsRev$cellB = whichCells$cellA = whichCells$cell1
  whichCellsTwice=rbind(whichCellsRev, whichCells)

  whichCellsD = NULL
  while(nrow(whichCellsTwice)){
    thisCell = as.integer(names(which.max(table(whichCellsTwice$cellA))))
    SthisCell = whichCellsTwice$cellA == thisCell
    whichCellsD = rbind(whichCellsD, whichCellsTwice[SthisCell,])
    whichCellsTwice = whichCellsTwice[!( SthisCell | (whichCellsTwice$cellB==thisCell)),]
  }

  res = list()
  for(Dcell1 in unique(whichCellsD[,'cell1'])){
    Scell2 = whichCellsD[whichCellsD[,'cell1']==Dcell1,'cell2']

    res[[as.character(Dcell1)]] =
      smoothingFinalStepEntries(
        cell1=Dcell1,
        cell2 = Scell2,
        focalList=kernelArray,
        coarse=coarse,
        fine=fine,
        offsetRaster=offsetRaster
      )
  }

  res
}

# Generates the smoothing matrix for all cells in the fine raster at the final iteration
smoothingFinalMat = function(
  lemObjects,
  bw,
  ncores
) {

  rasterCoarse = lemObjects$rasterCoarse
  rasterFine = lemObjects$rasterFine

  focalBw = lemObjects$focal$bw$bw %in% bw
  focal = lemObjects$focal[1:3]
  focal$focal = focal$focal[focalBw]
  focal$extended = focal$extended[focalBw]

  smoothingArray = lemObjects$smoothingArray[,,names(focal$focal)]
  if(length(dim(smoothingArray)) == 2) {
    smoothingArray = array(smoothingArray,
                           dim = c(dim(smoothingArray),1),
                           dimnames = list(dimnames(smoothingArray)[[1]], dimnames(smoothingArray)[[2]], names(focal$focal))
    )
  }

  rasterOffset = lemObjects$offset[[c("offset",
                                      paste("offset.", names(focal$focal), sep="")
  )]]

  theFinalStep = smoothingFinalStepDiag(
    rasterCoarse=rasterCoarse,
    rasterFine=rasterFine,
    focalList=focal,
    offsetRaster=rasterOffset,
    ncores=ncores)

  finalStepOffDiag = parallel::mcmapply(
    smoothingFinalStepOneDist,
    x=theFinalStep$uniqueDist,
    MoreArgs=list(
      allCells=theFinalStep$cells,
      focal=focal,
      coarse=rasterCoarse,
      fine=rasterFine,
      offsetRaster=rasterOffset
    ),
    mc.cores=ncores, SIMPLIFY=FALSE
  )

  #off-diagonals are positive
  if(length(finalStepOffDiag) > 0) {

    for(Ddist in 1:length(finalStepOffDiag)) {
      for(Dcell1 in 1:length(finalStepOffDiag[[Ddist]])) {
        for(Dcell2 in 1:length(finalStepOffDiag[[Ddist]][[Dcell1]])) {

          inDOne = dimnames(finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[2]][[1]])
          inPOne = match(finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[1]][[1]],
                         theFinalStep$pointsStepArray)

          inDRev = dimnames(finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[2]][[2]])
          inPRev = match(finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[1]][[2]],
                         theFinalStep$pointsStepArray)

          theFinalStep$smoothingStepArray[inPOne,inDOne[[2]],] =
            finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[2]][[1]]

          theFinalStep$smoothingStepArray[inPRev,inDRev[[2]],] =
            finalStepOffDiag[[Ddist]][[Dcell1]][[Dcell2]][[2]][[2]]
        }
      }
    }
  }

  #fill in the zeros in the smoothing matrix where distance is beyond largest bandwidth
  theFinalStep$smoothingStepArray[
    dimnames(theFinalStep$smoothingStepArray)[[1]] != "NA" &
      is.na(theFinalStep$smoothingStepArray)
    ] = 0

  result = list(
    pointsStepArray = theFinalStep$pointsStepArray,
    smoothingStepArray = theFinalStep$smoothingStepArray,
    smoothingArray = smoothingArray,
    regionMat = lemObjects$regionMat,
    offsetMat = lemObjects$offsetMat,
    offsetRaster = rasterOffset[["offset"]],
    rasterFine = rasterFine,
    polyCoarse = lemObjects$polyCoarse
  )

  return(result)
}

# Computes the relative risk estimation on the raster of fine polygons
riskEst = function(x,
                  lemObjects,
                  bw,
                  tol,
                  maxIter
) {

  regionMat = lemObjects$regionMat
  offsetMat = lemObjects$offsetMat
  smoothingMat = lemObjects$smoothingArray[,,bw]

  idCoarse = lemObjects$polyCoarse$id

  #fine raster did not include all regions in the coarse shapefile
  if(length(idCoarse) != dim(regionMat)[[2]]) {

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

  #risk estimation for point-locations
  offsetRaster = lemObjects$offsetRaster
  offsetRaster[values(offsetRaster) == 0] = NA

  pointsMat = cbind(Re(lemObjects$pointsStepArray), Im(lemObjects$pointsStepArray))
  smoothingStepMat = lemObjects$smoothingStepArray[,,bw]

  denom = t(regionMat) %*% offsetMat %*% Lambda
  em = t(t(obsCounts/denom) %*% t(regionMat) %*% offsetMat) * Lambda
  em[as.vector(!is.finite(em))] = 0
  em = as.matrix(em[match(dimnames(smoothingStepMat)[[2]], dimnames(em)[[1]])])

  pointsRisk = as.numeric(smoothingStepMat %*% em)

  pointsCells = !is.na(pointsRisk)
  pointsRisk = pointsRisk[pointsCells]
  pointsMat = pointsMat[pointsCells,]
  dimnames(pointsMat)[[1]] = 1:dim(pointsMat)[1]

  estRisk = SpatialPointsDataFrame(
    coords = pointsMat,
    data = data.frame(risk = pointsRisk),
    proj4string = CRS(proj4string(offsetRaster))
  )

  result = rasterize(estRisk, offsetRaster, field = "risk")
  result = raster::mask(result, offsetRaster)
  result = raster::trim(result)

  return(result)
}
