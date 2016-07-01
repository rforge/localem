# Retrieves the cell information of the rasters of the coarse and fine polygons, and offsets
getCellInfo = function(
  cell1, coarse,
  fine, offsetRaster
) {

  lims=drop(xyFromCell(coarse,cell1))
  resHalf = raster::res(coarse)/2

  theextent = extent(c(
    lims['x']+c(-1,1)*resHalf[1],
    lims['y']+c(-1,1)*resHalf[2]
  ))

  Scells = cellsFromExtent(fine,theextent)

  Scols = sort(unique(
    colFromCell(fine, Scells)
  ))
  Mcol = min(Scols)
  Ncols = length(Scols)
  Srows = sort(unique(
    rowFromCell(fine, Scells)
  ))
  Mrow = min(Srows)
  Nrows = length(Srows)

  id = getValuesBlock(fine,
                      row=Mrow, nrows=Nrows,col=Mcol,ncols=Ncols)

  theoffset = getValuesBlock(offsetRaster,
                             row=Mrow, nrows=Nrows,col=Mcol,ncols=Ncols)

  list(
    extent = theextent,
    Scells=Scells,
    Mcol=Mcol,Ncols=Ncols,Mrow=Mrow,Nrows=Nrows,
    Ncells = Ncols*Nrows,
    id=id,offset=theoffset,
    idCoarse=sort(unique(id[,'idCoarse'])),
    points = xyFromCell(fine, Scells)
  )
}

# Computes the entries of the smoothing matrix for the partitions
smoothingMatrixEntries = function(
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
        xoryaxis = c('y','x')[1 +
                                (abs(thisDist[1]) >= abs(thisDist[2]))
                              ]
        quadrant = 1+
          (thisDist[2]<0) +
          (thisDist[1]<0) +
          2*(thisDist[1]<0 & thisDist[2]>=0)
        thisKernelArray = kernelArray(
          paste('q', quadrant, xoryaxis, sep='')
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

    # the smoothing matrix (or part thereof)
    theDimnames = list(
      paste("c", cell1, "p", cellInfo1$idCoarse, sep=''),
      paste("c", Dcell2, "p", cellInfo2$idCoarse, sep=''),
      Sbw,
      c('straightup', 'transpose')
    )

    theDim = c(
      length(cellInfo1$idCoarse),
      length(cellInfo2$idCoarse),
      length(Sbw),2)

    smoothingMatrix = array(NA, theDim,
                            dimnames=theDimnames
    )
    smoothingMatrixOne = array(NA, theDim[-4],
                               dimnames=theDimnames[-4]
    )

    for(D1 in 1:length(Spart1)){
      inD1 = which(cellInfo1$id[,'idCoarse']==Spart1[D1])
      smoothedOffset1 = cellInfo1$offset[, offsetBwNames, drop=FALSE]
      smoothedOffset1 = array(smoothedOffset1,
                              c(dim(smoothedOffset1),cellInfo2$Ncells))
      smoothedOffset1 = aperm(smoothedOffset1,c(1,3,2))

      for(D2 in 1:length(Spart2)){
        inD2 = which(cellInfo2$id[,'idCoarse']==Spart2[D2])
        numerator = thisKernelArray[inD1, inD2,,drop=FALSE]
        smoothingMatrixOne[D1,D2,] =
          apply(numerator/smoothedOffset2[1:length(inD1),inD2,,drop=FALSE], 3, sum) * Scw2 / Sarea1[D1]
        if(!cellsAreEqual){
          smoothingMatrix[D1,D2,,'transpose'] =
            apply(numerator/smoothedOffset1[inD1,1:length(inD2),,drop=FALSE], 3, sum) * Scw2 / Sarea2[D2]
        } #cells Equal
      } #loop Spart22
    } # loop Spart1

    if(cellsAreEqual) {
      smoothingMatrix = smoothingMatrixOne
    } else {
      smoothingMatrix[,,,'straightup'] = smoothingMatrixOne
    }
    res[[as.character(Dcell2)]] = smoothingMatrix
  } # loop cell2

  if(cellsAreEqual) {
    res = res[[1]]
  }
  res
}

# Computes all of the diagonal entries of the smoothing matrix for the partitions
# Specifically, the computations are done for the same partitions
# Generates also the matrix that matches the IDs of the coarse polygons with the partitions
smoothingMatrixDiag = function(
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
    smoothingMatrixEntries,
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
  Spartitions = unlist(lapply(diagBlocks, function(qq) dimnames(qq)[[1]]))
  Sbw = dimnames(diagBlocks[[1]])[[3]]
  smoothingArray = array(NA,
                         c(length(Spartitions), length(Spartitions), length(Sbw)),
                         dimnames=list(Spartitions,Spartitions, Sbw)
  )

  for(D in 1:length(diagBlocks)){
    partHere = dimnames(diagBlocks[[D]])[[1]]
    smoothingArray[partHere,partHere,] =
      diagBlocks[[D]]
  }

  # matrices for local-EM M bit
  Scells = ifelse(is.na(values(rasterFine[["cellCoarse"]])) | is.na(values(rasterFine[["idCoarse"]])),
                  NA,
                  paste("c", values(rasterFine[["cellCoarse"]]), "p", values(rasterFine[["idCoarse"]]), sep = "")
  )
  meanOffsets = tapply(values(offsetRaster[['offset']]),
                       list(Scells),
                       mean, na.rm=TRUE)
  meanOffsets = meanOffsets[match(Spartitions, names(meanOffsets))]
  offsetMat = Diagonal(length(Spartitions), meanOffsets)

  regions = gsub("^c[[:digit:]]+p", "", Spartitions)
  regions = as.integer(regions)
  regionMat = outer(regions, regions, '==')
  regionMat = Matrix(regionMat)

  dimnames(regionMat)=dimnames(offsetMat) =
    list(Spartitions,Spartitions)

  expandCountMat = regionMat
  expandCountMat = expandCountMat[,!duplicated(regions),drop=FALSE]
  colnames(expandCountMat) = gsub("^c[[:digit:]]+p", "", colnames(expandCountMat) )
  expandCountMat = expandCountMat[,
                                  as.character(sort(as.integer(colnames(expandCountMat)))),
                                  drop=FALSE
                                  ]

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

  allCells =allCells[!cellsAllZeros,]
  allCells = allCells[order(allCells[,'colDiff'], allCells[,'rowDiff']),]
  allCells = allCells[, c('cell1','cell2','dist','distOf8', 'posY','YgtX')]

  list(smoothingArray=smoothingArray,
       regionMat=expandCountMat,
       offsetMat = offsetMat,
       cells=allCells,
       uniqueDist = unique(allCells$dist))
}

# Computes all of the off-diagonal entries of the smoothing matrix for the partitions
# Specifically, the computations are done for neighbouring partitions
smoothingMatrixOneDist = function(
  x, allCells,
  focal, coarse,
  fine, offsetRaster
) {
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
    whichCellsTwice = whichCellsTwice[!( SthisCell | (whichCellsTwice$cellB==thisCell) ),]
  }


  res = list()
  for(Dcell1 in unique(whichCellsD[,'cell1'])){
    Scell2 = whichCellsD[
      whichCellsD[,'cell1']==Dcell1,
      'cell2'
      ]

    res[[as.character(Dcell1)]] =
      smoothingMatrixEntries(
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
