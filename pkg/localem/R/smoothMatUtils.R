# Retrieves the cell information of the rasters of the coarse and fine polygons, and offsets
getCellInfo = function(
    cell1, coarse,
    fine, offsetRaster
) {
  
  lims=drop(raster::xyFromCell(coarse,cell1))
  resHalf = raster::res(coarse)/2
  
  theextent = raster::extent(c(
          lims['x']+c(-1,1)*resHalf[1],
          lims['y']+c(-1,1)*resHalf[2]
      ))
  
  Scells = raster::cellsFromExtent(fine,theextent)
  
  Scols = sort(unique(
          raster::colFromCell(fine, Scells)
      ))
  Mcol = min(Scols)
  Ncols = length(Scols)
  Srows = sort(unique(
          raster::rowFromCell(fine, Scells)
      ))
  Mrow = min(Srows)
  Nrows = length(Srows)
  
  id = raster::getValuesBlock(fine,
      row=Mrow, nrows=Nrows,
      col=Mcol,ncols=Ncols)
  id = cbind(
      raster::factorValues(fine, id), 
      index = 1:length(id)
  )
  
  theoffset = raster::getValuesBlock(offsetRaster,
      row=Mrow, nrows=Nrows,col=Mcol,ncols=Ncols)
  
  
  partitionDf = as.data.frame(table(
          fine=id[,'idFine'], 
          coarse=id[,'idCoarse'])
  )
  partitionDf = partitionDf[partitionDf$Freq > 0, , drop=FALSE]
  
  list(
      extent = theextent,
      Scells=Scells,
      Mcol=Mcol,Ncols=Ncols,Mrow=Mrow,Nrows=Nrows,
      Ncells = Ncols*Nrows,
      id=id,offset=theoffset,
      idCoarse=sort(unique(id[,'idCoarse'])),
      partitions = partitionDf,
      points = raster::xyFromCell(fine, Scells)
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
  
  
  offsetBwNames = grep("^bw([[:digit:]]|[.])+(xv[[:digit:]]+)?$", names(offsetRaster), value=TRUE)
  # extract data for extent of coarse cells from fine rasters
  cellInfo1=getCellInfo(
      cell1, coarse,
      fine, offsetRaster)
  
  if(!length(cellInfo1$idCoarse)){
    return(NULL)
  }
  
  cellsAreEqual = all(cell2==cell1)
  
  
  res = list()
  
  Scw2 = prod(raster::res(fine))
  
#  Spart1 = cellInfo1$idCoarse
#  SNcells1 = tapply(cellInfo1$id[,"idCoarse"], list(cellInfo1$id[,"idCoarse"]), length)
  Sarea1 = cellInfo1$partition$Freq * Scw2
  
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
        thisDist = raster::xyFromCell(coarse, Dcell2) -
            raster::xyFromCell(coarse, cell1)
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
    
    #Spart2 = cellInfo2$idCoarse
    #SNcells2 = tapply(cellInfo2$id[,"idCoarse"], list(cellInfo2$id[,"idCoarse"]), length)
    #Sarea2 = SNcells2 * Scw2
    Sarea2 = cellInfo2$partition$Freq * Scw2
    
    if(!length(Sarea2)){
      res[[as.character(Dcell2)]]=NULL
    }
    
    smoothedOffset2 = cellInfo2$offset[, offsetBwNames, drop=FALSE]
    smoothedOffset2 = array(smoothedOffset2,
        c(dim(smoothedOffset2),cellInfo1$Ncells),
        dimnames = c(dimnames(smoothedOffset2), 
            list(1:cellInfo1$Ncells)))
    smoothedOffset2 = aperm(smoothedOffset2, c(3,1,2))
    
    # the smoothing matrix (or part thereof)
    theDimnames = list(
        paste("c", cell1, "p", 
            cellInfo1$partition$coarse, '.', 
            cellInfo1$partition$fine, sep=''),
        paste("c", Dcell2, "p",						
            cellInfo2$partition$coarse, '.', 
            cellInfo2$partition$fine, sep=''),
        dimnames(smoothedOffset2)[[3]],
        c('straightup', 'transpose')
    )
    
    
    smoothingMatrix = array(NA, 
        unlist(lapply(theDimnames, length)),
        dimnames=theDimnames
    )
    smoothingMatrixOne = array(NA, 
        dim(smoothingMatrix)[-4],
        dimnames=theDimnames[-4]
    )
    
    for(D1 in 1:dim(smoothingMatrix)[1]){
      inD1 = which(
          (cellInfo1$id[,'idCoarse']==cellInfo1$partition[D1,'coarse']) &
              (cellInfo1$id[,'idFine']==cellInfo1$partition[D1,'fine'])
      )
      smoothedOffset1 = cellInfo1$offset[, offsetBwNames, drop=FALSE]
      smoothedOffset1 = array(smoothedOffset1,
          c(dim(smoothedOffset1),cellInfo2$Ncells))
      smoothedOffset1 = aperm(smoothedOffset1,c(1,3,2))
      
      for(D2 in 1:dim(smoothingMatrix)[2]){
        inD2 = which(
            (cellInfo2$id[,'idCoarse']==cellInfo2$partition[D2,'coarse']) &
                (cellInfo2$id[,'idFine']==cellInfo2$partition[D2,'fine'])
        )
        
        numerator = thisKernelArray[inD1, inD2,,drop=FALSE]
        smoothingMatrixOne[D1,D2,] = apply(
            numerator/smoothedOffset2[1:length(inD1),inD2,,drop=FALSE], 
            3, sum) * Scw2 / Sarea1[D1]
        if(!cellsAreEqual){
          smoothingMatrix[D1,D2,,'transpose'] =	apply(
              numerator/smoothedOffset1[inD1,1:length(inD2),,drop=FALSE], 
              3, sum) * Scw2 / Sarea2[D2]
        } #cells Equal
      } #loop Spart22
    } # loop Spart1
    
    if(cellsAreEqual) {
      smoothingMatrix = smoothingMatrixOne
    } else {
      smoothingMatrix[,,,'straightup'] = smoothingMatrixOne
    }
    res[[as.character(Dcell2)]] = smoothingMatrix
    # TO DO: write values into smoothing raster with spatial.tools::binary.image.write
  } # loop cell2
  
  if(cellsAreEqual) {
    res = res[[1]]
  }
  res
}


# function for clusterMap
oneBlockFun = function(Dcell1,     
    kernelArrayD, 
    rasterCoarse, 
    rasterFine, 
    offsetRaster,
    Spartitions, 
    layerSeq, 
    image_x, 
    image_y, 
    image_z,
    smoothingRaster, 
    theType,
    verbose
) {
  thisBlock = try(
      smoothingMatrixEntries(cell1 = Dcell1,
          focalList=kernelArrayD,
          coarse=rasterCoarse,
          fine=rasterFine,
          offsetRaster=offsetRaster
      ))
  if(class(thisBlock) == 'try-error') {
    warning(paste("Error in smoothingMatrixEntries", Dcell1))
    thisBlock = NULL
  }
  
  if(!is.null(thisBlock)) {
    
    partHere = dimnames(thisBlock)[[1]]
    matchPartHere = match(partHere, Spartitions)
    
    theOrder = order(matchPartHere)
    # from spatial.tools::binary_image_write          
    data_position = t(expand.grid(
            matchPartHere[theOrder], 
            matchPartHere[theOrder], 
            layerSeq
        ))
# turns out cell_position doesn't need to be an integer
    
    cell_position=
        ( (data_position[2,]-1)*image_x)+
        (data_position[1,])+
        ((data_position[3,]-1)*(image_x*image_y))
    
    haveWritten = FALSE
    writeCounter1 = 0
    
    while( (!haveWritten) & (writeCounter1 < 20) ) {
      
      out = mmap::mmap(
          smoothingRaster,
          mode=theType
      )	
      
      
      if(TRUE) {
        .Call("mmapReplaceReal", as.double(cell_position), 
            as.double(thisBlock[theOrder,theOrder,]), 
            out)#, PACKAGE='localEM') 
      } else {
        out[cell_position] = as.numeric(thisBlock[theOrder,theOrder,])        
      }
      
      haveWrittenFirst = mmap::munmap(out)  
      haveWritten = identical(haveWrittenFirst, 0L)
      if(!haveWritten & verbose){
        warning(paste("error in write attempt", writeCounter1, "of cell", Dcell1, haveWrittenFirst))
      }
      writeCounter1 = writeCounter1 + 1
    } # end while 
  } # end not null
  dim(thisBlock)
} # end oneBlockFun


oneBlockOffdiagFun = function(
    x,
    theMat, 
    rasterObjects,
    image_x,
    image_y,
    image_z,
    Spartitions,
    theType,
    smoothingRasterFile,
    layerSeq,
    verbose
) {
  
  thisBlock = smoothingMatrixOneDist(
      x=x,
      allCells=theMat$cells,
      focal=rasterObjects$focal,
      coarse=rasterObjects$rasterCoarse,
      fine=rasterObjects$rasterFine,
      offsetRaster=rasterObjects$offset
  )
  if(class(thisBlock) != 'try-error') {
    for(Dcell1 in 1:length(thisBlock)) {
      for(Dcell2 in 1:length(thisBlock[[Dcell1]])) {
        
        partHere = thisBlock[[Dcell1]][[Dcell2]]
        s1 = dimnames(thisBlock[[Dcell1]][[Dcell2]])[[1]]
        s2 = dimnames(thisBlock[[Dcell1]][[Dcell2]])[[2]]
        partHere = try(aperm(thisBlock[[Dcell1]][[Dcell2]][,,,'transpose', drop=FALSE], 
                c(2,1,3,4)))
        matchPartHere = lapply(dimnames(partHere)[1:2], match, 
            table=Spartitions)
        
        theOrder = lapply(matchPartHere, order)
        
        
        data_position = t(expand.grid(
                as.vector(matchPartHere[[1]][theOrder[[1]] ]), 
                as.vector(matchPartHere[[2]][theOrder[[2]] ]), 
                layerSeq)
        )
        
        cell_position=
            ( (data_position[2,]-1)*image_x)+
            (data_position[1,])+
            ((data_position[3,]-1)*(image_x*image_y))
        
        
        haveWritten = FALSE
        writeCounter1 = 0
        
        while(!haveWritten & (writeCounter1 < 20)) {
          out = mmap::mmap(
              smoothingRasterFile,
              mode=theType
          )	
          
          if(TRUE) {
            .Call("mmapReplaceReal", as.double(cell_position), 
                as.double(partHere[theOrder[[1]], theOrder[[2]],,] ),
                out)#, PACKAGE='localEM') 
          } else {
            out[cell_position] = as.double(partHere[theOrder[[1]], theOrder[[2]],,] )      
          }
          
          
          haveWrittenFirst = mmap::munmap(out)  
          
          haveWritten = identical(haveWrittenFirst, 0L)
          writeCounter1 = writeCounter1 + 1
        }
        if(writeCounter1 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
        
        partHere = try(thisBlock[[Dcell1]][[Dcell2]][,,,'straightup', drop=FALSE])
        matchPartHere = lapply(dimnames(partHere)[1:2], match, 
            table=Spartitions)
        theOrder = lapply(matchPartHere, order)
        
        data_position = t(expand.grid(
                as.vector(matchPartHere[[1]][theOrder[[1]] ]), 
                as.vector(matchPartHere[[2]][theOrder[[2]] ]), 
                layerSeq)
        )
        
        cell_position=
            ( (data_position[2,]-1)*image_x)+
            (data_position[1,])+
            ((data_position[3,]-1)*(image_x*image_y))
        
        
        haveWritten = FALSE
        writeCounter2 = 0
        
        while(!haveWritten & (writeCounter2 < 20)) {
          
          out = mmap::mmap(
              smoothingRasterFile,
              mode=theType
          )	
          dim(out)
          
          if(TRUE) {
            .Call("mmapReplaceReal", as.double(cell_position), 
                as.double(partHere[theOrder[[1]], theOrder[[2]],,] ),
                out)#, PACKAGE='localEM') 
          } else {
            out[cell_position] = as.double(partHere[theOrder[[1]], theOrder[[2]],,] )      
          }
          
          haveWrittenFirst = mmap::munmap(out)  
          
          haveWritten = identical(haveWrittenFirst, 0L)
          writeCounter2 = writeCounter2 + 1
        }
        if(writeCounter2 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
      } # end Dcell2  
    } # end Dcell1
    res = c(writeCounter1, writeCounter1)
  }  else { # end try error
    res =thisBlock
  }
  res
} # oneBlockOffdiag


# Computes all of the diagonal entries of the smoothing matrix for the partitions
# Specifically, the computations are done for the same partitions
# Generates also the matrix that matches the IDs of the coarse polygons with the partitions
smoothingMatrixDiag = function(
    rasterCoarse,rasterFine,
    focalList,offsetRaster, 
    filename,
    cl = NULL,
    verbose=FALSE
) {
  
  kernelArrayD = kernMat(
      cellDist=c(0,0),
      focalList=focalList,
      fineRelativeRes=rasterFine,
      cell1=c(0,0),
      coarse=rasterCoarse,
      reorder=FALSE,
      xv = names(offsetRaster))
  
  
  # order the partitions
  partitionIdMat = raster::levels(rasterFine)[[1]][,c('cellCoarse','idCoarse','partition')]
  partitionOrder = order(partitionIdMat[,1], partitionIdMat[,2],partitionIdMat[,3])
  Spartitions = partitionIdMat[partitionOrder,'partition']
  
  Npartitions = length(Spartitions)
  Nsmooths = dim(kernelArrayD)[3]
  layerSeq = 1:Nsmooths
  Sbw = dimnames(kernelArrayD)[[3]]
  
  smoothingRasterTemplate = brick(
      nrows=Npartitions, ncols=Npartitions, 
      xmn=0,xmx=1,ymn=0,ymx=1,
      nl = Nsmooths
  )
  names(smoothingRasterTemplate) = Sbw		
  
  firstFile = filename
  
  # write zeros in the smoothing matrix
  toWrite = matrix(as.double(0), nrow(smoothingRasterTemplate), nlayers(smoothingRasterTemplate))
  
  if(verbose) cat("creating raster to hold smoothing matrix\n")
  
  firstRaster = writeStart(
      smoothingRasterTemplate,
      filename = firstFile, 
      format = 'raster',
      datatype = 'FLT8S', 
      bandorder = 'BSQ',
      overwrite=file.exists(firstFile))
  
  for(Drow in 1:nrow(firstRaster)) {
    firstRaster = writeValues(firstRaster, toWrite, Drow) 
  }  
  firstRaster = writeStop(firstRaster)
  smoothingRaster = gsub("grd$", "gri", filename(firstRaster))
  
  # theType = mmap::real64();dput(theType, '')
  theType = structure(numeric(0), bytes = 8L, signed = 1L, class = c("Ctype", 
          "double"))
  
  
  if(verbose) cat("looping through diagonal cells\n")
  
  # dimensions, from spatial.tools::binary_image_write
  image_dims = dim(smoothingRasterTemplate)
  image_x=image_dims[1]
  image_y=image_dims[2]
  image_z=image_dims[3]
  
  forMoreArgs = list(
      kernelArrayD = kernelArrayD, 
      rasterCoarse = rasterCoarse, 
      rasterFine = rasterFine, 
      offsetRaster =offsetRaster,
      Spartitions =Spartitions, 
      layerSeq =layerSeq, 
      image_x =image_x, 
      image_y =image_y, 
      image_z =image_z,
      smoothingRaster = smoothingRaster, 
      theType = theType,
      verbose=verbose
  )
  
  if(!is.null(cl)) {
    
    parallel::clusterExport(cl, 
        varlist = c('smoothingMatrixEntries',
            'kernMat',
            'reorderCellsTranslate',
            'getCellInfo'), 
        envir = environment()
    )
    
    diagBlocks = parallel::clusterMap(
        cl, oneBlockFun, 
        Dcell1 = 1:ncell(rasterCoarse), 
        MoreArgs = forMoreArgs)
    
  } else {
    diagBlocks = mapply(
        oneBlockFun,
        Dcell1 = 1:ncell(rasterCoarse), 
        MoreArgs = forMoreArgs)
  }
  
  cellsWithData = which(unlist(lapply(diagBlocks, length))>0)
  
  
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
  
  outArray = brick(smoothingRaster)
  names(outArray) = Sbw
  list(
      smoothingArray = outArray,
      partitions = Spartitions,
      bw = Sbw,
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
      reorder=TRUE,
      xv = names(offsetRaster))
  
  
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

