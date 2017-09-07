# Sets the size of sigma for Gaussian density kernel
# The default is 3 times sigma
focalWeightWithSize = function(x, bw, size=NULL) {
  
  if(!length(size))
    size = 3*bw
  
  raster::focalWeight(
    x=x, d=c(bw, size),
    type='Gauss'
  )
  
}

bigListFun = function(Dbw, focalList, bigCentreCell, bigMat) {
  theDim = dim(focalList[[Dbw]])
  topleft = bigCentreCell - (theDim-1)/2
  xseq = seq(topleft[1], len=theDim[1], by=1)
  xNotZero = (xseq > 0) & (xseq <= dim(bigMat)[1])
  yseq = seq(topleft[2], len=theDim[2], by=1)
  yNotZero = (yseq > 0) & (yseq <= dim(bigMat)[2])
  res = bigMat
  res[
    xseq[xNotZero], yseq[yNotZero]
  ] = focalList[[Dbw]][xNotZero,yNotZero]
  Matrix::Matrix(res)
}


# Computes focal weight matrix for a Gaussian density kernel with specified bandwidth and size of sigma
focalFromBw = function(
  bw, fine, minDim = Inf, focalSize = NULL, cl=NULL
){
  centreCell = 1+dim(fine)[1:2]
  
  theDim = 2*rep(max(centreCell),2)-1
  bigCentreCell = (theDim-1)/2
  
  bigMat = matrix(0,
    theDim[1],theDim[2]
  )
  names(bw) = gsub("[[:space:]]", "",
    format(bw, scientific=FALSE)
  )
  
  if(is.null(focalSize)) {
    focalSize = 3*max(bw)
  }
  
# define Dbw to make package check happy
  
  forMoreArgs = list(x=fine, size=focalSize)
  if(!is.null(cl)) {
    focalList = parallel::clusterMap(
      cl, focalWeightWithSize, 
      bw=bw,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  } else {
    focalList = mapply(
      focalWeightWithSize, 
      bw=bw,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  }

  
  names(focalList) = paste("bw",
    names(bw),
    sep="")

  
  forMoreArgs = list(
    focalList=focalList, bigCentreCell=bigCentreCell, bigMat=bigMat
  )

  
  
  if(!is.null(cl)) {
    bigList = parallel::clusterMap(
      cl, bigListFun, 
      Dbw=names(focalList),
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  } else {
    bigList = mapply(
      bigListFun, 
      Dbw=names(focalList),
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  }

  
  
  bw = data.frame(
    bw=bw,
    dim=unlist(lapply(focalList,
        function(qq) min(dim(qq))))
  )
  
  bw$factLog2 = floor(pmax(0,log2(bw$dim) - log2(minDim)))
  bw$fact = 2^bw$factLog2
  
  result = list(focal=focalList,
    extended=bigList,centreCell=bigCentreCell,
    bw=bw)
  
  for(D in sort(setdiff(unique(bw$fact),1))){
    fineAgg = raster::aggregate(fine, fact=D)
    bwHere = sort(bw[bw$fact==D,'bw'])
    focalListD =
      mapply(
        focalWeightWithSize,
        bw=bwHere,
        MoreArgs=list(
          x=fineAgg,
          size = focalSize
        )
      )
    names(focalListD) = paste("bw",
      bwHere,
      sep="")
    result[[paste('focalAgg', D,sep='')]] = focalListD
  }
  
  return(result)
}


# Reorders the focal weight matrix to estimate kernel smoothing function for specified bandwidth
# Utilises the symmetrical design of the focal weight matrix
reorderCellsTranslate = function(Ncell){
  
  Ncell = Ncell[1]
  NcellSq = Ncell^2
  cellSeq = 1:NcellSq
  cellMat = matrix(cellSeq, Ncell,Ncell)
  reorderRows = as.vector(cellMat[,Ncell:1])
  reorderBoth = rev(cellSeq)
  switchRowCol = as.vector(t(cellMat))
  
  list(
    q2=cellSeq[switchRowCol],
    q3=reorderRows[switchRowCol],
    q4=reorderRows
  )
}


# Computes the kernel smoothing function for specified bandwidth
kernMat = function(
  cellDist, focalList,
  fineRelativeRes,
  cell1=c(0,0),
  coarse=NULL,
  reorder=FALSE,
  xv = NULL
) {
  
  
  if(length(cell1)==2 & length(cellDist)==2){
    cellDist = cellDist - cell1
  } else if(length(cellDist)== 1 & length(cell1) == 1) {
    if(is.null(coarse)) {
      warning("supply coarse if cellDist is a cell ID")
    }
    cellCoord = rowColFromCell(coarse, c(cell1,cellDist))
    cellDist = apply(cellCoord,2,diff)
  } else {
    warning('cellDist and cel1 should be both either length 1 (cell id) or length 2 (x,y)')
  }
  
  if(reorder){
    cellDist = c(max(abs(cellDist)), min(abs(cellDist)))
  }
  
  Sbw = names(focalList$extended) # names of bandwidths
  
  if(!is.vector(fineRelativeRes)) {
    if(is.null(coarse)) {
      warning("supply coarse if fineRelativeRes is a raster")
    }
    fineRelativeRes = round(raster::res(coarse)/raster::res(fineRelativeRes))
  }
  cellDistFine = cellDist * fineRelativeRes
  
  # construct array to hold result
  Ncell = round(prod(fineRelativeRes))
  kernelArray = array(NA,
    c(Ncell, Ncell, length(Sbw)),
    dimnames=list(NULL, NULL, Sbw)
  )
  
  # loop through fine cells in coarse cell 1
  
  for(Dcol in 1:(fineRelativeRes[2])){
    firstCellInCol = fineRelativeRes[2]*(Dcol-1)
    
    Scol = seq(
      focalList$centreCell[2]+cellDistFine[2]-Dcol+1,
      len= fineRelativeRes[2],by=1)[1:fineRelativeRes[2]]
    
    for(Drow in 1:(fineRelativeRes[1])){
      Dcell = firstCellInCol + Drow
      Srow = seq(
        focalList$centreCell[1]+cellDistFine[1]-Drow+1,
        len= fineRelativeRes[1],by=1)[1:fineRelativeRes[1]]
      for(Dbw in Sbw){
        kernelArray[Dcell,,Dbw] = as.vector(
          focalList$extended[[Dbw]][Srow, Scol])
      } # for Dbw
    } #  for Drow
  } # for Dcol
  
  if(length(xv)) {
    offsetBwNames = grep("^bw([[:digit:]]|[.])+(xv[[:digit:]]+)?$", xv, value=TRUE)
    Sbw = gsub("xv[[:digit:]]+", "", offsetBwNames)
    kernelArray = kernelArray[,,Sbw]
    dimnames(kernelArray)[[3]] = offsetBwNames
  } 
  
  
  
  if(reorder){
    resEnv = new.env()
    with(resEnv, kernelArray <- kernelArray)
    
    with(resEnv,
      reorder <- reorderCellsTranslate(fineRelativeRes)
    )
    
    resFun = function(qx) {
      if(qx=='q1x') { #x>0, y=0
        resk = kernelArray
      } else if(qx=='q2x') { #|x|>|y| or |x|=|y|; x>0,y<0
        resk = kernelArray
      } else if(qx=='q2y') { #|y|>|x|; x>0 or x=0, y<0
        resk = kernelArray[reorder$q2,reorder$q2,,drop=FALSE]
      } else if(qx == 'q3y') { #|y|>|x|; x<0 or x=0, y<0
        resk = kernelArray[reorder$q3,reorder$q3,,drop=FALSE]
      } else if(qx == 'q3x') { #|x|>|y| or |x|=|y|; x<0 , y<0
        resk = aperm(kernelArray[reorder$q4,reorder$q4,,drop=FALSE],c(2,1,3))
      } else {
        warning('qx cant be ', qx)
        resk = NULL
      }
      resk
    }
    
    environment(resFun) = resEnv
    resHere = resFun
    
  } else {
    resHere = kernelArray
  }
  resHere
}

# Computes the smoothed offsets as a RasterLayer
smoothedOffsetMapFun = function(
  x=x, offsetTempFileAgg, 
  focalArray, 
  forSmooth, Soutfile
) {
  offsetHere = raster::brick(offsetTempFileAgg)[[ forSmooth[x,'layer'] ]]
  res = raster::focal(
    offsetHere,
    w = focalArray[,,forSmooth[x,'bw'] ],
    na.rm=TRUE, pad=TRUE,
    filename = Soutfile[x],
    overwrite = file.exists(Soutfile[x])
  )
  filename(res)
}


# Implements the focal matrix (with edge correction) to raster
focalFunction = function(x, focalArray, Scvsets)  {
  apply(focalArray*x[,,Scvsets,drop=FALSE], 
    3, sum, na.rm=TRUE)
}

focalMultOneRowOuter = function(x, focalMat, focalDim, out, row, mmap = FALSE) {
  valuesHere = raster::getValuesFocal(x, row, 1L, focalDim, array=FALSE, padValue = NA)
  valuesNa = is.na(valuesHere[[1]])
  
  toWrite = simplify2array(lapply(valuesHere, function(xx) {
        xx[valuesNa] = 0
        xx %*% focalMat
      }
    ))
  toWriteOnes =  (!valuesNa) %*% focalMat
  
  naPoints = is.na(raster::getValues(x[[1]], row, 1L))
  toWriteOnes[naPoints & (toWriteOnes < 0.5)] = NA
  
  toWrite = abind::abind(
    toWrite, 
    ones= toWriteOnes,
    along=3
  )  
  toWrite = matrix(
    toWrite, 
    nrow = dim(toWrite)[1],
    ncol = prod(dim(toWrite)[-1])
  )
  
  if(mmap) {
    
    outMap = mmap::mmap(
      gsub("[.]grd$", ".gri", raster::filename(out)),
      mode=structure(numeric(0), bytes = 8L, signed = 1L, class = c("Ctype", 
          "double"))
    )	
    
    # band ordering BSQ
    indexSeq = mapply(seq, 
      from=seq(
        from= raster::cellFromRowCol(out, row, 1), 
        len=raster::nlayers(out), by=raster::ncell(out)),
      MoreArgs = list(
        by = 1, len = raster::ncol(out)
      )
    )
    
    outMap[indexSeq] = as.vector(toWrite)
    
    res = mmap::munmap(outMap)  
    
  } else {
    
    raster::writeValues(out, toWrite, row)
    res = TRUE
  }
  res
}


focalMultOneRow = function(x, focalMat, focalDim, out, row, bw, uniqueBw=NULL, mmap = FALSE) {
  
  valuesHere = raster::getValuesFocal(x, row, 1L, focalDim, array=FALSE, padValue = NA)
  if(is.matrix(valuesHere)) valuesHere = list(valuesHere)
  valuesNa = is.na(valuesHere[[1]])
  colnames(valuesNa) = colnames(focalDim)
  
  toWrite = list()
  for(D in 1:length(valuesHere)) {
    valuesHere[[D]][valuesNa] = 0
    toWrite[[D]] = valuesHere[[D]] %*% focalMat[, bw[D],drop=FALSE ]
  }
  toWrite = do.call(cbind, toWrite)
  colnames(toWrite) = paste(colnames(toWrite), names(valuesHere), sep='')
  
  onesMult = (!valuesNa) %*% focalMat[,uniqueBw, drop=FALSE]
  naPoints = is.na(raster::getValues(x[[1]], row, 1L))
  onesMult[naPoints & (onesMult < 0.5)] = NA
  
  colnames(onesMult) = paste(colnames(onesMult), 'ones', sep='')
  
  toWrite = cbind(toWrite, onesMult) 
  if(mmap) {
    
    outMap = mmap::mmap(
      gsub("[.]grd$", ".gri", raster::filename(out)),
      mode=structure(numeric(0), bytes = 8L, signed = 1L, class = c("Ctype", 
          "double"))
    )	
    
    # band ordering BSQ
    indexSeq = mapply(seq, 
      from=seq(
        from= raster::cellFromRowCol(out, row, 1), 
        len=raster::nlayers(out), by=raster::ncell(out)),
      MoreArgs = list(
        by = 1, len = raster::ncol(out)
      )
    )
    
    outMap[indexSeq] = as.vector(toWrite)
    
    res = mmap::munmap(outMap)  
    
  } else {
    raster::writeValues(out, toWrite, row)
    res = 0
  }
  
  res
}


focalMult = function(
  x, 
  w, 
  bw = NULL,
  filename = paste(tempfile(), '.grd', sep=''), 
  edgeCorrect=FALSE, 
  cl = NULL
) {
  
  # check if bwXXXX is in names of x
  bwExpr = "^bw([[:digit:]]|[.])+_?"
  if( all(grepl(bwExpr, names(x))) ) {
    bw = regmatches(names(x), regexpr(bwExpr, names(x)))
    bw = gsub("[[:punct:]]$", "", bw)
    names(x) = gsub(bwExpr, "", names(x))
  }
  
  # if outer, smooths every layer of x with every focal in w
  
  if(edgeCorrect) {
    firstFile = tempfile("forEdge",tempdir(), '.grd')
  } else {
    firstFile = filename
  }
  
  if(length(bw) ) {
    if(is.numeric(bw)) bw = paste("bw", bw, sep='')
    if(!all(bw %in% names(w)))
      warning("cant find bandwidth in w")
    uniqueBw = unique(bw)
    w = w[uniqueBw]
    outNames = 
      paste(
        c(bw, uniqueBw), 
        c(names(x), rep('ones', length(uniqueBw))),
        sep=''
      )
  } else {
    # every layer gets each bandwidth
    outNames = apply(
      expand.grid(names(w), c(names(x), 'ones')), 
      1, FUN=paste, collapse='')
  }
  
  outBrick = brick(raster(x), nl = length(outNames))  
  names(outBrick) = outNames
  
  # get rid of outer parts of the focal array
  # if the values are all very small
  wArray = do.call(abind::abind, c(w, list(along=3)))
  
  Dseq = seq(1, pmax(1,floor(dim(wArray)[1]/2)-2))
  intSeq = NULL
  for(D in Dseq) {
    innerSeq = seq(D, dim(wArray)[1]+1-D)
    intSeq = cbind(intSeq,
        apply(wArray[innerSeq,innerSeq,,drop=FALSE],3,sum)
    )
  }
  
  toCrop = Dseq[max(which(apply(intSeq, 2, min) > 0.99))]
  innerSeq = seq(toCrop, dim(wArray)[1]+1-toCrop)
  wArray = wArray[
      innerSeq, innerSeq,,drop=FALSE
      ] 
  wArray = wArray / rep(intSeq[,toCrop], each=prod(dim(wArray)[1:2]))
  
  
  focalMat = matrix(
    wArray,
    ncol = dim(wArray)[3], 
    dimnames = list(NULL, dimnames(wArray)[[3]]))
  
  outBrick = writeStart(
    outBrick, 
    firstFile, 
    format = 'raster',
    datatype = 'FLT8S', 
    bandorder = 'BSQ',
    overwrite=file.exists(firstFile))
  
  toWrite = matrix(as.double(0), prod(dim(x)[1:2]), length(outNames))
  outBrick = writeValues(outBrick, toWrite, 1)
  rm(toWrite)
  
  forMoreArgs = list(
    x=x,
    focalMat = focalMat,
    focalDim = dim(wArray)[1],
    out = outBrick
  )
  
  if(length(bw) ) {
    
    forMoreArgs = c(
      forMoreArgs,
      list(
        bw = bw,
        uniqueBw = uniqueBw
      )
    )
    focalFunUse = focalMultOneRow
  } else {
    focalFunUse = focalMultOneRowOuter
  }
  
  if(!is.null(cl)) {
    outBrick = writeStop(outBrick)
    
    junk = parallel::clusterMap(
      cl,
      focalFunUse,
      row = 1:nrow(outBrick), 
      MoreArgs = c(forMoreArgs, list(mmap=TRUE))
    )
  } else {
    junk = mapply(
      focalFunUse,
      row = 1:nrow(outBrick), 
      MoreArgs = forMoreArgs
    )
    outBrick = writeStop(outBrick)
  }  
  
  
  
  if(edgeCorrect) {
    Sresult = grep("^bw([[:digit:]]|[.])+ones$", names(outBrick), invert=TRUE, value=TRUE)
    Sbw = regmatches(Sresult, regexpr("^bw([[:digit:]]|[.])+", Sresult))
    Sones = paste(Sbw, 'ones', sep='')
    
    resInd = match(Sresult, names(outBrick))
    edgeInd = match(Sones, names(outBrick))
    
    if(any(is.na(resInd)) | any(is.na(edgeInd))) {
      warning('cant find bw for edge correction')
    }
    
    res = raster::calc(
      outBrick,
      function(xx) {
        res = xx[resInd]/xx[edgeInd]
        res
      },
      filename= filename,
      overwrite = file.exists(filename)
    )
  } else {
    res = outBrick
  }
  
  setMinMax(res)
}

