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
    firstFile = paste(tempfile(), '.grd', sep='')
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
  
  focalMat = matrix(
    do.call(abind::abind, c(w, list(along=3))),
    ncol = length(w), 
    dimnames = list(NULL, names(w)))
  
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
    focalDim = dim(w[[1]]),
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

