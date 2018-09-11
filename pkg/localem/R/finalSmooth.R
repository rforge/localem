
# Computes the risk estimation to the fine raster by smoothing the risk of the partitions at the final iteration
finalSmooth = function(
  x,
  Slayers,
  fact=1,
  ncores=1,
  filename=tempfile("final", tempdir(), ".grd"),
  focalMat = NULL
  ) {

  toSmooth = x$riskAll
  levels(toSmooth)[[1]] = levels(toSmooth)[[1]][, c("ID", Slayers)]

  toSmooth = deratify(toSmooth, 
    filename = tempfile("deratify", tempdir(), ".grd"))

  if(is.null(focalMat)) {
    theFocal = x$smoothingMatrix$focal$focal
  } else {
    theFocal = focalMat
  }

  if(fact > 1) {
    toSmooth = raster::aggregate(
      toSmooth, fact=fact,fun=mean,
      filename = tempfile("deratify", tempdir(), ".grd")
      )
    theFocal = lapply(theFocal,
      function(xx) {
        dim1 = round( (dim(xx)-1)/2)
        seqHere = mapply(
          seq,
          to = dim1,
          MoreArgs = list(from=0, by=fact)
          )
        seqHere = rbind(
          -seqHere, seqHere
          )
        seqHere = apply(seqHere, 2, function(xxx) sort(unique(xxx)))
        seqHere = seqHere + 1 + 
        matrix(dim1, nrow=nrow(seqHere), ncol=2, byrow=TRUE)
        res = xx[seqHere[,1], seqHere[,2]]  
        res/sum(res)
      })

  }

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

  theFinalEst = focalMult(
    x=toSmooth,
    w=theFocal,
    edgeCorrect = TRUE,
    filename = tempfile('finalsmoothed',
      tempdir(), '.grd'),
    cl = theCluster
    )
  if(any(grepl("deratify", filename(toSmooth)))) {
    unlink(gsub("[[:alpha:]]$", "*", filename(toSmooth)))
  }

  names(theFinalEst) = Slayers

  # done with the cluster
  if(endCluster)
    parallel::stopCluster(theCluster)

  if(fact > 1) {
    result = raster::projectRaster(
      theFinalEst,
      x$riskAll,
      method = 'bilinear',
      filename = filename,
      overwrite = file.exists(filename)
      )
  } else {
    result = raster::writeRaster(theFinalEst,
     filename = filename,
     overwrite = file.exists(filename)
     )
  }
  # remove temporary raster files
  unlink(gsub("[[:alpha:]]$", "*", filename(theFinalEst)))

  return(result)
}

# Computes the risk estimation to the fine raster by smoothing the risk of the partitions at the final iteration
finalSmoothOldQuestionmark = function(
  x,
  Slayers,
  fact=1,
  ncores=1,
  filename=tempfile("final", tempdir(), ".grd")
  ) {

  toSmooth = x$riskAll
  levels(toSmooth)[[1]] = levels(toSmooth)[[1]][, c("ID", Slayers)]

  toSmooth = deratify(toSmooth, 
    filename = tempfile("deratify", tempdir(), ".grd"))

  theFocal = x$smoothingMatrix$focal$focal

  if(fact > 1) {
    toSmooth = raster::aggregate(
      toSmooth, fact=fact,fun=mean,
      filename = tempfile("deratify", tempdir(), ".grd")
      )
    theFocal = lapply(theFocal,
      function(xx) {
        dim1 = round( (dim(xx)-1)/2)
        seqHere = mapply(
          seq,
          to = dim1,
          MoreArgs = list(from=0, by=fact)
          )
        seqHere = rbind(
          -seqHere, seqHere
          )
        seqHere = apply(seqHere, 2, function(xxx) sort(unique(xxx)))
        seqHere = seqHere + 1 + 
        matrix(dim1, nrow=nrow(seqHere), ncol=2, byrow=TRUE)
        res = xx[seqHere[,1], seqHere[,2]]  
        res/sum(res)
      })

  }

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
  
  tempFileFocal = tempfile('finalsmoothed',
   tempdir(), '.grd')
  theFinalEst = focalMult(
    x=toSmooth,
    w=theFocal,
    edgeCorrect = TRUE,
    filename = tempFileFocal,
    cl = theCluster
    )

  names(theFinalEst) = Slayers

  # done with the cluster
  if(endCluster)
    parallel::stopCluster(theCluster)

  if(fact > 1) {
    result = raster::projectRaster(
      theFinalEst,
      x$riskAll,
      method = 'bilinear',
      filename = filename,
      overwrite = file.exists(filename)
      )
  } else {
    result = raster::writeRaster(theFinalEst,
     filename = filename,
     overwrite = file.exists(filename)
     )
  }
  # remove temporary raster files
  unlink(gsub("[[:alpha:]]$", "*", filename(theFinalEst)))
  unlink(gsub("grd$", "gri", tempFileFocal))
  return(result)
}
