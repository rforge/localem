
# Computes the risk estimation to the fine raster by smoothing the risk of the partitions at the final iteration
finalSmooth = function(
    x, 
    Slayers, 
    ncores, 
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
  theFinalEst = focalMult(
      x=toSmooth, 
      w=xOrig$smoothingMatrix$focal$focal, 
      edgeCorrect = TRUE,
      filename = paste(tempfile(), '.grd', sep = ''), 
      cl = theCluster
  )

  names(theFinalEst) = Slayers
  
  # done with the cluster
  if(endCluster)  
    parallel::stopCluster(theCluster)
  
  result = raster::writeRaster(theFinalEst, 
			filename = filename, 
			overwrite = file.exists(filename)
			)
  
  return(result)
}
