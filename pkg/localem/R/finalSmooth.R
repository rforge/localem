#' @title Computes final smoothed LEM estimate
#'
#' @description Takes a partition-level LEM estiamte and smoothis it to a continuous risk estimate 
#'
#' @param x result from lemXV
#' @param counts observed counts, defaults to all counts for which cross validation was done
#' @param bw smoothing bandwidth, defaults to the optimal bandwidths for each set of counts
#' @param ncores Number of cores/threads for parallel processing
#' @param filename name of final raster file
#' @param verbose if TRUE, print progress information
#'
#' @details The optimal bandwidth for each layer is 
#' @return A raster brick
#'
#' @export
lemFinal = function(
  x, 
  counts = colnames(x$xv)[-1],
  bw = x$xv[apply(x$xv[,counts],2,which.min),'bw'], 
  ncores=1, 
  filename = paste(tempfile(), '.grd', sep=''), 
  verbose=FALSE) {
  
  finalBw = rep_len(bw, length(counts))
  
  Scounts = counts
  Slayers = paste("bw", finalBw, "_", Scounts, sep='')
  Slayers = gsub("^bwbw", "bw", Slayers)
  toSmooth = x$riskAll
  levels(toSmooth)[[1]] = levels(toSmooth)[[1]][, c("ID", Slayers)]
  toSmooth = deratify(toSmooth)
  
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
  
  
  result = focalMult(
    x=toSmooth, 
    w=x$smoothingMatrix$focal$focal, 
    edgeCorrect = TRUE,
    filename = filename,
    cl = theCluster
  )
  
  if(endCluster)  
    parallel::stopCluster(theCluster)
    
  result
  
}
