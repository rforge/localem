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
    ncores=1, filename = tempfile(), verbose=FALSE) {
  

  finalBw = rep_len(bw, length(counts))
  
  Scounts = counts
  Slayers = paste("bw", finalBw, "_", Scounts, sep='')
  
  xFocal = x$smoothingMatrix$focal$array[,,paste('bw', finalBw, sep=''), drop=FALSE]
 
  # smooth the risk and integrate kernel over non-NA area
  focalFunction = function(x, fa)  {
        apply(fa*x, 3, sum, na.rm=TRUE) / 
        apply(fa*(!is.na(x)), 3, sum)
  }
  
  toSmooth = x$smoothingMatrix$rasterFine
  levels(toSmooth)[[1]] = levels(x$smoothingMatrix$rasterFine)[[1]][,
     c("ID", Slayers)]
  toSmooth = deratify(toSmooth)
  
  if(ncores>1) spatial.tools::sfQuickInit(
        ncores, methods = TRUE,        
        .packages = c('Matrix','raster'))
  
  
  suppressWarnings(
      smoothedRisk <- spatial.tools::rasterEngine(
          x=toSmooth, fun=focalFunction, 
          args = list(fa=xFocal),
          window_dims = dim(xFocal),
          outbands=dim(xFocal)[3],
          outfiles = 1,
          processing_unit = 'single',
          chunk_format = 'array',
          filename = gsub("[.]gr(d|i)$", "", filename), overwrite=TRUE,
          verbose=(verbose>2)
      ))
  if(ncores>1) spatial.tools::sfQuickStop()
  names(smoothedRisk) = Slayers
  
  smoothedRisk
  
}
