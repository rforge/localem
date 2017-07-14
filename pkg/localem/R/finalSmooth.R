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
  
  xFocal = x$smoothingMatrix$focal$focal[paste("bw", finalBw, sep='')]
  xFocal = do.call(abind::abind, c(xFocal, list(along=3)))
  
  # smooth the risk and integrate kernel over non-NA area
  focalFunction = function(x, fa)  {
    apply(fa*x, 3, sum, na.rm=TRUE) / 
      apply(fa*(!is.na(x)), 3, sum)
  }
  
  toSmooth = x$riskAll
  levels(toSmooth)[[1]] = levels(toSmooth)[[1]][, c("ID", Slayers)]
  toSmooth = deratify(toSmooth)

  Soutfile = file.path(tempdir(), paste('finalSmooth', Slayers, '.grd', sep='')) 
  names(Soutfile) = Slayers
  
  endCluster = FALSE
  theCluster = NULL
  if(length(grep("cluster", class(ncores))) ) {
    if(verbose) cat("using existing cluster\n")
    theCluster = ncores
  } else if(!is.null(ncores)) {
    if(ncores > 1) {
      if(verbose) cat("starting new cluster\n")
      theCluster = parallel::makeCluster(spec=ncores, type='PSOCK', methods=TRUE)
      parallel::setDefaultCluster(theCluster)
      endCluster = TRUE
    }
  }
  
  
  oneFinalFun = function(Dsmooth, toSmooth, xFocal, Soutfile) {
    res = raster::focal(
      toSmooth[[Dsmooth]],
      w = xFocal[,,gsub("_.*$", "", Dsmooth)],
      na.rm=TRUE, pad=TRUE,
      filename = Soutfile[Dsmooth],
      overwrite = file.exists(Soutfile[Dsmooth])
    )
    filename(res)
  }

  forMoreArgs = list(
    toSmooth = toSmooth, xFocal = xFocal, Soutfile = Soutfile
    )
  
  if(!is.null(theCluster)) {
    foreachResult = parallel::clusterMap(
      theCluster, 
      oneFinalFun,
      Dsmooth = Slayers,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  } else {
    foreachResult = mapply(
      oneFinalFun,
      Dsmooth = Slayers,
      MoreArgs = forMoreArgs,
      SIMPLIFY=FALSE)
  }

    
  
#  foreachResult = foreach(
#      Dsmooth = Slayers, .packages='raster'
#    ) %dopar% { 
#      res = raster::focal(
#        toSmooth[[Dsmooth]],
#        w = xFocal[,,gsub("_.*$", "", Dsmooth)],
#        na.rm=TRUE, pad=TRUE,
#        filename = Soutfile[Dsmooth],
#        overwrite = file.exists(Soutfile[Dsmooth])
#      )
#      filename(res)
#    }
  
  
  if(endCluster) parallel::stopCluster(theCluster)
  
  smoothedStack = raster::stack(foreachResult)
  names(smoothedStack) = Scounts
  
  result = raster::brick(
    smoothedStack, filename = filename, overwrite = file.exists(filename)
  )
  names(result) = Scounts
  
  if(FALSE) {
    suppressWarnings(
      result <- spatial.tools::rasterEngine(
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
  }
  
  result
  
}
