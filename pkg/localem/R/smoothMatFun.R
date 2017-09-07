#' @title Generates the smoothing matrix for the partitions and bandwidths
#'
#' @description The \code{smoothingMatrix} function computes the entries of the smoothing matrix for the local-EM algorithm for the input bandwidths based on the partitions created by rasterizing the spatial polygons of case and population data, and their smoothed offsets created by applying the kernel smoothing function with input bandwidths. 
#'
#' @param rasterObjects Raster stacks for partitions and smoothed offsets
#' @param ncores Number of cores/threads for parallel processing
#' @param path Folder to store raster data
#' @param filename Filename (must have .grd extension) of the entries of smoothing matrix
#' @param verbose Verbose output
#'
#' @details After using the \code{smoothingMatrix} function, the smoothing matrix is an array containing the integrated kernel smoothing entries of the partitions divided by the integrated kernel smoothing entries of the study region for each specified bandwidth. 
#'  
#' @return The \code{smoothingMatrix} function returns a list containing an array of smoothing matrix entries, and the input rasters of partitions and smoothed offsets. 
#'  
#' @examples 
#' \dontrun{ 
#' # case and population data
#' data('kentuckyCounty')
#' data('kentuckyTract')
#'
#' # parameters
#' ncores = 2
#' cellsCoarse = 8
#' cellsFine = 100
#' bw = c(10, 15, 17.5, 20) * 1000
#' path = 'example'
#' 
#' # rasters of case and population data
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, 
#'								polyFine = kentuckyTract, 
#'								cellsCoarse = cellsCoarse, 
#'								cellsFine = cellsFine, 
#'								bw = bw, 
#'								ncores = ncores, 
#'								path = path, 
#'								idFile = 'lemId.grd', 
#'								offsetFile = 'lemOffsets.grd', 
#'								verbose = TRUE)
#'
#' # smoothing matrix
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'									ncores = ncores, 
#'									path = path, 
#'									filename = 'lemSmoothMat.grd', 
#'									verbose = TRUE)
#'}
#'
#' @export
smoothingMatrix = function(
    rasterObjects,
    ncores = 1,
    path = getwd(), 
    filename, 
    verbose = FALSE
){
  
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  if(missing(filename)) {
    filename = paste(tempfile('lemSmoothMat', path), '.grd', sep = '')
  }
  if(!length(grep('/', filename))) {
    filename = file.path(path, filename)
  }
  if(!(length(grep("\\.gr[id]$", filename)))){
    warning("filename should have .grd extension")
  }
  
  if(verbose) {
    cat(date(), "\n")
    cat("diagonal blocks of smoothing matrix\n")
  }
  
  
  endCluster = FALSE
  theCluster = NULL
  if(length(grep("cluster", class(ncores))) ) {
    theCluster = ncores
    parallel::clusterEvalQ(theCluster, library('raster'))
    parallel::clusterEvalQ(theCluster, library('Matrix'))
  } else if(!is.null(ncores)) {
    if(ncores > 1) {
      theCluster = parallel::makeCluster(spec=ncores, type='PSOCK', methods=TRUE)
      parallel::setDefaultCluster(theCluster)
      parallel::clusterEvalQ(theCluster, library('raster'))
      parallel::clusterEvalQ(theCluster, library('Matrix'))
      endCluster = TRUE
    }
  }
  
  theMat = smoothingMatrixDiag(
      rasterCoarse=rasterObjects$rasterCoarse,
      rasterFine=rasterObjects$rasterFine,
      focalList=rasterObjects$focal,
      offsetRaster=rasterObjects$offset,
      filename = filename,
      cl = theCluster,
      verbose=verbose)
  
  smoothingRaster = theMat$smoothingArray
  smoothingRasterFile = gsub("grd$", "gri", filename(smoothingRaster))
  layerSeq = 1:nlayers(smoothingRaster)
  Sbw = theMat$bw
  Spartitions = theMat$partitions
  
  if(verbose) {
    cat(date(), "\n")
    cat('off-diagonals of smoothing matrix', "\n")
  }
  
  
  
  # dimensions, from spatial.tools::binary_image_write
  image_dims = dim(smoothingRaster)
  
  
  forMoreArgs = list(
      theMat = theMat, 
      rasterObjects = rasterObjects,
      image_x=image_dims[1],
      image_y=image_dims[2],
      image_z=image_dims[3],
      Spartitions = Spartitions,
      # theType = mmap::real64();dput(theType, '')
      theType = structure(numeric(0), bytes = 8L, signed = 1L, class = c("Ctype", 
              "double")),
      smoothingRasterFile = smoothingRasterFile,
      layerSeq = layerSeq,
      verbose=verbose
  )
  
  if(!is.null(theCluster)) {
    
    parallel::clusterExport(theCluster, 
        varlist = c('smoothingMatrixOneDist',
            'smoothingMatrixEntries',
            'kernMat',
            'reorderCellsTranslate',
            'getCellInfo'), 
        envir = environment()
    )
    
    offDiag = parallel::clusterMap(
        theCluster,
        oneBlockOffdiagFun,
        x = as.vector(theMat$uniqueDist),
        MoreArgs = forMoreArgs
    )
  } else {
    offDiag = mapply(
        oneBlockOffdiagFun,
        x = as.vector(theMat$uniqueDist),
        MoreArgs = forMoreArgs
    )
    
  }  
  
  if(any(unlist(offDiag) >= 20)) warning("problem writing smoothing matrix to disk")
  
  if(endCluster)  
    parallel::stopCluster(theCluster)
  
  
  if(any(unlist(lapply(offDiag, class)) == 'try-error') ) {
    warning("errors in smoothing matrix construction")
    return(c(list(offDiag = offDiag), 
            theMat, rasterObjects[setdiff(names(rasterObjects), names(theMat))]))
  }
  
  result = c(theMat, rasterObjects[setdiff(names(rasterObjects), names(theMat))])
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
  
  return(result)
}

