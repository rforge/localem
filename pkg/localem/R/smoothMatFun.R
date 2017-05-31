#' @title Generates the smoothing matrix for the partitions and bandwidths
#'
#' @description The \code{smoothingMatrix} function computes the entries of the smoothing matrix for the local-EM algorithm for the input bandwidths based on the partitions created by rasterizing the coarse and fine polygons, and their smoothed offsets created by applying the kernel smoothing function with input bandwidths. 
#'
#' @param rasterObjects Raster objects of partitions and smoothed offsets
#' @param ncores Number of cores/threads for parallel processing
#' @param verbose Verbose output
#'
#' @details After using the \code{smoothingMatrix} function, the smoothing matrix is an array containing the integrated kernel smoothing entries of the partitions divided by the integrated kernel smoothing entries of the study region for each specified bandwidth. 
#'  
#' @return The \code{smoothingMatrix} function returns a list containing an array of smoothing matrix entries, and the input rasters of partitions and smoothed offsets. 
#'  
#' @examples 
#' \dontrun{ 
#' data(kentuckyCounty)
#' data(kentuckyTract)
#' 
#' ncores = 1 + (.Platform$OS.type == 'unix')
#' 
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, 
#'                            polyFine = kentuckyTract, 
#'                            cellsCoarse = 6, 
#'                            cellsFine = 100, 
#'                            bw = c(10, 12, 15, 17, 20, 25) * 1000, 
#'                            ncores = ncores, 
#'                            idFile = 'id.grd', 
#'                            offsetFile = 'offset.grd', 
#'                            verbose = TRUE)
#'
#'
#' lemSmoothMat = smoothingMatrix(rasterObjects = lemRaster, 
#'                                ncores = ncores, 
#'                                verbose = TRUE)
#'}
#'
#' @export
smoothingMatrix = function(
  rasterObjects,
  ncores = 1,
  verbose = FALSE
){
	
  if(verbose) {
    cat(date(), "\n")
    cat("diagonal blocks of smoothing matrix for partitions\n")
  }
  theMat = smoothingMatrixDiag(
    	rasterCoarse=rasterObjects$rasterCoarse,
    	rasterFine=rasterObjects$rasterFine,
    	focalList=rasterObjects$focal,
    	offsetRaster=rasterObjects$offset,
    	ncores=ncores)
	
  if(verbose) {
    cat(date(), "\n")
    cat("off-diagonal blocks of smoothing matrix\n")
  }
  offDiag = parallel::mcmapply(
    	smoothingMatrixOneDist,
    	x=theMat$uniqueDist,
    	MoreArgs=list(
      		allCells=theMat$cells,
      		focal=rasterObjects$focal,
      		coarse=rasterObjects$rasterCoarse,
      		fine=rasterObjects$rasterFine,
      		offsetRaster=rasterObjects$offset
    	),
    	mc.cores=ncores, SIMPLIFY=FALSE,
      mc.preschedule=FALSE
  )
	
  if(verbose) {
    cat(date(), "\n")
    cat("assembling smoothing matrix\n")
  }
	
  for(Ddist in seq(1, len=length(offDiag), by=1)) {
    for(Dcell1 in 1:length(offDiag[[Ddist]])) {
      for(Dcell2 in 1:length(offDiag[[Ddist]][[Dcell1]])) {
				
        s1 = dimnames(offDiag[[Ddist]][[Dcell1]][[Dcell2]])[[1]]
        s2 = dimnames(offDiag[[Ddist]][[Dcell1]][[Dcell2]])[[2]]
				
        if(length(s1) == 1 || length(s2) == 1) {
          theMat$smoothingArray[s2,s1,] = 
							offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'transpose']
        } else {
          theMat$smoothingArray[s2,s1,] = 
							aperm(offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'transpose'], c(2,1,3))
        }
        theMat$smoothingArray[s1,s2,] = 
						offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'straightup']
				
      }
    }
  }
	
  #fill in the zeros in the smoothing matrix where distance is beyond largest bandwidth
  theMat$smoothingArray[is.na(theMat$smoothingArray)] = 0
	
  bw = gsub("^bw", "", dimnames(theMat$smoothingArray)[[3]])
  smoothingArrayInf = apply(!is.finite(theMat$smoothingArray), 3, any)
	
  if(any(smoothingArrayInf)) {
    cat("excluding bandwidths (b/c infinite values in smoothing matrix): ", bw[smoothingArrayInf],  "\n")
		
    theMat$smoothingArray = theMat$smoothingArray[,,!smoothingArrayInf,drop=FALSE]
  }
	
  result = c(theMat, rasterObjects[setdiff(names(rasterObjects), names(theMat))])

  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
  
  return(result)
  
}
