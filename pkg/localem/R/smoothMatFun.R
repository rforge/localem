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
  filename = 'smoothingMatrix.gri',
  verbose = FALSE
){
  
  if(verbose) {
    cat(date(), "\n")
    cat("diagonal blocks of smoothing matrix\n")
    cat("if there are errors stop the cluster with spatial.tools::sfQuickStop()\n")
  }
  
  if(ncores > 1) spatial.tools::sfQuickInit(ncores, methods = FALSE)
  
  
  theMat = smoothingMatrixDiag(
    rasterCoarse=rasterObjects$rasterCoarse,
    rasterFine=rasterObjects$rasterFine,
    focalList=rasterObjects$focal,
    offsetRaster=rasterObjects$offset,
    filename = filename,
    verbose=verbose)
  
  
  smoothingRaster = theMat$smoothingArray
  smoothingRasterFile = gsub("grd$", "gri", filename(smoothingRaster))
  layerSeq = 1:nlayers(smoothingRaster)
  Sbw = theMat$bw
  Spartitions = theMat$partitions
  
  if(verbose) {
    cat('off-diagonals of smoothing matrix', "\n")
    cat(date(), "\n")
  }
  
  # theType = mmap::real64();dput(theType, '')
  theType = structure(numeric(0), bytes = 8L, signed = 1L, class = c("Ctype", 
      "double"))
  
  offDiag = foreach::foreach(
      x = as.vector(theMat$uniqueDist), .packages='localEM', .export = 'smoothingMatrixOneDist'
    ) %dopar% {
      
      thisBlock = try(smoothingMatrixOneDist(x,
          allCells=theMat$cells,
          focal=rasterObjects$focal,
          coarse=rasterObjects$rasterCoarse,
          fine=rasterObjects$rasterFine,
          offsetRaster=rasterObjects$offset
        ), silent=TRUE)
      
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
            
            haveWritten = FALSE
            writeCounter1 = 0
            while(!haveWritten & (writeCounter1 < 20)) {
              haveWritten = tryCatch(spatial.tools::binary_image_write(
                  smoothingRasterFile,
                  mode = theType, 
                  image_dims = dim(smoothingRaster),
                  data=as.double(partHere[theOrder[[1]], theOrder[[2]],,] ),
                  data_position = list(
                    as.vector(matchPartHere[[1]][theOrder[[1]] ]), 
                    as.vector(matchPartHere[[2]][theOrder[[2]] ]), 
                    layerSeq)), 
                error = function(err) {warning(err);-1} )
              haveWritten = (haveWritten != -1)
              writeCounter1 = writeCounter1 + 1
            }
            if(writeCounter1 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
            
            partHere = try(thisBlock[[Dcell1]][[Dcell2]][,,,'straightup', drop=FALSE])
            matchPartHere = lapply(dimnames(partHere)[1:2], match, 
              table=Spartitions)
            theOrder = lapply(matchPartHere, order)
            
            haveWritten = FALSE
            writeCounter2 = 0
            while(!haveWritten & (writeCounter2 < 20)) {
              haveWritten = tryCatch(spatial.tools::binary_image_write(
                  smoothingRasterFile, 
                  mode = theType, 
                  image_dims = dim(smoothingRaster),
                  data=as.double(partHere[theOrder[[1]], theOrder[[2]],,] ),
                  data_position = list(
                    as.vector(matchPartHere[[1]][theOrder[[1]] ]), 
                    as.vector(matchPartHere[[2]][theOrder[[2]] ]), 
                    layerSeq)), 
                error = function(err) {warning(err);-1} )
              haveWritten = (haveWritten != -1)
              writeCounter2 = writeCounter2 + 1
            }
            if(writeCounter2 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
          } # end Dcell2  
        } # end Dcell1
        res = c(writeCounter1, writeCounter1)
      }  else { # end try error
        res =thisBlock
      }
    } # end foreach
  
  if(any(unlist(offDiag) >= 20)) warning("problem writing smoothing matrix to disk")
  
  if(ncores > 1)  spatial.tools::sfQuickStop()
  
  
  if(verbose) {
    cat(date(), "\n")
    cat("replacing NA's with zeros\n")
  }
  
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
  
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



