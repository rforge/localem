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
    cat("diagonal blocks of smoothing matrix for partitions\n")
  }
  
  spatial.tools::sfQuickInit(ncores, methods = FALSE)
  
  theMat = smoothingMatrixDiag(
      rasterCoarse=rasterObjects$rasterCoarse,
      rasterFine=rasterObjects$rasterFine,
      focalList=rasterObjects$focal,
      offsetRaster=rasterObjects$offset,
      filename = filename,
      verbose=verbose)
  
  
  smoothingRaster = theMat$smoothingArray
  smoothingRasterFile = filename(smoothingRaster)
  layerSeq = 1:nlayers(smoothingRaster)
  Sbw = theMat$bw
  
  if(verbose) {
    cat(date(), "\n")
  }
  myBar = raster::pbCreate(length(theMat$uniqueDist),
      c('','text')[1+verbose], 
      label='off-diagonals of smoothing matrix')
  offDiag = foreach::foreach(
          x = theMat$uniqueDist, .packages='localEM', .export = 'smoothingMatrixOneDist'
      ) %dopar% {
        
#        parallel::mcmapply(
        thisBlock = smoothingMatrixOneDist(x,
#    	x=theMat$uniqueDist,
#    	MoreArgs=list(
            allCells=theMat$cells,
            focal=rasterObjects$focal,
            coarse=rasterObjects$rasterCoarse,
            fine=rasterObjects$rasterFine,
            offsetRaster=rasterObjects$offset
#    	),
#    	mc.cores=ncores, SIMPLIFY=FALSE,
#      mc.preschedule=FALSE
        )
        
        for(Dcell1 in 1:length(thisBlock)) {
          for(Dcell2 in 1:length(thisBlock[[Dcell1]])) {
            partHere = thisBlock[[Dcell1]][[Dcell2]]
            s1 = dimnames(thisBlock[[Dcell1]][[Dcell2]])[[1]]
            s2 = dimnames(thisBlock[[Dcell1]][[Dcell2]])[[2]]
            
            partHere = aperm(thisBlock[[Dcell1]][[Dcell2]][,,,'transpose', drop=FALSE], c(2,1,3,4))
            matchPartHere = lapply(dimnames(partHere)[1:2], match, table=theMat$partitions)
            
            haveWritten = FALSE
            writeCounter1 = 0
            while(!haveWritten & (writeCounter1 < 20)) {
              haveWritten = tryCatch(spatial.tools::binary_image_write(
                      gsub("d$", "i", smoothingRasterFile), 
                      image_dims = dim(smoothingRaster),
                      data=partHere,
                      data_position = list(
                          matchPartHere[[1]], 
                          matchPartHere[[2]], 
                          layerSeq)), 
                  error = function(err) {FALSE} )
              writeCounter1 = writeCounter1 + 1
              
            }
            if(writeCounter1 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
            
            partHere = thisBlock[[Dcell1]][[Dcell2]][,,,'straightup', drop=FALSE]
            matchPartHere = lapply(dimnames(partHere)[1:2], match, table=theMat$partitions)
            
            haveWritten = FALSE
            writeCounter2 = 0
            while(!haveWritten & (writeCounter1 < 20)) {
              haveWritten = tryCatch(spatial.tools::binary_image_write(
                      gsub("d$", "i", smoothingRasterFile), 
                      image_dims = dim(smoothingRaster),
                      data=partHere, 
                      data_position = list(
                          matchPartHere[[1]], 
                          matchPartHere[[2]], 
                          layerSeq)), 
                  error = function(err) {FALSE} )
              writeCounter2 = writeCounter2 + 1
            }
            if(writeCounter2 >= 20) warning(paste("dist", x, "cells", Dcell1, Dcell2))
            
          }
        }
        raster::pbStep(myBar)
        c(writeCounter1, writeCounter1)
      } # end foreach
      raster::pbClose(myBar)
  spatial.tools::sfQuickStop()
  
  
  
  if(FALSE) { # array for debugging
    if(verbose) {
      cat(date(), "\n")
      cat("assembling smoothing matrix\n")
    }
    
    smoothingRaster = theMat$smoothingArray
    smoothingRasterFile = filename(smoothingRaster)
    layerSeq = 1:nlayers(smoothingRaster)
    Spartitions = theMat$partitions
    Sbw = theMat$bw
    
    smoothingArray = aperm(as.array(smoothingRaster), c(2,1,3))
    dimnames(smoothingArray) = list(theMat$partitions, theMat$partitions, names(smoothingRaster))
    smoothingArray[is.na(smoothingArray)] = 0
    print(range(aperm(as.array(smoothingRaster), c(2,1,3)) - smoothingArray))
    
    # old code
    for(Ddist in seq(1, len=length(offDiag), by=1)) {
      for(Dcell1 in 1:length(offDiag[[Ddist]])) {
        for(Dcell2 in 1:length(offDiag[[Ddist]][[Dcell1]])) {
          partHere = offDiag[[Ddist]][[Dcell1]][[Dcell2]]
          s1 = dimnames(offDiag[[Ddist]][[Dcell1]][[Dcell2]])[[1]]
          s2 = dimnames(offDiag[[Ddist]][[Dcell1]][[Dcell2]])[[2]]
          
          if(FALSE) { # array for debugging
            
            if(length(s1) == 1 || length(s2) == 1) {
              smoothingArray[s2,s1,] = 
                  offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'transpose']
              
              
            } else { # neither length 1
              smoothingArray[s2,s1,] = 
                  aperm(offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'transpose'], c(2,1,3))
              
            } # done the transpose part
            smoothingArray[s1,s2,] = 
                offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'straightup']
            
          } # array for debugging
          
          partHere = aperm(offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'transpose', drop=FALSE], c(2,1,3,4))
          matchPartHere = lapply(dimnames(partHere)[1:2], match, table=Spartitions)
          
          spatial.tools::binary_image_write(
              gsub("d$", "i", smoothingRasterFile), 
              image_dims = dim(smoothingRaster),
              data=partHere,
              data_position = list(
                  matchPartHere[[1]], 
                  matchPartHere[[2]], 
                  layerSeq)
          )
          
          
          partHere = offDiag[[Ddist]][[Dcell1]][[Dcell2]][,,,'straightup', drop=FALSE]
          matchPartHere = lapply(dimnames(partHere)[1:2], match, table=Spartitions)
          
          spatial.tools::binary_image_write(
              gsub("d$", "i", smoothingRasterFile), 
              image_dims = dim(smoothingRaster),
              data=partHere, 
              data_position = list(
                  matchPartHere[[1]], 
                  matchPartHere[[2]], 
                  layerSeq)
          )
          
#        print(Dcell2)
#        print(range(aperm(as.array(smoothingRaster), c(2,1,3)) - smoothingArray))
        }
#      print(Dcell1)
#      print(range(aperm(as.array(smoothingRaster), c(2,1,3)) - smoothingArray))
      }
#    print(Ddist)
  }
#    print(range(aperm(as.array(smoothingRaster), c(2,1,3)) - smoothingArray))   
      
#fill in the zeros in the smoothing matrix where distance is beyond largest bandwidth
#  smoothingArray[is.na(smoothingArray)] = 0
      
#  range(aperm(as.array(smoothingRaster), c(2,1,3)) - smoothingArray)
      
      
#  bw = gsub("^bw", "", dimnames(theMat$smoothingArray)[[3]])
#  smoothingArrayInf = apply(!is.finite(theMat$smoothingArray), 3, any)
      smoothingArrayInf = FALSE
      if(any(smoothingArrayInf)) {
        cat("excluding bandwidths (b/c infinite values in smoothing matrix): ", bw[smoothingArrayInf],  "\n")
        
        theMat$smoothingArray = theMat$smoothingArray[,,!smoothingArrayInf,drop=FALSE]
      }
    } # end debugging
    
    
    result = c(theMat, rasterObjects[setdiff(names(rasterObjects), names(theMat))])
    
    if(verbose) {
      cat(date(), "\n")
      cat("done\n")
    }
    
    return(result)
    
  }
