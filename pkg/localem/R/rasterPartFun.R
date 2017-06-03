#' @title Generates the partitions by overlaying rasters of the coarse and fine polygons
#'
#' @description The \code{rasterPartition} function first rasterizes the coarse and fine spatial polygons based on their respective input resolutions, and then, overlays these rasters to generate the raster of partitions of the local-EM algorithm. It also applies the kernel smoothing function with input bandwidths to the expected counts of the fine polygons to obtain the smoothed offsets (i.e., smoothed expected counts / cell area) of the partitions. 
#' 
#' @param polyCoarse Spatial polygons of case data
#' @param polyFine Spatial polygons of population data
#' @param bw Vector of bandwidths
#' @param focalSize Distance to truncate Gaussian kernel bandwidth (default is 3 times of the bandwidth)
#' @param cellsCoarse Horizontal/vertical resolution of raster applied to coarse polygons
#' @param cellsFine Horizontal/vertical resolution of raster applied to fine polygons
#' @param ncores Number of cores/threads for parallel processing
#' @param idFile Filename (must have .grd extension) of the raster of partitions
#' @param offsetFile Filename (must have .grd extension) of the rasters of smoothed offsets
#' @param verbose Verbose output
#' 
#' 
#' @details After using the \code{rasterPartition} function, the fine raster is a raster stack containing the IDs for the partitions created by overlaying the coarse and fine rasters. The offset raster is a raster stack containing the offsets of the partitions smoothed with the specified bandwidths. These values represent the denominator of the kernel smoothing matrix. 
#'
#' @return The \code{rasterPartition} function returns a list containing the raster of the coarse polygons, raster stacks of the partitions and offsets, focal weight matrix of the Gaussian kernel, and the input coarse polygons. 
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
#'}
#'
#' @export
rasterPartition = function(
    polyCoarse, 
    polyFine, 
    cellsCoarse, 
    cellsFine, 
    bw,
    focalSize = NULL,
    xv = NULL, 
    ncores = 1, 
    path = tempdir(),
    idFile = paste(tempfile('lemId', path), '.grd', sep=''), 
    offsetFile = paste(tempfile('lemOffset', path), '.grd', sep=''),
    verbose = FALSE
){
  
  if(verbose) {
    cat(date(), "\n")
    cat("obtaining rasters\n")
  }
  
  if(!(length(grep("\\.grd$", offsetFile)))){
    warning("offsetFile should have .grd extension")
  }
  if(!(length(grep("\\.grd$", idFile)))){
    warning("idFile should have .grd extension")
  }
  if(is.numeric(cellsCoarse)) {
    rasterCoarse=geostatsp::squareRaster(polyCoarse, cellsCoarse)
    values(rasterCoarse) = seq(1, ncell(rasterCoarse))
    names(rasterCoarse) = "cellCoarse"
  } else {
    rasterCoarse = cellsCoarse
  }
  
  if(is.numeric(cellsFine)) {
    rasterFine = disaggregate(rasterCoarse,
        ceiling(cellsFine/ncol(rasterCoarse)))
    names(rasterFine) = 'cellCoarse'
  } else {
    rasterFine = cellsFine
  }
  
  # coarse poly ID's for fine raster
  
  polyCoarseIdCol = grep("^id$", names(polyCoarse), value=TRUE)
  if(!length(polyCoarseIdCol)) {
    polyCoarseIdCol = names(polyCoarse)[1]
  }
  
  idCoarse = 1:length(polyCoarse)
  names(idCoarse) = polyCoarse@data[[polyCoarseIdCol]]
  polyCoarse$idCoarse = idCoarse
  
  rasterIdCoarse = rasterize(
      polyCoarse,
      rasterFine,
      field='idCoarse')
  names(rasterIdCoarse) = 'idCoarse'
  
  polyFine@data[is.na(polyFine@data[,'expected']),'expected'] = 0
  
  if(length(Nxv)==1) {
    xvMat = getXvMat(polyCoarse$idCoarse, Nxv)
  } else {
    if(is.null(Nxv)) {
      xvMat = matrix()
    } else {
      xvMat = Nxv
    }
  }
  
  
  rasterOffset = geostatsp::spdfToBrick(
      x=polyFine,
      template=rasterFine,
      pattern='^expected$')
  names(rasterOffset) = 'offset'
  
  # fine ID (iffset
  # scale the offsets to cases per cell
  # times 10 (roughly)
  # then create fine ID's with homogeneous offsets
  maxOffset = 10^5 / maxValue(rasterOffset)
  ratifyOffset = ratify(round(rasterOffset*maxOffset), count=FALSE)
  stuff= levels(ratifyOffset)[[1]]
  stuff$idFine = seq(1, nrow(stuff))
  levels(ratifyOffset) = stuff
  
  rasterIdFine = deratify(ratifyOffset, 'idFine')
  
  for(Dxv in 1:ncol(xvMat)) {
    xvHere = which(xvMat[,Dxv])
    maskHere = calc(rasterIdCoarse, 
        function(x) ! x %in% xvHere )
    rasterOffset = addLayer(
        rasterOffset, 
        rasterOffset[[1]] * maskHere
    )
  }
  names(rasterOffset) = c("offset", 
      paste('xvOffset', colnames(xvMat), sep=''))
  
  offsetTempFile = paste(tempfile(), '.grd', sep='')
  
  rasterOffset = raster::mask(rasterOffset, rasterIdCoarse[['idCoarse']],
      filename=offsetTempFile, overwrite = file.exists(offsetTempFile))
  
  
  rasterFineId = brick(rasterIdCoarse, rasterIdFine, rasterFine)	
  
  rasterFineId = writeRaster(rasterFineId, idFile,
      overwrite=file.exists(idFile))
  
  if(verbose) {
    cat(date(), "\n")
    cat("smoothing offsets\n")
  }
  
  theFocal = focalFromBw(
      bw = bw, 
      fine=rasterFine, 
      focalSize=focalSize, 
      ncores=ncores)
  
  
  forSmooth = expand.grid(
      bw=names(theFocal$focal),
      layer=names(rasterOffset))
  
  focalArray = array(
      unlist(theFocal$focal),
      c(dim(theFocal$focal[[1]]), length(theFocal$focal)),
      dimnames = list(NULL,NULL, names(theFocal$focal))
  )
  focalArray = focalArray[,,forSmooth[,'bw']]
  dimnames(focalArray)[[3]] = paste(
      forSmooth[,'bw'], 
      gsub("offset", "", forSmooth[,'layer'], ignore.case=TRUE),
      sep=''
  )
  theFocal$array = focalArray   
  
  Scvsets = match(forSmooth[,'layer'],names(rasterOffset))
  
  focalFunction = function(x, focalArray, Scvsets)  {
    
    apply(focalArray*x[,,Scvsets], 
        3, sum, na.rm=TRUE)
  }
  
  rasterOffset = setMinMax(rasterOffset)
  
  spatial.tools::sfQuickInit(ncores, methods = FALSE)
  suppressWarnings(
      smoothedOffset <- spatial.tools::rasterEngine(
          rasterOffset, focalFunction, 
          args = list(Scvsets=Scvsets, focalArray=focalArray),
          window_dims = dim(theFocal$focal[[1]]),
          outbands=nrow(forSmooth),
          outfiles = 1,
          processing_unit = 'single',
          setMinMax = TRUE,
      ))
  spatial.tools::sfQuickStop()
  names(smoothedOffset) = dimnames(focalArray)[[3]]
  
  
  offsetStack = writeRaster(addLayer(rasterOffset, smoothedOffset),
      filename = offsetFile, overwrite = file.exists(offsetFile))
  
  
# create list of partitions
  partitions = as.data.frame(na.omit(unique(rasterFineId)))
  partitions$partition = paste('c', partitions$cellCoarse, 'p', partitions$idCoarse,
      '.', partitions$idFine, sep='')
  partitions = cbind(ID = 1:nrow(partitions), partitions)
# raster with partition ID's
  
  partitionRaster = raster::calc(rasterFineId, 
      function(x) which(
            x[1]==partitions[,'idCoarse'] &
                x[2] == partitions[,'idFine'] &
                x[3] == partitions[,'cellCoarse']
        )[1])
  levels(partitionRaster)[[1]] = partitions
  
# offsetMat is the value of the offset at all points in the partition
  
  partitionOffsets = as.data.frame(zonal(
          offsetStack[[grep("^bw", names(offsetStack), invert=TRUE)]],
          partitionRaster,
          'mean', na.rm=TRUE
      ))
  partitionOffsets$partition = partitions[match(
          partitionOffsets$zone, partitions$ID
      ), 'partition']    
  
  
  offsetMat = apply(partitionOffsets[,grep("[oO]ffset", colnames(partitionOffsets))], 2, 
      function(x) {
        res = Matrix::Diagonal(nrow(partitionOffsets), x*prod(res(offsetStack)))
        dimnames(res) = list(partitionOffsets$partition, partitionOffsets$partition) 
        res
      })
  
  regions = gsub("^c[[:digit:]]+p|\\.[[:digit:]]+$", "", partitionOffsets$partition)
  regions = as.integer(regions)
  regionMat = outer(regions, regions, '==')
  regionMat = Matrix::Matrix(regionMat)
  
  dimnames(regionMat)=
      list(partitionOffsets$partition, partitionOffsets$partition)
  
  partitionFreq = raster::freq(partitionRaster)
  rownames(partitionFreq) = levels(partitionRaster)[[1]][partitionFreq[,'value'],'partition']
  partitionAreas = partitionFreq[,'count'][rownames(regionMat)] * prod(res(partitionRaster))  
  
  expandCountMat = regionMat
  expandCountMat = expandCountMat[,!duplicated(regions),drop=FALSE]
  colnames(expandCountMat) = gsub("^c[[:digit:]]+p|\\.[[:digit:]]+$", 
      "", colnames(expandCountMat) )
  expandCountMat = expandCountMat[,
      as.character(sort(as.integer(colnames(expandCountMat)))),
      drop=FALSE
  ]
  
    
  partitionRaster = writeRaster(partitionRaster, file=idFile, overwrite=file.exists(idFile))
  
  
  result = list(
      rasterCoarse=rasterCoarse,
      rasterFine=partitionRaster,
      focal=theFocal,
      offset=offsetStack,
      partitionAreas=partitionAreas,
      xv = xvMat
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
  
  return(result)
  
}
