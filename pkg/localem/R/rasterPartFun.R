#' @title Generates the partitions by overlaying rasters of spatial polygons
#'
#' @description The \code{rasterPartition} function first rasterizes the spatial polygons of the case and population data based on their respective input resolutions, and then, overlays these rasters to generate the raster of partitions of the local-EM algorithm. It also applies the kernel smoothing function with input bandwidths to the expected counts of the fine polygons to obtain the smoothed offsets (i.e., smoothed expected counts per cell area) of the partitions. If cross-validation is specified, smoothed offsets are computed for each cross-validation set.  
#' 
#' @param polyCoarse Spatial polygons of case data
#' @param polyFine Spatial polygons of population data
#' @param cellsCoarse Minimum resolution for rasterization of case data for numerical accuracy of smoothing matrix
#' @param cellsFine Minimum resolution for rasterization of population data for numerical integration of smoothing matrix
#' @param bw Vector of bandwidths
#' @param focalSize Distance to truncate Gaussian kernel, default is 2.5 times largest bandwidth
#' @param xv (Optional) Number of cross-validation sets, or matrix where rows are coarse polygons and columns are cross-validation sets
#' @param ncores Number of cores/threads for parallel processing
#' @param path Folder to store raster data
#' @param idFile Filename (must have .grd extension) of the raster of partitions
#' @param offsetFile Filename (must have .grd extension) of the rasters of smoothed offsets
#' @param verbose Verbose output
#' 
#' 
#' @details After using the \code{rasterPartition} function, the partition raster is a raster stack containing the IDs for the partitions created by overlaying the rasterizations of the spatial polygons of case and population data. The offset raster is a raster stack containing the offsets of the partitions smoothed with the specified bandwidths. These values represent the denominator of the kernel smoothing matrix. 
#'
#' @return The \code{rasterPartition} function returns a list containing the raster of the case data, raster stacks of the partitions and offsets, focal weight matrix of the Gaussian kernel, and the input coarse polygons. 
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
  path = getwd(),
  idFile, 
  offsetFile,
  verbose = FALSE
){
  
	dir.create(path, showWarnings = FALSE, recursive = TRUE)

	
	if(missing(idFile)) {
		idFile = paste(tempfile('lemId', path), '.grd', sep = '')
	}
	if(!length(grep('/', idFile))) {
		idFile = file.path(path, idFile)
	}
	
	if(missing(offsetFile)) {
		offsetFile = paste(tempfile('lemOffset', path), '.grd', sep = '')
	}
	if(!length(grep('/', offsetFile))) {
		offsetFile = file.path(path, offsetFile)
	}


	
  if(verbose) {
    cat(date(), "\n")
    cat("obtaining rasters\n")
  }
  
  if(!(length(grep("\\.gr[id]$", offsetFile)))){
    warning("offsetFile should have .grd extension")
  }
  if(!(length(grep("\\.gr[id]$", idFile)))){
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
  
  if(length(xv)==1) {
    xvMat = getXvMat(polyCoarse$idCoarse, xv)
  } else {
    if(is.null(xv)) {
      xvMat = matrix()
	  dimnames(xvMat) = list(1, 1)
    } else {
      xvMat = xv
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
  stuff= raster::levels(ratifyOffset)[[1]]
  if(!is.null(stuff)) {
    stuff$idFine = seq(1, nrow(stuff))
    levels(ratifyOffset) = list(stuff)
  } else {
    warning('offset problems, cant find levels')
  }
  rasterIdFine = deratify(ratifyOffset, 'idFine')
  
  for(Dxv in 1:ncol(xvMat)) {
    xvHere = which(xvMat[,Dxv])
    maskHere = raster::calc(rasterIdCoarse, 
      function(x) ! x %in% xvHere )
    rasterOffset = addLayer(
      rasterOffset, 
      rasterOffset[[1]] * maskHere
    )
  }
  names(rasterOffset) = c("offset", 
    paste('xvOffset', colnames(xvMat), sep = ''))
  
  offsetTempFile = file.path(path, 'offsetTemp.grd')
  
  rasterOffset = raster::mask(rasterOffset, rasterIdCoarse[['idCoarse']],
    filename=offsetTempFile, overwrite = file.exists(offsetTempFile))
  
  
  rasterFineId = brick(rasterIdCoarse, rasterIdFine, rasterFine)	
  
  rasterFineId = writeRaster(rasterFineId, 
    filename = idFile,
    overwrite=file.exists(idFile))
  
  if(verbose) {
    cat(date(), "\n")
    cat("computing focal array\n")
  }
  
  rasterOffset = setMinMax(rasterOffset)
  
  
  if(is.null(focalSize))
    focalSize = 2.5*max(bw)
  
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
  
  if(verbose) cat("creating arrays\n")
  
  # focal used for smoothing matrix
  theFocal = focalFromBw(
    bw = bw, 
    fine=rasterFine, 
    focalSize=focalSize,
    cl=theCluster)
  
  forSmooth = expand.grid(
    bw=names(theFocal$focal),
    layer=names(rasterOffset))
  
  focalArray = array(
    unlist(theFocal$focal),
    c(dim(theFocal$focal[[1]]), length(theFocal$focal)),
    dimnames = list(NULL,NULL, names(theFocal$focal))
  )
  focalArray = focalArray[,,forSmooth[,'bw'], drop=FALSE]
  dimnames(focalArray)[[3]] = paste(
    forSmooth[,'bw'], 
    gsub("offset", "", as.character(forSmooth[,'layer']), ignore.case=TRUE),
    sep = ''
  )
  theFocal$array = focalArray   
  
  if(verbose) cat(".\n")
  
  
  if(verbose) {
    cat(date(), "\n")
    cat("smoothing offsets\n")
  }
  
  
  smoothedOffset = focalMult(
    x= rasterOffset, 
    w = theFocal$focal,
    filename = paste(tempfile(), '.grd', sep = ''), 
    edgeCorrect=FALSE, 
    cl = NULL
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("smoothing offsets done\n")
  }
  
  names(smoothedOffset) = gsub("[oO]ffset", "", names(smoothedOffset)) 
  offsetStack = writeRaster(
    addLayer(rasterOffset, 
      smoothedOffset[[grep('ones$', names(smoothedOffset), invert=TRUE)]]),
    filename = offsetFile, overwrite = file.exists(offsetFile))  
  
  
  if(endCluster) parallel::stopCluster(theCluster)
  
# create list of partitions
  partitions = as.data.frame(stats::na.omit(raster::unique(rasterFineId)))
  partitions$partition = paste('c', partitions$cellCoarse, 'p', partitions$idCoarse,
    '.', partitions$idFine, sep = '')
  partitions = cbind(ID = 1:nrow(partitions), partitions)
# raster with partition ID's
  
  partitionRaster = raster::calc(rasterFineId, 
    function(x) which(
        x[1]==partitions[,'idCoarse'] &
          x[2] == partitions[,'idFine'] &
          x[3] == partitions[,'cellCoarse']
      )[1])
  levels(partitionRaster) = list(partitions)
  partitionRaster = writeRaster(partitionRaster, filename=idFile, overwrite=file.exists(idFile))
  
  
# offsetMat is the value of the offset at all points in the partition
  
  partitionOffsets = as.data.frame(zonal(
      offsetStack[[grep("^bw", names(offsetStack), invert=TRUE)]],
      partitionRaster,
      'mean', na.rm=TRUE
    ))
  partitionOffsets$partition = partitions[match(
      partitionOffsets$zone, partitions$ID
    ), 'partition']    
  
  
# offsetMat is diagonal matrix with element the integral of offset in the partition 
# different for each CV set
  offsetMat = apply(partitionOffsets[,grep("[oO]ffset", colnames(partitionOffsets))], 2, 
    function(x) {
      res = Matrix::Diagonal(nrow(partitionOffsets), x*prod(res(offsetStack)))
      dimnames(res) = list(partitionOffsets$partition, partitionOffsets$partition) 
      res
    })
  
  regions = gsub("^c[[:digit:]]+p|\\.[[:digit:]]+$", "", partitionOffsets$partition)
  regions = as.integer(regions)
  regionMat = outer(polyCoarse$idCoarse, regions, '==')
  regionMat = Matrix::Matrix(regionMat)
  
  dimnames(regionMat)=
    list(polyCoarse$idCoarse, partitionOffsets$partition)
  
  partitionFreq = raster::freq(partitionRaster)
  rownames(partitionFreq) = raster::levels(partitionRaster)[[1]][partitionFreq[,'value'],'partition']
  partitionAreas = partitionFreq[,'count'][colnames(regionMat)] * prod(raster::res(partitionRaster))  
  
  
  
  result = list(
    rasterCoarse=rasterCoarse,
    polyCoarse = polyCoarse,
    rasterFine=partitionRaster,
    focal=theFocal,
    offset=offsetStack,
    offsetMat = offsetMat,
    regionMat = regionMat,
    partitionAreas=partitionAreas,
    xv = xvMat
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
  
  return(result)
  
}

