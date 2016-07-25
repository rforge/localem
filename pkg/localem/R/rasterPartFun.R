#' @title Generates the partitions by overlaying rasters of the coarse and fine polygons
#'
#' @description The \code{rasterPartition} function first rasterizes the coarse and fine spatial polygons based on their respective input resolutions, and then, 
#'  overlays these rasters to generate the raster of partitions of the local-EM algorithm. 
#'  It also applies the kernel smoothing function with input bandwidths to the expected counts of the fine polygons to obtain the smoothed offsets (i.e., smoothed expected counts / cell area) of the partitions. 
#' 
#' @param polyCoarse Spatial polygons of case data
#' @param polyFine Spatial polygons of population data
#' @param bw Vector of bandwidths
#' @param cellsCoarse Horizontal/vertical resolution of raster applied to coarse polygons
#' @param cellsFine Horizontal/vertical resolution of raster applied to fine polygons
#' @param ncores Number of cores/threads for parallel processing
#' @param idFile Filename (must have .grd extension) of the raster of partitions
#' @param offsetFile Filename (must have .grd extension) of the rasters of smoothed offsets
#' @param verbose verbose output
#' 
#' 
#' @details After using the \code{rasterPartition} function, the fine raster is a raster stack containing the IDs for the partitions created by overlaying the coarse and fine rasters. 
#'  The offset raster is a raster stack containing the offsets of the partitions smoothed with the specified bandwidths. These values represent the denominator of the kernel smoothing matrix. 
#'
#'@return The \code{rasterPartition} function returns a list containing the raster of the coarse polygons, 
#'  raster stacks of the partitions and offsets, 
#'  focal weight matrix, and 
#'  the coarse spatial polygons. 
#'  
#' @examples 
#' data(kentuckyCounty)
#' 
#' data(kentuckyTract)
#' 
#' \dontrun{
#' lemRaster = rasterPartition(polyCoarse = kentuckyCounty, polyFine = kentuckyTract, 
#'                    cellsCoarse = 40, cellsFine = 400, 
#'                    bw = c(10, 15, 20, 25) * 1000, 
#'                    ncores = 4, 
#'                    idFile = 'id.grd', offsetFile = 'offset.grd')
#'}
#'
#' @export
rasterPartition = function(
  	polyCoarse, 
  	polyFine, 
  	cellsCoarse = 40, 
  	cellsFine = 400, 
  	bw, 
  	ncores = 2, 
  	idFile = paste(tempfile(), 'Id.grd', sep = ''), 
  	offsetFile = paste(tempfile(), 'Offset.grd', sep = ''), 
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
	
  rasterCoarse=geostatsp::squareRaster(polyCoarse, cellsCoarse)
  values(rasterCoarse) = seq(1, ncell(rasterCoarse))
  names(rasterCoarse) = "cellCoarse"
	
  rasterFine = disaggregate(rasterCoarse,
      ceiling(cellsFine/ncol(rasterCoarse)))
	names(rasterFine) = 'cellCoarse'
	
	# coarse poly ID's for fine raster
	idCoarse = 1:length(polyCoarse)
  names(idCoarse) = polyCoarse$id
  polyCoarse$idCoarse = idCoarse
	
	rasterIdCoarse = rasterize(
			polyCoarse,
      rasterFine,
      field='idCoarse')
	names(rasterIdCoarse) = 'idCoarse'
	
  polyFine@data[is.na(polyFine@data[,'expected']),'expected'] = 0
	
	rasterOffset = geostatsp::spdfToBrick(
			x=polyFine,
      template=rasterFine,
      pattern='^expected$')
  names(rasterOffset) = 'offset'
	
  rasterOffset = raster::mask(rasterOffset, rasterIdCoarse[['idCoarse']])
	
  rasterOffset = writeRaster(rasterOffset,
      paste(tempfile(), ".grd"))
	
	
	
	
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
	
	rasterFineId = brick(rasterIdCoarse, rasterIdFine, rasterFine)	
	
	
  rasterFineId = writeRaster(rasterFineId, idFile,
      overwrite=file.exists(idFile))
	
	
  theFocal = focalFromBw(bw = bw, rasterFine, ncores=ncores)
	
  Sagg = sort(setdiff(unique(theFocal$bw$fact), 1))
  offsetAgg = parallel::mcmapply(
    	aggregate,
    	fact = Sagg,
    	MoreArgs=list(
      		x=rasterOffset[['offset']])
  )
  for(D in names(offsetAgg)){
    offsetAgg[[D]] = addLayer(offsetAgg[[D]], raster(offsetAgg[[D]]))
    names(offsetAgg[[D]]) = 'offset'
    offsetAgg[[D]] = writeRaster(offsetAgg[[D]],
        paste(tempfile(), ".grd"))
  }
	
  forSmooth = expand.grid(bw=names(theFocal$focal),
      layer='offset')
	
  smoothedOffset = parallel::mcmapply(
    	Dbw=forSmooth[,'bw'],
    	Dlayer = forSmooth[,'layer'],
    	function(Dbw,Dlayer){
      	res = raster::focal(
        		x=rasterOffset[[Dlayer]],
        		w=theFocal$focal[[Dbw]],
        		na.rm=TRUE,pad=TRUE)
      	names(res) = paste(Dlayer,Dbw,sep='.')
      	res
    	},
    	mc.cores=ncores)
	
  offsetStack = rasterOffset
  for(D in 1:length(smoothedOffset))
    offsetStack = addLayer(offsetStack, smoothedOffset[[D]])
	
  theNames = names(offsetStack)
  offsetStack = raster::mask(offsetStack, rasterFineId[['idCoarse']])
  offsetStack = writeRaster(offsetStack, offsetFile,overwrite=TRUE)
  names(offsetStack) = theNames
	
  result = list(
    	rasterCoarse=rasterCoarse,
    	rasterFine=rasterFineId,
    	focal=theFocal,
    	offset=offsetStack,
    	polyCoarse = polyCoarse[,grep("^id", names(polyCoarse))]
  )
  
  return(result)
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }
}
