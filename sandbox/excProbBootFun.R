#' Computes the exceedance probabilities with pre-created risk estimation for input parameters
excProbBoot = function(
  lemObjects, 
  lemBootData, 
  path = getwd(), 
  filename = 'lemExcProb.grd', 
  verbose = FALSE
){
  
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  # if(missing(filename)) {
  # filename = paste(tempfile('lemExcProb', path), '.grd', sep = '')
  # }
  if(!length(grep('/', filename))) {
    filename = file.path(path, filename)
  }
  if(!(length(grep("\\.gr[id]$", filename)))){
    warning("filename should have .grd extension")
  }
  
  # warning messages
  if(missing(lemObjects)) {
    stop("smoothing matrix and rasters not supplied")
  }
  
  # bandwidths of interest
  bwString = lemObjects$bw
  bw = as.numeric(gsub('^bw', '', bwString))
  
  # risk estimate of interest
  theEstRisk = lemObjects$riskEst
  theEstNames = names(theEstRisk)
  
  # exceedance probabilities
  if(verbose) {
    cat(date(), "\n")
    cat("computing exceedance probabilities with input thresholds\n")
  }
  
  # estimated risk from bootstrap cases
  load(file = lemBootData)
  
  theExcProbList = list()
  for(inB in 1:length(bw)) {
    
    theBootEstRisk = theBootRiskList[[bwString[inB]]]
  
    theProbEstRisk = raster::overlay(x = theBootEstRisk, y = theEstRisk[[inB]], 
      fun = function(x,y) return(x < y), 
      filename = paste(tempfile('probBootTemp', path), '.grd', sep = ''), 
      overwrite = TRUE)
      
    indexT = gsub('count[[:digit:]]+_', '', names(theBootEstRisk))
      
    theExcProb = raster::stackApply(theProbEstRisk, 
      indices = indexT, 
      fun = mean, 
      filename = paste(tempfile('excProbTemp', path), '.grd', sep = ''), 
      overwrite = TRUE)
      
      theExcProbList[[inB]] = theExcProb
  }
  
  theExcProbStack = raster::writeRaster(
    raster::stack(theExcProbList), 
    filename = filename, 
    overwrite = file.exists(filename))
  names(theExcProbStack) = paste(rep(theEstNames, each = length(threshold)), 
    '_threshold', rep(threshold, length(theEstNames)), 
    sep = '')
  
  result = list(
    riskEst = theEstRisk, 
    excProb = theExcProbStack
  )
  
  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }	
  
  # remove temporary raster files
  unlink(file.path(path, 'probBootTemp*'))
  unlink(file.path(path, 'excProbTemp*'))
  
  return(result)
}
