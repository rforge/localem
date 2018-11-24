#' @title Computes the exceedance probabilities
#'
#' @description The \code{excProb} function first bootstraps cases with the input risk thresholds and expected counts from the rasterization of the spatial polygons of population data, and then, computes the exceedance probabilities with the same bandwidth as the risk estimation.
#'
#' @param lemObjects List of arrays for the smoothing matrix, and raster stacks for the partition, smoothed offsets and risk estimation
#' @param threshold Vector of risk thresholds
#' @param Nboot Number of bootstraps
#' @param bw Bandwidth for smoothing bootstrap samples
#' @param fact Aggregation factor prior to 'final step' smoothing (set to zero to skip final step)
#' @param ncores Number of cores/threads for parallel processing
#' @param iterations List of convergence tolerance, number of iterations, and use of gpuR package for running local-EM recursions
#' @param path Folder for storing rasters
#' @param filename Filename (must have .grd extension) of the exceedance probabilities
#' @param verbose Verbose output
#'
#' @details After using the \code{excProb} function, the exceedance probabilities are computed on a fine resolution based on the rasterization of the spatial polygons of population data.
#'
#' @return The \code{excProb} function returns a raster brick of exceedance probabilities of input risk thresholds.
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
#' threshold = c(1, 1.1, 1.25, 1.5)
#' Nboot = 100
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
#'
#' # risk estimation for bandwidth of 10km
#' lemRisk = riskEst(cases = kentuckyCounty[,c('id','count')],
#'						lemObjects = lemSmoothMat,
#'						bw = bw[1],
#'						fact = 1,
#'						ncores = ncores,
#'						path = path,
#'						filename = 'lemRisk.grd',
#'						verbose = TRUE)
#'
#' # exceedance probabilities for bandwidth of 10km
#' lemExcProb = excProb(lemObjects = lemRisk,
#'					  	threshold = threshold,
#'					  	Nboot = Nboot,
#'						bw = bw[1],
#'						fact = 1,
#'              		ncores = ncores,
#'              		path = path,
#'					  	filename = 'lemExcProb.grd',
#'              		verbose = TRUE)
#'
#' # plot exceedance probabilities
#' pCol = mapmisc::colourScale(lemExcProb$excProb,
#'						breaks = c(0,0.2,0.8,0.95,1),
#'						style = 'fixed', dec = 2,
#'						col = c('green','yellow','orange','red'))
#'
#' par(mfrow = c(2,2), mar = c(0.5,0.5,3,0.5))
#' for(inT in 1:length(threshold)) {
#'		plot(lemExcProb$excProb[[inT]],
#'			main = paste('Exc Prob, t=', threshold[inT], ' (bw=10km)', sep = ''),
#'			col = pCol$col, breaks = pCol$breaks,
#'			axes = FALSE, box = FALSE, legend = FALSE,
#'			add = FALSE)
#' }
#' mapmisc::legendBreaks('topright', pCol)
#' }
#'
#' @export
excProb = function(
    lemObjects,
    threshold = 1,
    Nboot = 100,
    bw,
  	fact = 1,
	ncores = 1,
	iterations = list(tol = 1e-5, maxIter = 1000, gpu = FALSE),
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
  if(missing(bw)) {
    bw = as.numeric(gsub("^bw|_.*", "", names(lemObjects$riskEst)[1]))
  }
  bwString = grep(
    paste0("^(bw)?",format(bw[1],scientific=FALSE),"$"),
    names(lemObjects$smoothingMatrix$smoothingArray),
    value=TRUE)[1]
  bw = as.numeric(gsub('^bw|xv[[:digit:]]+', '', bwString))

  # risk estimate of interest
  theEstRisk = lemObjects$riskEst
  theEstNames = names(theEstRisk)

  if(verbose) {
    cat(date(), "\n")
    cat("generating bootstrap cases for input thresholds\n")
  }

  # offsets of spatial polygons of case data based on population data
  idCoarse = 1:length(lemObjects$smoothingMatrix$polyCoarse)

  offsetRaster = raster::stack(lemObjects$smoothingMatrix$offset$offset,
      raster::deratify(lemObjects$smoothingMatrix$rasterFine))
  offsetDf = stats::aggregate(x = values(offsetRaster$offset) * prod(res(offsetRaster)),
      by = list(idCoarse = values(offsetRaster$idCoarse)),
      FUN = sum)
  colnames(offsetDf) = c('idCoarse','offset')
  offsetDf = merge(data.frame(idCoarse = idCoarse), offsetDf, by = 'idCoarse', all = TRUE)
  offsetDf$offset[is.na(offsetDf$offset)] = 0

  # bootstrap cases
  offsetT = outer(offsetDf$offset, threshold)

  bootCountsDf = matrix(
      data = stats::rpois(
          length(offsetT) * Nboot,
          rep(offsetT, Nboot)),
      nrow = nrow(offsetDf),
      ncol = length(threshold) * Nboot,
      dimnames = list(rownames(offsetDf),
          paste('count', rep(1:Nboot, rep(length(threshold), Nboot)),
              '_threshold', rep(threshold, Nboot),
              sep = '')
      )
  )

  # estimate risk from bootstrap cases
  if(verbose) {
    cat(date(), "\n")
    cat("running local-EM estimation for bootstrap cases\n")
  }

  theBootRiskList = list()
#  for(inB in 1:length(bw))
  inB = 1
#  {

#	indexBw = which(bw[1:inB] %in% bw[inB])

	# generate results for input bw
 #   if(length(indexBw) == 1) {

      bootLemRisk = riskEst(
			cases = bootCountsDf,
            lemObjects = lemObjects,
            bw = bw[inB],
			fact=fact,
			ncores=ncores, 
			iterations=iterations, 
            path = path,
            filename = tempfile(
				paste0('riskBoot', bwString),
				path, '.grd'),
            verbose = verbose)


        bootEstRisk = bootLemRisk$riskEst

#        theBootRiskList[[inB]] = bootEstRisk

    # use previous results if same bw was used before
 #   } else {
 #      theBootRiskList[[inB]] = theBootRiskList[[indexBw[1]]]
 #   }
 # }

    # # first bw
    # if(inB == 1) {

      # bootLemRisk = riskEst(
          # cases = bootCountsDf,
          # lemObjects = lemObjects$smoothingMatrix,
          # bw = bw[inB],
          # ncores = ncores,
          # iterations = iterations,
          # path = path,
          # filename = paste('riskBootTempBw', bw[inB], '.grd', sep = ''),
          # verbose = verbose)
      # bootEstRisk = bootLemRisk$riskEst

      # theBootRiskList[[inB]] = bootEstRisk

      # # additional bw
    # } else {

      # indexBw = which(bw[1:inB] %in% bw[inB])

      # if(length(indexBw) == 1){

        # bootLemRisk = riskEst(
            # cases = bootCountsDf,
            # lemObjects = lemObjects$smoothingMatrix,
            # bw = bw[inB],
            # ncores = ncores,
            # iterations = iterations,
            # path = path,
            # filename = paste('riskBootTempBw', bw[inB], '.grd', sep = ''),
            # verbose = FALSE)
        # bootEstRisk = bootLemRisk$riskEst

        # theBootRiskList[[inB]] = bootEstRisk

        # # use previous results if same bw was used before
      # } else {
        # theBootRiskList[[inB]] = theBootRiskList[[indexBw[1]]]
      # }
    # }
  # }

  # exceedance probabilities
  if(verbose) {
    cat(date(), "\n")
    cat("computing exceedance probabilities with input thresholds\n")
  }

#  theExcProbList = list()
#  for(inB in 1:length(bw)) {

 #   theBootEstRisk = theBootRiskList[[inB]]
#stuff1 <<- bootEstRisk
#stuff2 <<- inB

bwHere = grep(
  paste0(bwString, '_'),
  names(lemObjects$riskEst), value=TRUE)

if(!length(bwHere))
  warning("can't find risk estimate for bandwidth supplied")
# TO DO:  compute risk estimate for requested bandwidth


Sthreshold = factor(
  gsub("bw[[:digit:]]+_|count[[:digit:]]+_", "", names(bootEstRisk))
  )
thresholdMat = outer(
  Sthreshold, levels(Sthreshold), '=='
  ) / Nboot
colnames(thresholdMat) = levels(Sthreshold)

theExcProbList = parallel::mcmapply(
  function(Dy) {
    raster::overlay(
        x = bootEstRisk,
        y = lemObjects$riskEst[[Dy]],
        fun = function(x,y) {
          (x<y) %*% thresholdMat
#          bob <<- x<y
#        apply(x < y,1,function(xx) {
#            tapply(xx, Sthreshold, mean)
          },
        filename = tempfile(
          paste0('probBoot', Dy),
          tmpdir=tempdir(), fileext='.grd')
        )
  },
  Dy = bwHere,
  mc.cores=ncores,
  SIMPLIFY=FALSE
  )

#    indexT = gsub('count[[:digit:]]+_', '', names(theBootEstRisk))

 #   theExcProb = raster::stackApply(theProbEstRisk,
 #       indices = indexT,
 #       fun = mean,
 #       filename = paste(tempfile('excProbTemp', path), '.grd', sep = ''),
 #       overwrite = TRUE)

#    theExcProbList[[inB]] = theExcProb
#  }
theExcProbStack = raster::stack(theExcProbList)
names(theExcProbStack) = as.vector(t(outer(
  names(theExcProbList), levels(Sthreshold), paste, sep='_'
  )))

  theExcProbBrick = raster::writeRaster(
      theExcProbStack,
      filename = filename,
      overwrite = file.exists(filename))


  result = list(
      bootEst = bootEstRisk,
      excProb = theExcProbBrick
  )

  if(verbose) {
    cat(date(), "\n")
    cat("done\n")
  }

  # remove temporary raster files
#  unlink(file.path(path, 'riskBootTempBw*'))
#  unlink(file.path(path, 'probBootTemp*'))
#  unlink(file.path(path, 'excProbTemp*'))

  return(result)
}


#' @export
excProbMulti = function(
  lemObjects,
  threshold = 1,
  Nboot = 100,
  bw = lemObjects$bw[1],
  est,
  fact = 1,
  ncores = 1,
  iterations = list(tol = 1e-5, maxIter = 1000, gpu = FALSE),
  path = getwd(),
  filename = 'lemExcProb.grd',
  verbose = FALSE
  ){

  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  
  # warning messages
  if(missing(lemObjects)) {
    stop('smoothing matrix and rasters not supplied')
  }

  if(!length(grep('/', filename))) {
    filename = file.path(path, filename)
  }
  if(!(length(grep("\\.gr[id]$", filename)))){
    warning("filename should have .grd extension")
  }
  
  # bandwidths of interest
  bwString = grep(
    paste0('^(bw)?', bw[1], '$'), 
    names(lemObjects$smoothingMatrix[[1]]$smoothingMatrix$smoothingArray), 
    value=TRUE)[1]
  bw = as.numeric(gsub('^bw|xv[[:digit:]]+', '', bwString))

  # risk estimate of interest
  if(missing(est)) {
    theEstRisk = lemObjects$riskEst[[
      grep(paste0('^', bwString, '(_|$)'), names(lemObjects$riskEst))
    ]]
  } else {
    theEstRisk = est
  }
  theEstNames = names(theEstRisk)

  # bootstrap cases
  if(verbose) {
    cat(date(), '\n')
    cat('generating bootstrap cases for input thresholds\n')
  }

  bootCountsList = list()
  for(inM in 1:length(lemObjects$smoothingMatrix)) {
    
  	lemBootMap = lemObjects$smoothingMatrix[[inM]]
  	
  	# offsets
  	idCoarse = 1:length(lemBootMap$smoothingMatrix$polyCoarse)
  
  	rasterFine = lemBootMap$smoothingMatrix$rasterFine
  	levels(rasterFine) = data.frame(levels(rasterFine)[[1]], 
  	                                offset = diag(lemBootMap$smoothingMatrix$offsetMat$offset))

    offsetRaster = raster::deratify(rasterFine)
    offsetDf = stats::aggregate(x = values(offsetRaster$offset) * prod(res(offsetRaster)), 
      by = list(idCoarse = values(offsetRaster$idCoarse)), 
      FUN = sum)
    colnames(offsetDf) = c('idCoarse','offset')
    offsetDf = merge(data.frame(idCoarse = idCoarse), offsetDf, by = 'idCoarse', all = TRUE)
    offsetDf$offset[is.na(offsetDf$offset)] = 0
  
    # bootstrap cases
    offsetT = outer(offsetDf$offset, threshold)

    bootCountsDf = matrix(
      data = stats::rpois(
        length(offsetT) * Nboot, 
        rep(offsetT, Nboot)), 
      nrow = nrow(offsetDf), 
      ncol = length(threshold) * Nboot, 
      dimnames = list(rownames(offsetDf), 
        paste('count', rep(1:Nboot, rep(length(threshold), Nboot)), 
          '_threshold', rep(threshold, Nboot),
          sep = '')
        )
      )
    
    bootCountsList[[inM]] = bootCountsDf
  }

  # estimate risk from bootstrap cases
  if(verbose) {
    cat(date(), '\n')
    cat('running local-EM estimation for bootstrap cases\n')
  }
  
  bootLemRisk = riskEstMulti(
    cases = bootCountsList, 
    lemObjects = lemObjects$smoothingMatrix, 
    bw = bw, 
    fact = fact, 
    ncores = ncores, 
	iterations = iterations, 
    path = path, 
    filename = tempfile(paste0('riskBoot', bwString), path, '.grd'), 
    verbose = verbose)
  
  bootEstRisk = bootLemRisk$riskEst
  
  # exceedance probabilities
  if(verbose) {
    cat(date(), '\n')
    cat('computing exceedance probabilities with input thresholds\n')
  }
  
  Fthreshold = factor(
    gsub('bw[[:digit:]]+_|count[[:digit:]]+_', '', names(bootEstRisk))
  )
  
  thresholdMat = outer(Fthreshold, levels(Fthreshold), '==') / Nboot
  colnames(thresholdMat) = levels(Fthreshold)
  
  theExcProbList = parallel::mcmapply(
    function(Dy) {
      raster::overlay(
        x = bootEstRisk,
        y = theEstRisk[[Dy]],
        fun = function(x,y) {
          (x < y) %*% thresholdMat
        },
        filename = tempfile(paste0('probBoot', Dy), tmpdir = tempdir(), fileext = '.grd')
      )
    }, 
    Dy = names(theEstRisk), 
    mc.cores = ncores, 
    SIMPLIFY = FALSE
  )
  
  theExcProbStack = raster::stack(theExcProbList)
  names(theExcProbStack) = as.vector(t(outer(
    names(theExcProbList), levels(Fthreshold), paste, sep = '_'
  )))
  
  theExcProbBrick = raster::writeRaster(
    theExcProbStack, 
    filename = filename, 
    overwrite = file.exists(filename))
  
  theFiles = unlist(lapply(theExcProbList, raster::filename))
  if(any(nchar(theFiles))) {
    unlink(gsub('grd$', 'gri', theFiles))
  }
  
  result = list(
    bootEst = bootEstRisk, 
    excProb = theExcProbBrick
  )
  
  if(verbose) {
    cat(date(), '\n')
    cat('done\n')
  }
  
  return(result)
}
