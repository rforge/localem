options(scipen = 22)

#+ thePackages
library('localEM')
library('geostatsp')
#'

#+ theData
cat(date(), '\n')
cat('loading offset data', '\n')

##GTA files
load('gtaFsaCt2007.RData')

gtaSelectFsa = gtaSelectFsa[substr(gtaSelectFsa$id, 1, 1) == 'M',]

##local-EM parameters
cellsCoarse = 15
cellsFine = 165

##simulation parameters
cellsSimulate = 4 * cellsFine
Nsim = 20
theRangeM = c(2, 4, 8, 16) * 1000

##log-offsets
theLogExpectedRaster = geostatsp::spdfToBrick(
  gtaOrigCt, 
  geostatsp::squareRaster(gtaOrigCt, cellsSimulate), 
  pattern = '^expected$',
  logSumExpected = TRUE)
theLogExpectedRaster = raster::trim(raster::mask(theLogExpectedRaster, gtaSelectFsa))
#'

#+ theSim
cat(date(), '\n')
cat('running simulations', '\n')

for(inR in 1:length(theRangeM)) {
  
  cat(date(), '\n')
  cat('simulating case data for range: ', theRangeM[inR], '\n')
  
  ##model parameters
  set.seed(201710)
  theModel = c(mean = 0, variance = 0.4^2, range = theRangeM[inR], shape = 2)

  ##simulate surfaces and cases from LGCP
  theSimData = geostatsp::simLgcp(
    param = theModel, 
    covariates = list(logExpected = theLogExpectedRaster), 
    offset = 'logExpected', 
    n = Nsim)
  
  ##relative intensity surfaces
  theRasterNames = names(theSimData$raster)
  
  mySimSurfaceGTA = theSimData$raster[[grep('^relativeIntensity', theRasterNames)]]
  
  ##case locations
  mySimSpGTA = theSimData[grep('^events[[:digit:]]+?', names(theSimData))]

  ##aggregated cases
  theSimCases = lapply(
    theSimData[grep('^events[[:digit:]]+?', names(theSimData))], 
    function(qq) over(qq, gtaSelectFsa)[,'id'])
  
  mySimCasesGTA = as.data.frame(lapply(
    theSimCases, 
    function(xx) as.vector(table(xx, exclude = NULL)[as.character(gtaSelectFsa$id)])))
  mySimCasesGTA[is.na(mySimCasesGTA)] = 0
  colnames(mySimCasesGTA) = gsub('^events', 'count', colnames(mySimCasesGTA))
  rownames(mySimCasesGTA) = as.character(gtaSelectFsa$id)
  
  ##log-offsets
  myLogExpectedGTA = theSimData$raster[['logExpected']]
  
  save(myLogExpectedGTA, 
       mySimSurfaceGTA, mySimSpGTA, mySimCasesGTA, 
       file = paste('simCasesR', theRangeM[inR], 'mGTA.RData', sep = '')
  )
}
#'

cat(date(), '\n')
cat('done\n')

quit('no')
