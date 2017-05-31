options(scipen = 22)

#+ thePackages
# devtools::install("../pkg/localem")
library('localEM')

source('../pkg/localem/R/lemXvFun.R')
if(Sys.info()['user'] == 'patrick') {
	cellsFine = 300
	ncores = 8
} else {
	ncores = 2
	cellsFine = 60
}
#'

#+ theData
data(kentuckyCounty)
data(kentuckyTract)

kentuckyTract$expected = (3000 / sum(kentuckyTract$expected, na.rm = TRUE)) * kentuckyTract$expected

bandwidth = seq(5, 40, by = 5) * 1000
#'

#+ theSmoothingMat
lemRasterKentucky = rasterPartition(
  polyCoarse = kentuckyCounty, 
  polyFine = kentuckyTract, 
  cellsCoarse = 10, 
  cellsFine = 120, 
  bw = bandwidth, 
  ncores = 8, 
  idFile = 'id.grd', 
  offsetFile = 'off.grd', 
  verbose = FALSE)

lemSmoothMatKentucky = smoothingMatrix(
  rasterObjects = lemRasterKentucky, 
  ncores = 2, 
  verbose = FALSE)
#'

#+ theSim
logOffsetKentucky = log(lemRasterKentucky$offset$offset)
names(logOffsetKentucky) = 'logOffset'

set.seed(0)
simCasesKentucky = geostatsp::simLgcp(
  param = c(mean = 0, variance = 0.4^2, 
    range = 120 * 1000, shape = 2),
  covariates = logOffsetKentucky, 
  offset = 'logOffset')

simCounts = table(over(simCasesKentucky$events, kentuckyCounty)$id)
kentuckyCounty$count = simCounts[kentuckyCounty$id]
kentuckyCounty$count[is.na(kentuckyCounty$count)] = 0

simSpKentucky = simCasesKentucky$events
simSurfaceKentucky = simCasesKentucky$raster[['relativeIntensity']]
#'


#+ theXv
lemCvOldKentucky = lemXv(
  x = kentuckyCounty, 
  lemObjects = lemSmoothMatKentucky, 
  Nxv = 5, 
  ncores = 2, 
  verbose = FALSE)
    
bestBwOldKentucky = lemCvOldKentucky$bw[which.min(lemCvOldKentucky$cv)]

lemCvNewKentucky = lemXvNew(
  polyCoarse = kentuckyCounty, 
  polyFine = kentuckyTract, 
  cellsCoarse = 10, 
  cellsFine = 120, 
  bw = bandwidth, 
  Nxv = 4, 
  ncores = 2, 
  tol = 1e-6, 
  maxIter = 2000, 
  verbose = FALSE)
    
bestBwNewKentucky = lemCvNewKentucky$bw[which.min(lemCvNewKentucky$cv)]
#'

#+ theLemRisk
for(inBw in 1:length(bandwidth)) {

  lemEstKentucky = riskEst(
    x = kentuckyCounty, 
    lemObjects = lemSmoothMatKentucky, 
    bw = bandwidth[inBw], 
    ncores = 1)
  
  if(inBw == 1) {
    lemRiskKentucky = lemEstKentucky
    
  } else {
    lemRiskKentucky = addLayer(lemRiskKentucky, lemEstKentucky)
    
  }
}
#'



