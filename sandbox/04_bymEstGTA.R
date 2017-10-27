options(scipen = 22)

#+ thePackages
library('localEM')
library('geostatsp')
library('INLA')
library('diseasemapping')

INLA:::inla.dynload.workaround()
#'

#+ theData
cat(date(), '\n')
cat('loading offset data', '\n')

##GTA files
load('gtaFsaCt2007.RData')

gtaSelectFsa = gtaSelectFsa[substr(gtaSelectFsa$id, 1, 1) == 'M',]

##local-EM parameters
Nboot = 140
threshold = c(1, 1.1, 1.25, 1.5)
iterations = list(tol = 1e-6, maxIter = 2000, gpu = FALSE)
Ncores = 15
path = 'lemLupusV20'

##simulation parameters
Nsim = 20
theRangeM = c(2, 4, 8, 16) * 1000
#'

#+ theSmoothingMat
cat(date(), '\n')
cat('loading smoothing matrix', '\n')

load('lemSetupGTA.RData')

##offsets of spatial polygons of case data based on population data
idCoarse = 1:length(lemSmoothMatGTA$polyCoarse)

offsetRaster = raster::stack(lemSmoothMatGTA$offset$offset, 
                             raster::deratify(lemSmoothMatGTA$rasterFine))
offsetDf = stats::aggregate(x = values(offsetRaster$offset), 
                            by = list(idCoarse = values(offsetRaster$idCoarse)), 
                            FUN = sum)
colnames(offsetDf) = c('idCoarse','offset')
offsetDf = merge(data.frame(idCoarse = idCoarse), offsetDf, by = 'idCoarse', all = TRUE)
offsetDf$offset[is.na(offsetDf$offset)] = 0
offsetDf$expected = offsetDf$offset * prod(res(offsetRaster$offset))
offsetDf$logExpected = log(offsetDf$expected)
offsetDf$logExpected[is.infinite(offsetDf$logExpected)] = NA
offsetDf$id = lemSmoothMatGTA$polyCoarse$id

gtaSelectFsa = merge(gtaSelectFsa[,c('id','FSA')], 
                     offsetDf, 
                     by = 'id')
#'

#+ theBymEst
for(inR in 1:length(theRangeM)) {
  
  cat(date(), '\n')
  cat('loading case data for range: ', theRangeM[inR], '\n')
  
  ##simulated cases
  load(paste('simCasesR', theRangeM[inR], 'mGTA.RData', sep = ''))
  
  gtaSelectFsa = merge(gtaSelectFsa[,c('id','FSA','offset','expected','logExpected')], 
                       mySimCasesGTA, 
                       by.x = 'id', by.y = 'row.names')

  ##risk estimation and exceedance probabilities
  load(paste('lemEstR', theRangeM[inR], 'mGTA.RData', sep = ''))  
  
  cat(date(), '\n')
  cat('computing risk estimation and exceedance probabilities', '\n')

  theBymRiskList = list()
  theBymExcProbList = list()
  for(inS in 1:Nsim) {

    ##raster surface of interest
    theRiskSurface = myLemRiskGTA[[inS]]
    
    ##model fit
    formula = as.formula(paste('count', inS, ' ~ offset(logExpected)', sep = ''))
    bymFitGTA = diseasemapping::bym(formula = formula, data = gtaSelectFsa)

    ##risk estimation
    theBymFit = bymFitGTA$data
    
    bymRiskGTA = raster::rasterize(theBymFit, y = theRiskSurface, field = 'fitted.exp')
    names(bymRiskGTA) = 'fitted.exp'
    
    theBymRiskList[[inS]] = bymRiskGTA

    ##exceedance probabilities
    theBymMarg = bymFitGTA$inla$marginals.fitted.bym
    
    for(inT in 1:length(threshold)) {
      
      theBymFit$excProb = geostatsp::excProb(theBymMarg, threshold = log(threshold[inT]))
      theBymExcProb = raster::rasterize(theBymFit, y = theRiskSurface, field = 'excProb')
      
      if(inT == 1) {
        bymExcProbGTA = theBymExcProb
      } else {
        bymExcProbGTA = raster::addLayer(bymExcProbGTA, theBymExcProb)
      }
    }
    names(bymExcProbGTA) = paste('count', inS, '_threshold', threshold, sep = '')

    theBymExcProbList[[inS]] = bymExcProbGTA
  }
 
  ##output results to external file
  fileRisk = file.path(path, paste('bymRiskR', theRangeM[inR], 'm.grd', sep = ''))
  fileExcProb = file.path(path, paste('bymExcProbR', theRangeM[inR], 'm.grd', sep = ''))

  myBymRiskGTA = raster::writeRaster(
    raster::stack(theBymRiskList), 
    filename = fileRisk, 
    overwrite = file.exists(fileRisk))
  names(myBymRiskGTA) = paste('fitted.exp_count', 1:Nsim, sep = '')
  
  myBymExcProbGTA = raster::writeRaster(
    raster::stack(theBymExcProbList), 
    filename = fileExcProb, 
    overwrite = file.exists(fileExcProb))
  names(myBymExcProbGTA) = paste('count', rep(1:Nsim, each = length(threshold)), 
                                 '_threshold', rep(threshold, Nsim), 
                                 sep = '')

  save(myBymRiskGTA, myBymExcProbGTA, 
       file = paste('bymEstR', theRangeM[inR], 'mGTA.RData', sep = '') 
  )
}  
#'

cat(date(), '\n')
cat('done', '\n')

quit('no')
