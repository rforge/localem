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

##specificity parameter
theSpec = seq(0, 1, by = 1/Nboot)
#'

#+ theRoc
for(inR in 1:length(theRangeM)) {
  
  cat(date(), '\n')
  cat('computing ROC for range: ', theRangeM[inR], '\n')
  
  ##simulation surfaces and cases
  load(paste('simCasesR', theRangeM[inR], 'mGTA.RData', sep = ''))
  
  ##risk estimation and exceedance probabilities
  load(paste('lemEstR', theRangeM[inR], 'mGTA.RData', sep = ''))

  load(paste('bymEstR', theRangeM[inR], 'mGTA.RData', sep = ''))
  
  myLemRocGTA = list()
  myBymRocGTA = list()
  for(inS in 1:Nsim) {

    ##simulation surface
    theSimSurface = mySimSurfaceGTA[[inS]]
    names(theSimSurface) = 'relativeIntensity'
    
    ##ROC from local-EM analysis
    lemExcProbGTA = 
      myLemExcProbGTA[[grep(paste('(count|case)', inS, '_', sep = ''), names(myLemExcProbGTA))]]
    names(lemExcProbGTA) = 
      gsub('^bw[[:digit:]]+\\_(count|case)[[:digit:]]+\\_', '', names(lemExcProbGTA))
    names(lemExcProbGTA) = gsub('^threshold', 'threshold.', names(lemExcProbGTA))

    lemRocGTA = try(
      geostatsp::spatialRoc(
        fit = lemExcProbGTA, 
        rr = threshold, 
        truth = theSimSurface, 
        random = FALSE, 
        spec = theSpec), 
      silent = TRUE)
    
  	if(class(lemRocGTA) != 'try-error') {
  	  myLemRocGTA[[inS]] = lemRocGTA[,]
  	}

    ##ROC from BYM analysis
    bymExcProbGTA = 
      myBymExcProbGTA[[grep(paste('(count|case)', inS, '_', sep = ''), names(myBymExcProbGTA))]]
    names(bymExcProbGTA) = 
      gsub('^(count|case)[[:digit:]]+\\_', '', names(bymExcProbGTA))
    names(bymExcProbGTA) = gsub('^threshold', 'threshold.', names(bymExcProbGTA))
    
    bymRocGTA = try(
      geostatsp::spatialRoc(
        fit = bymExcProbGTA, 
        rr = threshold, 
        truth = theSimSurface, 
        random = FALSE, 
        spec = theSpec), 
      silent = TRUE)
    
  	if(class(bymRocGTA) != 'try-error') {
  	  myBymRocGTA[[inS]] = bymRocGTA[,]
  	}
  }
  
  save(myLemRocGTA, myBymRocGTA, 
       file = paste('rocR', theRangeM[inR], 'mGTA.RData', sep = '') 
  )
}  
#'

cat(date(), '\n')
cat('done', '\n')

quit('no')
