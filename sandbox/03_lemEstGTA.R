options(scipen = 22)

#+ thePackages
library('localEM')
library('geostatsp')
#'

#+ theAddOns
source('excProbBootFun.R')
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
#'

#+ theLemEst
for(inR in 1:length(theRangeM)) {
  
  cat(date(), '\n')
  cat('loading case data for range: ', theRangeM[inR], '\n')
  
  ##simulated cases
  load(paste('simCasesR', theRangeM[inR], 'mGTA.RData', sep = ''))

  gtaSelectFsa = merge(gtaSelectFsa[,c('id','FSA')], 
                       mySimCasesGTA, 
                       by.x = 'id', by.y = 'row.names')
  
  cat(date(), '\n')
  cat('computing CV scores', '\n')

  ##CV scores
  lemXvGTA = localEM::lemXv(
    cases = gtaSelectFsa@data[,grep('^count[[:digit:]]', names(gtaSelectFsa))], 
    lemObjects = lemSmoothMatGTA, 
    ncores = Ncores, 
    iterations = iterations, 
    randomSeed = theRangeM[inR], 
    path = path, 
    filename = paste('lemRiskR', theRangeM[inR], 'm.grd', sep = ''), 
    verbose = TRUE)

  myLemXvGTA = lemXvGTA$xv
  myLemRiskGTA = lemXvGTA$riskEst
  
  cat(date(), '\n')
  cat('computing exceedance probabilities', '\n')
  
  ##risk estimation and exceedance probabilities
  # excProbGTA = localEM::excProb(
  #   lemObjects = lemXvGTA, 
  #   threshold = threshold, 
  #   Nboot = Nboot, 
  #   ncores = min(Ncores, 7), 
  #   iterations = iterations, 
  #   path = path, 
  #   filename = paste('lemExcProbR', theRangeM[inR], 'm.grd', sep = ''), 
  #   verbose = TRUE)

  excProbGTA = excProbBoot(
    lemObjects = lemXvGTA, 
    lemBootData = 'lemBootGTA.RData', 
    path = path, 
    filename = paste('lemExcProbR', theRangeM[inR], 'm.grd', sep = ''), 
    verbose = TRUE)
  
  myLemExcProbGTA = excProbGTA$excProb

  save(myLemXvGTA, myLemRiskGTA, myLemExcProbGTA, 
       file = paste('lemEstR', theRangeM[inR], 'mGTA.RData', sep = '')
      )
}
#'

cat(date(), '\n')
cat('done', '\n')

quit('no')
