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
Sbw = c(seq(0.5, 2.5, by = 0.5), seq(3, 7, by = 1), seq(8, 16, by = 2)) * 1000
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

#+ theLemBoot
cat(date(), '\n')
cat('generating bootstrap cases for input thresholds\n')

# offsets of spatial polygons of case data based on population data
idCoarse = 1:length(lemSmoothMatGTA$polyCoarse)
  
offsetRaster = raster::stack(lemSmoothMatGTA$offset$offset, 
  raster::deratify(lemSmoothMatGTA$rasterFine))
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

theBootCounts = bootCountsDf

# estimate risk from bootstrap cases
theBootRiskList = list()
for(inB in 1:length(Sbw)) {
  
	cat(date(), '\n')
	cat('running local-EM estimation for bootstrap cases for bw: ', Sbw[inB], '\n')
      
	bootLemRisk = localEM::riskEst(
    cases = bootCountsDf, 
    lemObjects = lemSmoothMatGTA, 
    bw = Sbw[inB], 
    ncores = min(Ncores, 7), 
    iterations = iterations, 
    path = path, 
    filename = paste('lemBootBw', Sbw[inB], '.grd', sep = ''), 
    verbose = FALSE)
  bootEstRisk = bootLemRisk$riskEst
      
  theBootRiskList[[inB]] = bootEstRisk

	if(inB %% 4 == 0) {
		save(Nboot, threshold, iterations, 
		     theBootCounts, theBootRiskList, 
		     file = 'lemBootGTA.RData'
		)
	}

}
names(theBootRiskList) = paste('bw', Sbw, sep = '')

save(Nboot, threshold, iterations, 
     theBootCounts, theBootRiskList, 
     file = 'lemBootGTA.RData'
)
#'

cat(date(), '\n')
cat('done', '\n')

quit('no')
