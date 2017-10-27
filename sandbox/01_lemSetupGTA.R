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
Sbw = c(seq(0.5, 2.5, by = 0.5), seq(3, 7, by = 1), seq(8, 16, by = 2)) * 1000
Ncores = 15
Nxv = 4
path = 'lemLupusV20'
#'

#+ theSmoothingMat
cat(date(), '\n')
cat('computing rasters and smoothing matrix', '\n')

lemRasterGTA = localEM::rasterPartition(
  polyCoarse = gtaSelectFsa, 
  polyFine = gtaOrigCt, 
  cellsCoarse = cellsCoarse, 
  cellsFine = cellsFine, 
  bw = Sbw, 
  ncores = Ncores, 
  xv = Nxv, 
  path = path, 
  idFile = 'lemId.grd', 
  offsetFile = 'lemOffsets.grd', 
  verbose = TRUE)

lemSmoothMatGTA = localEM::smoothingMatrix(
  rasterObjects = lemRasterGTA, 
  ncores = Ncores, 
  path = path, 
  filename = 'lemSmoothMat.grd', 
  verbose = TRUE)

save(lemSmoothMatGTA, 
     file = 'lemSetupGTA.RData'
     )
#'

cat(date(), '\n')
cat('done', '\n')

quit('no')
