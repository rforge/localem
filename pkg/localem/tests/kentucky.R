library('localEM')

data('kentuckyCounty') 
data('kentuckyTract')

ncores = 1+ (.Platform$OS.type=='unix')

lemRaster = rasterPartition(
		polyCoarse = kentuckyCounty, 
		polyFine = kentuckyTract, 
    cellsCoarse = 20, cellsFine = 200, 
    bw = c(10, 20) * 1000, 
    ncores = ncores, 
    idFile = 'id.grd', 
		offsetFile = 'offset.grd'
)


lemSmoothMat = smoothingMatrix(
		rasterObjects = lemRaster, 
    ncores = ncores) 
		
