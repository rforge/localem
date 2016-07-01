library('localEM')

data('kentuckyCounty') 
data('kentuckyTract')


lemRaster = rasterPartition(
		polyCoarse = kentuckyCounty, 
		polyFine = kentuckyTract, 
    cellsCoarse = 20, cellsFine = 200, 
    bw = c(10, 20) * 1000, 
    ncores = 2, 
    idFile = 'id.grd', 
		offsetFile = 'offset.grd'
)
