library('localEM')
library('mapmisc')

data('kentuckyCounty') 
data('kentuckyTract')

if(interactive()) {
	kMapOrig = openmap(
			kentuckyCounty,
			path='stamen-toner')
	kMap = tonerToTrans(kMapOrig)
} else {
	data('kMap')
}

ncores = 1+ (.Platform$OS.type=='unix')

lemRaster = rasterPartition(
		polyCoarse = kentuckyCounty, 
		polyFine = kentuckyTract, 
    cellsCoarse = 5, 
		cellsFine = 100,
    bw = c(20,50, 100) * 1000, 
    ncores = ncores, 
    idFile = 'id.grd', 
		offsetFile = 'offset.grd'
)

sum(kentuckyTract$expected)
sum(values(lemRaster$offset$offset), na.rm=TRUE)*prod(res(lemRaster$offset))

if(!interactive()) {
	pdf("krasters.pdf")
}

oCol = colourScale(
		lemRaster$offset$offset,
		breaks=10, style='equal',
		dec=7, transform='sqrt'
		)
map.new(kentuckyTract[35,], buffer=10000)
plot(lemRaster$offset$offset,
		col=oCol$col, breaks=oCol$breaks,
		legend=FALSE, add=TRUE)
plot(kentuckyTract,add=TRUE)
plot(kMap, add=TRUE)

legendBreaks("topleft", 
		oCol$breaks*10^6, col=oCol$col,
		title='cases/km2', bg='white')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)

if(!interactive()) {
	dev.off()
}

lemSmoothMat = smoothingMatrix(
		rasterObjects = lemRaster, 
    ncores = ncores) 

dim(lemSmoothMat$smoothingArray)
length(lemSmoothMat$polyCoarse)
length(kentuckyTract)
length(unique(values(lemRaster$rasterFine$idFine)))
length(unique(values(lemRaster$offset$offset)))
