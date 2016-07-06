#' # Local-EM example with Kentucky

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

#' # Rasters

lemRaster = rasterPartition(
		polyCoarse = kentuckyCounty, 
		polyFine = kentuckyTract, 
    cellsCoarse = 5, 
		cellsFine = 120,
    bw = c(5, 10,15, 20, 30) * 1000, 
    ncores = ncores
)

sum(kentuckyTract$expected)
sum(values(lemRaster$offset$offset), na.rm=TRUE)*prod(res(lemRaster$offset))

if(!interactive()) {
	pdf("kentucky.pdf")
}

oCol = colourScale(
		lemRaster$offset$offset,
		breaks=10, style='equal',
		dec=7, transform='sqrt'
)
map.new(kentuckyTract)
plot(lemRaster$offset$offset,
		col=oCol$col, breaks=oCol$breaks,
		legend=FALSE, add=TRUE)
plot(kMap, add=TRUE)

legendBreaks("topleft", 
		oCol$breaks*10^6, col=oCol$col,
		title='cases/km2', bg='white')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)

#' # Simulate cases

kLogOffset = log(lemRaster$offset$offset)
names(kLogOffset) = 'logOffset'

set.seed(0)
kCases = geostatsp::simLgcp(
		param=c(mean=-1, variance=0.4^2, 
        range=120*1000, shape=2),
		covariates = kLogOffset,
		offset='logOffset')

length(kCases$events)
sum(values(kCases$raster$intensity), na.rm=TRUE)*prod(res(kCases$raster))


iCol = colourScale(
		kCases$raster$relativeIntensity,
		breaks=8, dec=1, style='equal'
)
map.new(kentuckyTract)
plot(kCases$raster$relativeIntensity,
		col=iCol$col, breaks=iCol$breaks,
		legend=FALSE, add=TRUE)
plot(kMap, add=TRUE)

	legendBreaks("topleft", 
			iCol, title='relative risk',
			bg='white')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)
		
		
map.new(kentuckyTract)
plot(kMap, add=TRUE)
points(kCases$events, col='#FF000030')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)

countyCounts = table(over(kCases$events, kentuckyCounty)$id)
kentuckyCounty$count = countyCounts[kentuckyCounty$id]
kentuckyCounty$count[is.na(kentuckyCounty$count)]= 0

#' # Smoothing Matrix

lemSmoothMat = smoothingMatrix(
		rasterObjects = lemRaster, 
    ncores = ncores,
		verbose=TRUE) 

dim(lemSmoothMat$smoothingArray)
dimnames(lemSmoothMat$smoothingArray)[[3]]


lemCv = lemXv(x = kentuckyCounty, 
    lemObjects = lemSmoothMat, 
    ncores = ncores) 

par(mar=c(5,5,1,1))
plot(do.call(cbind, lemCv), xlab='bw', ylab='cv', log='x')

#' # Risk estimates

lemRisk = riskEst(x = kentuckyCounty,
    lemObjects = lemSmoothMat,
    bw = lemCv$bw[which.min(lemCv$cv)]
)


rCol = colourScale(
		lemRisk,
		breaks=10, style='equal',
		dec=1
)

map.new(kentuckyTract)

plot(lemRisk,
		col=rCol$col, breaks=rCol$breaks,
		legend=FALSE, add=TRUE)

plot(kMap, add=TRUE)

legendBreaks("topleft", 
		rCol,
		title='rr', bg='white')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)


#' # Exceedance probabilities

#lemExcProb = excProb(x = kentuckyCounty, 
#                    lemObjects = lemRisk, 
#                    threshold = c(1, 1.1), 
#                    Nboot = 10, 
#                    ncores = ncores) 

if(!interactive()) {
	dev.off()
}
