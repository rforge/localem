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
		cellsFine = 100,
    bw = c(5, 10, 20, 40) * 1000, 
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

do.call(rbind, lemCv)

#' # Risk estimates

lemRisk = riskEst(x = kentuckyCounty,
    lemObjects = lemSmoothMat,
    bw = lemCv$bw[2]
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
