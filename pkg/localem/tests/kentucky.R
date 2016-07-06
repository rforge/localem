#' # Local-EM example with Kentucky

#+ setup

ncores = 1+ (.Platform$OS.type=='unix')

library('localEM')
library('mapmisc')

data('kentuckyCounty') 
data('kentuckyTract')

data('kMap')

kentuckyTract$expected = kentuckyTract$expected/5
#'

#+ forFigures, eval=FALSE, include=FALSE
if(!interactive()) {
	pdf("kentucky.pdf")
}
#'





#' # Rasters

#+ rasters, cache=TRUE
lemRaster = rasterPartition(
		polyCoarse = kentuckyCounty, 
		polyFine = kentuckyTract, 
    cellsCoarse = 5, 
		cellsFine = 150,
    bw = c(8, 10,15, 20, 30) * 1000, 
    ncores = ncores
)

sum(kentuckyTract$expected)
sum(values(lemRaster$offset$offset), na.rm=TRUE)*prod(res(lemRaster$offset))
#'


#+ plotOffset, echo=FALSE
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
#'

#' 
#' # Simulate cases
#' 

#+ simcases, cache=TRUE
kLogOffset = log(lemRaster$offset$offset)
names(kLogOffset) = 'logOffset'

set.seed(0)
kCases = geostatsp::simLgcp(
		param=c(mean=0, variance=0.4^2, 
        range=120*1000, shape=2),
		covariates = kLogOffset,
		offset='logOffset')

length(kCases$events)
sum(values(kCases$raster$intensity), na.rm=TRUE)*prod(res(kCases$raster))


countyCounts = table(over(kCases$events, kentuckyCounty)$id)
kentuckyCounty$count = countyCounts[kentuckyCounty$id]
kentuckyCounty$count[is.na(kentuckyCounty$count)]= 0

#'


#+ plotCases, fig.cap = 'simulated', fig.subcap = c('intensity','events'), echo=FALSE
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
#+		


#' # Smoothing Matrix

#+ smoothingMatrix, cache=TRUE
lemSmoothMat = smoothingMatrix(
		rasterObjects = lemRaster, 
    ncores = ncores,
		verbose=TRUE) 

dim(lemSmoothMat$smoothingArray)
dimnames(lemSmoothMat$smoothingArray)[[3]]
#'

#+ cv, cache=TRUE
lemCv = lemXv(x = kentuckyCounty, 
    lemObjects = lemSmoothMat, 
    ncores = ncores) 
bestBw = lemCv$bw[which.min(lemCv$cv)]
#'

#+ plotCv, echo=FALSE
par(mar=c(5,5,1,1))
plot(do.call(cbind, lemCv), xlab='bw', ylab='cv', log='x', pch=16, col='red')
abline(v=bestBw, lty=3)
#'

#' # Risk estimates

#+ riskExt, cache=TRUE
lemRisk = riskEst(x = kentuckyCounty,
    lemObjects = lemSmoothMat,
    bw = bestBw
)
#'

#+ plotRisk, echo=FALSE
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
#'

#' # Exceedance probabilities

#+ excProb, cache=TRUE
lemExcProb = excProb(x = kentuckyCounty, 
    lemObjects = lemSmoothMat, 
		estimate=lemRisk,
		bw=bestBw,
    threshold = c(1.2, 1.5), 
    Nboot = 10, 
    ncores = ncores) 
#'

#+ plotExcProb, echo=FALSE
pCol = colourScale(
		lemExcProb$threshold.1.5,
		breaks=5, style='equal'
)								

map.new(kentuckyTract)

plot(lemExcProb$threshold.1.5,
		col=pCol$col, breaks=pCol$breaks,
		legend=FALSE, add=TRUE)

plot(kMap, add=TRUE)

legendBreaks("topleft", 
		pCol,
		title='prob', bg='white')

scaleBar(kentuckyCounty,
		'topleft',
		inset=c(0.3, 0.1),
		bg='white'
)
#'


#' # ROC


#+ lemRoc
lemRoc = geostatsp::spatialRoc(
 fit= lemExcProb, 
 truth = kCases, 
		border=NULL, random=FALSE,
		prob = NULL, 
 spec = seq(0,1,by=0.01)
)
#'


#+ plotRoc, fig.cap='ROC', echo=FALSE

matplot(lemRoc[,1], lemRoc[,c(2,3)], ylab='sens', xlab='1-spec', type='l',
		col=2:3)
legend("bottomright", lty=1, col=2:3, legend=colnames(lemRoc)[-1], 
		title='threshold')

#'

#+ forFiguresEnd, eval=FALSE, include=FALSE
if(!interactive()) {
	dev.off()
}

#'


