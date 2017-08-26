## ----setup, message = FALSE----------------------------------------------

# specify number of grid cells and number of cores for computations in parallel 
ncores = 2
cellsFine = 80
cellsCoarse = 8
nsim = 4
nxv = 4
Sbw = seq(10, 35, by = 5) * 1000
threshold = c(1, 1.25, 1.5, 2)
Nboot = 100 

# on unix systems results will be saved in /tmp/localEMusername
path = file.path(dirname(tempdir()),  paste('localEM', Sys.info()['user'], sep = ''))

# specify number of grid cells for simulation
cellsSimulate = 200
if(!requireNamespace('RandomFields', quietly = TRUE)) {
    cellsSimulate = 100
}

dir.create(path, showWarnings = FALSE, recursive = TRUE)
dir.create('cache', showWarnings = FALSE, recursive = TRUE)


## ----data----------------------------------------------------------------

library('mapmisc')
library('rgdal')

data('kentuckyCounty', package = 'localEM') 
data('kentuckyTract', package = 'localEM') 
data('kMap', package = 'localEM')


## ----simcases, message = FALSE, warning = FALSE, cache = (basename(getwd()) == 'www'), cache.path = 'cache'----

kentuckyOffset = geostatsp::spdfToBrick(
    kentuckyTract,
    geostatsp::squareRaster(kentuckyTract, cellsSimulate),
    pattern = '^expected$',
    logSumExpected = TRUE)

set.seed(0)
kCases = geostatsp::simLgcp(
    param = c(mean = 0, variance = 0.4^2, range = 120 * 1000, shape = 2),
    covariates = list(logExpected = kentuckyOffset), 
    offset = 'logExpected', n = nsim)


## ----aggregateCases------------------------------------------------------

kCases$agg = lapply(
    kCases[grep('^events[[:digit:]]+?', names(kCases))],
    function(qq) over(qq, kentuckyCounty)[,'id'])
countyCounts = as.data.frame(lapply(
    kCases$agg,  
    function(xx) as.vector(table(xx, exclude = NULL)[as.character(kentuckyCounty$id)])))
countyCounts[is.na(countyCounts)] = 0
names(countyCounts) = gsub('^events', 'count', names(countyCounts))
rownames(countyCounts) = as.character(kentuckyCounty$id)
kentuckyCounty = merge(kentuckyCounty, countyCounts, by.x = 'id', by.y = 'row.names')


## ----plotOffset, echo = FALSE, warning = FALSE, fig.cap = 'Events for Simulation 1', fig.subcap = c('Offset', 'Relative Intensity', 'Events', 'Counts'), fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

# offset map
oCol = colourScale(
    c(0, exp(maxValue(kentuckyOffset))), 
    breaks = 10, style = 'equal', dec = 7, 
    transform = 'sqrt')

map.new(kentuckyTract, 
    mar = c(0, 0, 1, 0), 
	  main = 'Offsets')
plot(kentuckyOffset, 
    col = oCol$col, 
    breaks = pmax(-100,log(oCol$breaks)), 
    legend = FALSE, 
    add = TRUE)
plot(kMap, add = TRUE)
legendBreaks('topleft', 
    breaks = oCol$breaks * 10^6, col = oCol$col, 
    title = expression(plain('cases/km')^2), 
    bg = 'white')


# simulated risk map
iCol = colourScale(
    kCases$raster$relativeIntensity1, 
    breaks = 8, style = 'equal', dec = -log10(0.5))

map.new(kentuckyTract, 
    mar = c(0, 0, 1, 0), 
	  main = 'Simulated Risk')
plot(kCases$raster$relativeIntensity1, 
    col = iCol$col, breaks = iCol$breaks, 
    legend = FALSE, 
    add = TRUE)
plot(kMap, add = TRUE)
legendBreaks('topleft', 
    col = iCol$col, breaks = iCol$breaks, 
    title = 'RR', 
    bg = 'white')


# simulated events map
map.new(kentuckyTract, 
    mar = c(0, 0, 1, 0), 
	  main = 'Simulated Events')
plot(kMap, add = TRUE)
points(kCases$events1, col = '#FF000030', pch = 20, cex = 0.5)
scaleBar(kentuckyCounty, 
    pos = 'topleft', 
    bg = 'white')


# simulated counts
cCol = colourScale(
    kentuckyCounty$count1, 
    breaks = 10, style = 'quantile', dec = -1)

map.new(kentuckyTract, 
    mar = c(0, 0, 1, 0), 
    main = 'Aggregated Counts')
plot(kentuckyCounty, col = cCol$plot, add = TRUE)
plot(kMap, add = TRUE)
legendBreaks('topleft', cCol)
scaleBar(kentuckyCounty, 
    pos = 'topleft', 
    inset = c(0.2, 0), 
    bg = 'white')


## ----cvOne, warning = FALSE----------------------------------------------

library('localEM')

fileHere = file.path(path, 'xvKentucky.RData')

if(!file.exists(fileHere)) {

    xvKentucky = lemXv(
        cases = kentuckyCounty[,c('id','count1')], 
        population = kentuckyTract, 
        cellsCoarse = cellsCoarse,  
        cellsFine = cellsFine, 
        bw = Sbw, 
        xv = nxv, 
        ncores = ncores, 
        path = path, 
        verbose = TRUE)

    save(xvKentucky, file = fileHere)
} else {

    load(fileHere)
}


## ----cvPlotOne, echo = FALSE, fig.cap = 'Cross-validation Scores for Simulation 1', fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

plot(xvKentucky$xv[,1]/1000, xvKentucky$xv[,2], 
    main = 'CV Scores for Simulation 1', 
    xlab = 'km', ylab = '-log p', 
    type = 'o')


## ----cvAll, warning = FALSE, cache = (basename(getwd()) == 'www'), cache.path = 'cache'----

xvAllKentucky = lemXv(
    cases = kentuckyCounty@data[,grep('^count[[:digit:]]', names(kentuckyCounty))], 
    lemObjects = xvKentucky$smoothingMatrix, 
    ncores = ncores, 
    path = 'cache', 
    verbose = TRUE)


## ----cvPlotAll, echo = FALSE, fig.cap = 'Cross-validation Scores for All Simulations', fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

maxY = min(max(xvAllKentucky$xv[,-1]), 30)

matplot(xvAllKentucky$xv[,'bw']/1000, xvAllKentucky$xv[,-1], 
    main = 'CV Scores for All Simulations', 
    xlab = 'km', ylab = '-log p', 
    type = 'l', lty = 1, col = 1:4, 
    ylim = c(0, maxY))
legend('topright', paste('Simulation ', 1:4, sep = ''), 
    lty = 1, col = 1:4, 
    cex = 0.6)


## ----plotRiskAll, echo = FALSE, warning = FALSE, fig.cap = 'Risk Estimation with Optimal Bandwidths', fig.subcap = paste('Simulation ', 1:4, sep = ''), fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

toPlot = xvAllKentucky$riskEst

rCol = colourScale(
    toPlot, 
    breaks = 8, style = 'equal', dec = 2)

bw = as.numeric(gsub('^bw', '', xvAllKentucky$bw))

for(inR in 1:nlayers(toPlot)) {

    map.new(kentuckyCounty, 
        mar = c(0, 0, 1, 0), 
        main = paste('Simulation ', inR, ': Estimated Risk, bw=', bw[inR] / 1000, ' km', sep = ''))
    plot(toPlot[[inR]], 
        col = rCol$col, breaks = rCol$breaks, 
        legend = FALSE, 
        add = TRUE)
    plot(kMap, add = TRUE)
    legendBreaks('topleft', rCol, bg = 'white')
}


## ----riskOne, warning = FALSE--------------------------------------------

riskKentucky = riskEst(
    cases = kentuckyCounty[,c('id','count1')], 
    lemObjects = xvKentucky$smoothingMatrix, 
    bw = Sbw, 
    ncores = ncores, 
    path = path, 
    verbose = TRUE)


## ----plotRiskOne, echo = FALSE, warning = FALSE, fig.cap = 'Risk Estimation for Simulation 1', fig.subcap = paste('Risk, bw=', Sbw / 1000, ' km', sep = ''), fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

# estimated risk maps
toPlot = riskKentucky$riskEst

rCol = colourScale(
    toPlot, 
    breaks = 8, style = 'equal', dec = 2)

for(inR in 1:nlayers(toPlot)) {

    map.new(kentuckyCounty, 
        mar = c(0, 0, 1, 0), 
    	  main = paste('Simulation 1: Estimated Risk, bw=', Sbw[inR] / 1000, ' km', sep = ''))
    plot(toPlot[[inR]], 
        col = rCol$col, breaks = rCol$breaks, 
        legend = FALSE, 
        add = TRUE)
    plot(kMap, add = TRUE)
    legendBreaks('topleft', rCol, bg = 'white')
}


## ----excProbOne, warning = FALSE-----------------------------------------

excProbKentucky = excProb(
    lemObjects = xvKentucky, 
    threshold = threshold, 
    Nboot = Nboot, 
    ncores = ncores, 
    path = path, 
    verbose = TRUE)


## ----plotExcProbOne, echo = FALSE, warning = FALSE, fig.cap = 'Exceedance Probabilities for Simulation 1', fig.subcap = paste('Threshold ', threshold, sep = ''), fig.height = 4, fig.width = 6, out.width = '80%', fig.ncol = NULL----

toPlot = excProbKentucky$excProb

pCol = colourScale(
    toPlot, 
    breaks = c(0, 0.2, 0.8, 0.95, 1), style = 'fixed', 
    col = c('green', 'yellow', 'orange', 'red'))								

for(inT in 1:nlayers(toPlot)) {

    map.new(kentuckyCounty, 
        mar = c(0, 0, 1, 0), 
    	  main = paste('Simulation 1: Exceedance Probabilities, t=', threshold[inT], sep = ''))
    plot(toPlot[[inT]], 
        col = pCol$col, breaks = pCol$breaks, 
        legend = FALSE, 
        add = TRUE)
    plot(kMap, add = TRUE)
}
legendBreaks('topleft', pCol, bg = 'white')


