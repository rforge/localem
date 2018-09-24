
rasterPartitionRgcp = function(
	coarsePolyList,
	finePolyList,
	Ncoarse, Nfine,
	pathBase
	) {

	rasterList = list()

	DmapFirst = Dmap = names(coarsePolyList)[1]

	rasterList[[ Dmap ]] = 
	localEM::rasterPartition(
		polyCoarse = coarsePolyList[[Dmap]],
		polyFine = finePolyList[[Dmap]],
		Ncoarse,
		Nfine,
		NULL, NULL, NULL,
		path=file.path(pathBase, Dmap),
		verbose=TRUE
		)
	fineCells = deratify(rasterList[[DmapFirst]]$rasterFine, 'cellCoarse')
	coarseCells = rasterList[[DmapFirst]]$rasterCoarse

	for(Dmap in names(coarsePolyList)[-1]) {
		rasterList[[ Dmap ]] = 
		localEM::rasterPartition(
			polyCoarse = coarsePolyList[[Dmap]],
			polyFine = finePolyList[[Dmap]],
			cellsCoarse = coarseCells[['cellCoarse']],
			cellsFine = fineCells,
			NULL, NULL, NULL,
			path=file.path(pathBase, Dmap),
			verbose=TRUE)
	}

	toBrick = lapply(rasterList, function(xx) xx$offset[['offset']])
	names(toBrick)[1] = 'x'
	offsetStack = do.call(raster::stack, toBrick)
	offsetStack = offsetStack * prod(res(offsetStack))
	names(offsetStack) = names(rasterList)

# sum of offsets in intersection of
# coarse polygons and coarse cells
	Oijl = list()
	for(Dmap in names(rasterList)) {
		stuff = deratify(rasterList[[Dmap]]$rasterFine, 
			c('idCoarse','cellCoarse'))
		notNa= which(!is.na(values(offsetStack[[Dmap]])))
		offsetSumHere = cbind(
			offset=values(offsetStack[[Dmap]])[notNa], 
			as.data.frame(deratify(rasterList[[Dmap]]$rasterFine, 
				c('idCoarse','cellCoarse')))[notNa,])
		offsetAggHere = aggregate(
			offsetSumHere[,'offset'],
			offsetSumHere[,c('idCoarse','cellCoarse')],
			sum, na.rm=TRUE)
		offsetAggHere=offsetAggHere[offsetAggHere$x > 0, ]

		Oijl[[Dmap]] = sparseMatrix(
			offsetAggHere[,'cellCoarse'],
			offsetAggHere[,'idCoarse'],
			x = offsetAggHere[,'x'],
			dim = c(
				ncell(coarseCells),
				length(rasterList[[Dmap]]$polyCoarse)),
			dimnames = list(
				values(coarseCells),
				rasterList[[Dmap]]$polyCoarse$idCoarse
				)
			)
	}

	offsetIL = lapply(Oijl, function(xx) apply(xx,MARGIN=1, FUN=sum))
	offsetL = apply(do.call(cbind, offsetIL), 1, sum)
	offsetMat = Matrix::Diagonal(length(offsetL), offsetL)

	list(
		Oijl = Oijl,
		cells =  list(coarse=coarseCells, fine=fineCells),
		offset = list(raster = offsetStack, matrix = offsetMat)
		)
}
