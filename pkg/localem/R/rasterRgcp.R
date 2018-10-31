
rasterPartitionRgcp = function(
	coarsePolyList,
	finePolyList,
	Ncoarse, Nfine,
	pathBase, verbose=FALSE,
	mc.cores=1
	) {

	# add an X before integers in names
	# because raster will do it if we don't
	names(coarsePolyList) = 
		gsub("(^[[:digit:]])",
			"X\\1", names(coarsePolyList))

	names(finePolyList) = 
		gsub("(^[[:digit:]])",
			"X\\1", names(finePolyList))


	Dmap = names(coarsePolyList)[1]

	if(verbose) {
		cat('first rasterization ')
	}
	rasterList1 = 
	localEM::rasterPartition(
		polyCoarse = coarsePolyList[[Dmap]],
		polyFine = finePolyList[[Dmap]],
		Ncoarse,
		Nfine,
		NULL, NULL, NULL,
		path=file.path(pathBase, Dmap),
		verbose=verbose
		)
	fineCells = deratify(rasterList1$rasterFine, 'cellCoarse')
	coarseCells = rasterList1$rasterCoarse

	Smap = setdiff(names(coarsePolyList)[-1], Dmap)
	MoreArgs = list(
		cellsCoarse = coarseCells[['cellCoarse']],
		cellsFine = fineCells,
		bw=NULL, focalSize=NULL, xv=NULL,
		verbose=verbose
	)

	if(verbose) {
		cat('remaining rasterizations ')
	}

	if(mc.cores <= 1) {
		cl = NULL
		rasterList2 = Map(localEM::rasterPartition,
			polyCoarse = coarsePolyList[Smap],
			polyFine = finePolyList[Smap],
			path=file.path(pathBase, Dmap),
			MoreArgs = MoreArgs
			)
	} else {
		cl = parallel::makeCluster(mc.cores, type='PSOCK', methods=TRUE)
		parallel::clusterEvalQ(cl, require('raster'))
		parallel::clusterEvalQ(cl, require('Matrix'))

		rasterList2 = parallel::clusterMap(
			cl = cl,
			fun=localEM::rasterPartition,
			polyCoarse = coarsePolyList[Smap],
			polyFine = finePolyList[Smap],
			path=file.path(pathBase, Smap),
			MoreArgs = MoreArgs
			)
	}

	if(FALSE) { # for debugging
		rasterList2 = list()
		for(DDmap in Smap)
			rasterList2[[DDmap]] = localEM::rasterPartition(
			polyCoarse = coarsePolyList[[DDmap]],
			polyFine = finePolyList[[DDmap]],
			path=file.path(pathBase, DDmap),
		cellsCoarse = coarseCells[['cellCoarse']],
		cellsFine = fineCells,
		bw=NULL, focalSize=NULL, xv=NULL,
		verbose=verbose
				)
	}

	rasterList = c(
		list(x=rasterList1),
		rasterList2
		)

	toBrick = lapply(rasterList, function(xx) xx$offset[['offset']])
	offsetStack = do.call(raster::stack, toBrick)
	offsetStack = offsetStack * prod(res(offsetStack))
	names(offsetStack) = names(rasterList) = c(Dmap, Smap)

# sum of offsets in intersection of
# coarse polygons and coarse cells
	if(verbose) {
		cat(' intersections ')
	}

	if(!is.null(cl)) {
		Oijl = parallel::clusterMap(
			cl = cl,
			fun=rgcpOijlFun,
			rasterHere = rasterList,
			offsetStackHere = as.list(offsetStack),
			MoreArgs = list(coarseCells = coarseCells)
			)
		parallel::stopCluster(cl)
	} else {
		Oijl = Map(
			fun=rgcpOijlFun,
			rasterHere = rasterList,
			offsetStackHere = as.list(offsetStack),
			MoreArgs = list(coarseCells = coarseCells)
			)
	}

	offsetIL = lapply(Oijl, function(xx) apply(xx,MARGIN=1, FUN=sum))
	offsetL = apply(do.call(cbind, offsetIL), 1, sum)
	offsetMat = Matrix::Diagonal(length(offsetL), offsetL)

	if(verbose) {
		cat('\n')
	}
	
	list(
		Oijl = Oijl,
		cells =  list(coarse=coarseCells, fine=fineCells),
		offset = list(raster = offsetStack, matrix = offsetMat)
		)
}


rgcpOijlFun = function(rasterHere, offsetStackHere, coarseCells) {
	stuff = raster::deratify(rasterHere$rasterFine, 
		c('idCoarse','cellCoarse'))
	notNa= which(!is.na(raster::values(offsetStackHere)))
	offsetSumHere = cbind(
		offset=raster::values(offsetStackHere)[notNa], 
		as.data.frame(raster::deratify(rasterHere$rasterFine, 
			c('idCoarse','cellCoarse')))[notNa,])
	offsetAggHere = aggregate(
		offsetSumHere[,'offset'],
		offsetSumHere[,c('idCoarse','cellCoarse')],
		sum, na.rm=TRUE)
	offsetAggHere=offsetAggHere[offsetAggHere$x > 0, ]
	Matrix::sparseMatrix(
		offsetAggHere[,'cellCoarse'],
		offsetAggHere[,'idCoarse'],
		x = offsetAggHere[,'x'],
		dim = c(
			raster::ncell(coarseCells),
			length(rasterHere$polyCoarse)),
		dimnames = list(
			raster::values(coarseCells),
			rasterHere$polyCoarse$idCoarse
			)
		)
}