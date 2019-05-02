
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

  if(is.numeric(Ncoarse)) {
	elist = lapply(coarsePolyList, extent)
	if(length(elist) > 1) {
		names(elist)[1:2] = c('x','y')
		coarseExtent = do.call(raster::merge, elist)
	} else {
		coarseExtent = elist[[1]]
	}

    rasterCoarse=geostatsp::squareRaster(coarseExtent, Ncoarse)
    raster::projection(rasterCoarse) = 
    	raster::projection(coarsePolyList[[1]])
    values(rasterCoarse) = seq(1, ncell(rasterCoarse))
    names(rasterCoarse) = "cellCoarse"
  } else {
    rasterCoarse = Ncoarse
  }


  if(is.numeric(Nfine)) {
    rasterFine = disaggregate(rasterCoarse,
      ceiling(Nfine/ncol(rasterCoarse)))
    names(rasterFine) = 'cellCoarse'
  } else {
    rasterFine = Nfine
  }


	Smap = names(coarsePolyList)
	MoreArgs = list(
		cellsCoarse = rasterCoarse[['cellCoarse']],
		cellsFine = rasterFine,
		bw=NULL, focalSize=NULL, xv=NULL,
		verbose=verbose
	)

	if(verbose) {
		cat('rasterizations ')
	}

	if(FALSE) {
		# for debugging
			polyCoarse = coarsePolyList[[Smap[3]]]
			polyFine = finePolyList[[Smap[3]]]
			path=file.path(pathBase, Smap[3])

			rasterList = list()
			for(Dmap in Smap) {
				cat(Dmap)
				rasterList[[Dmap]] = localEM::rasterPartition(
					polyCoarse = coarsePolyList[[Dmap]],
					polyFine = finePolyList[[Dmap]],
					path = file.path(pathBase, Dmap),
					cellsCoarse = rasterCoarse[['cellCoarse']],
					cellsFine = rasterFine,
					bw=NULL, focalSize=NULL, xv=NULL,
					verbose=verbose
					)
			}
	}

#	dir.create(file.path(pathBase, Smap), showWarnings=FALSE)
	if(mc.cores <= 1) {
		cl = NULL
		rasterList = Map(localEM::rasterPartition,
			polyCoarse = coarsePolyList,
			polyFine = finePolyList,
			path=file.path(pathBase, Smap),
			MoreArgs = MoreArgs
			)
	} else {
		cl = parallel::makeCluster(mc.cores, type='PSOCK', methods=TRUE)
		parallel::clusterEvalQ(cl, require('raster'))
		parallel::clusterEvalQ(cl, require('Matrix'))

		rasterList = parallel::clusterMap(
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

	rasterList = c(
		list(x=rasterList1),
		rasterList2
		)
	}

	toBrick = lapply(rasterList, function(xx) xx$offset[['offset']])
	names(toBrick)[1] = 'x'
	offsetStack = do.call(raster::stack, toBrick)
	offsetStack = offsetStack * prod(res(offsetStack))
	names(offsetStack) = names(rasterList) = names(coarsePolyList)

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
			MoreArgs = list(coarseCells = rasterCoarse[['cellCoarse']])
			)
		parallel::stopCluster(cl)
	} else {
		Oijl = Map(
			f=rgcpOijlFun,
			rasterHere = rasterList,
			offsetStackHere = as.list(offsetStack),
			MoreArgs = list(coarseCells = rasterCoarse[['cellCoarse']])
			)
	}

	offsetIL = lapply(Oijl, function(xx) apply(xx, MARGIN=1, FUN=sum))
	offsetL = apply(do.call(cbind, offsetIL), 1, sum)
	offsetMat = Matrix::Diagonal(length(offsetL), offsetL)

	if(verbose) {
		cat('\n')
	}
	
	list(
		Oijl = Oijl,
		cells =  list(coarse=rasterCoarse[['cellCoarse']], fine=rasterFine),
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