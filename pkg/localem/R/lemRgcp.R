
emsRgcp = function(
	range,
	sd,
	data,
	shape,
	mc.cores = 1,
	maxIter = 200,
	tol = .Machine$double.eps^(1/3),
	reduce = TRUE,
	verbose=FALSE
	) {

	# do the last range first, because it takes the longest
	Srange = sort(unique(range), decreasing=TRUE)
	Ssd = sort(unique(sd))

	if(missing(shape)) {
		shape = data$shape
	} 

	nCoresOuter = mc.cores[1]
	if(length(mc.cores)<=1) {
		nCoresOuter = pmin(mc.cores, length(Srange))
		nCoresInner = max(c(1,round(mc.cores/nCoresOuter)))
	} else {
		nCoresInner = mc.cores[2]
	}
	if(verbose){
		cat("outer cores ", nCoresOuter, ', inner cores ', nCoresInner, '\n')
		cat('ranges ', length(Srange), ' variances ', length(Ssd), 
			' intercepts ', length(unique(data$obsIntercept$intercept)), 
			' datasets ',length(unique(data$obsIntercept$obs)),'\n')
	}




	if(nCoresOuter > 1) {
		outerCluster = parallel::makeCluster(nCoresOuter, type='PSOCK', methods=TRUE)
		parallel::clusterEvalQ(outerCluster, require('Matrix'))
		parallel::clusterEvalQ(outerCluster, require('data.table'))
		parallel::clusterEvalQ(outerCluster, require('localEM'))
		parallel::clusterExport(outerCluster, 
			varlist = 'nCoresInner',
			envir = environment())
		parallel::clusterEvalQ(outerCluster, 
			data.table::setDTthreads(pmax(1,nCoresInner)))

		if(verbose){
			cat('cluster pid: ', 
				unlist(parallel::clusterCall(outerCluster, function() Sys.getpid())), 
				'\n')
		} 
	} else {
		outerCluster = NULL
		require('data.table')
	}

	if(nCoresInner > 1) {
		if(nCoresOuter <= 1) {
			innerCluster <- parallel::makeCluster(nCoresInner, type='PSOCK', methods=TRUE)
			parallel::clusterEvalQ(innerCluster, require('Matrix'))
			parallel::clusterEvalQ(innerCluster, require('data.table'))
			parallel::clusterEvalQ(innerCluster, require('localEM'))
			parallel::clusterExport(innerCluster, 
				varlist = 'nCoresInner',
				envir = environment())
			data.table::setDTthreads(1)
			if(verbose){
				cat('inner pid ',
					unlist(parallel::clusterCall(innerCluster, Sys.getpid)),
					'\n')
			}
		} else {
			parallel::clusterEvalQ(outerCluster, 
				{innerCluster <- parallel::makeCluster(nCoresInner, type='PSOCK', methods=TRUE);
				parallel::clusterEvalQ(innerCluster, require('Matrix'))
				parallel::clusterEvalQ(innerCluster, require('data.table'))
				parallel::clusterEvalQ(innerCluster, require('localEM'))
				parallel::clusterExport(innerCluster, 
					varlist = 'nCoresInner',
					envir = environment())
				data.table::setDTthreads(1)
			})

			if(verbose){
				cat('inner pid\n')
				print(simplify2array(
					parallel::clusterEvalQ(outerCluster, 
					{
						parallel::clusterCall(innerCluster, Sys.getpid)
					})
					))
				cat('\n')
			} 
		}

	} else {
		innerCluster = NULL
		parallel::clusterEvalQ(outerCluster, {
			innerCluster = NULL
		})
	}



	if(!is.null(outerCluster)) {
		res = parallel::clusterMap(outerCluster, 
			emsOneRange,
			range = Srange,
			MoreArgs = list(
				Ssd = Ssd,
				data=data,
				maxIter = maxIter,
				tol = tol,
				cl = 'innerCluster',
				reduce = reduce,
				verbose = verbose + verbose*(nCoresOuter==1)
				)) 
	} else {

		res = Map(
			emsOneRange,
			range = Srange,
			MoreArgs = list(
				Ssd = Ssd,
				data=data,
				maxIter = maxIter,
				tol = tol,
				cl = innerCluster,
				reduce = reduce,
				verbose = 2*verbose
				)) 
	}

	if(nCoresOuter > 1) {
		if(nCoresInner > 1) {
			parallel::clusterEvalQ(outerCluster, parallel::stopCluster(innerCluster))
		}
		parallel::stopCluster(outerCluster)
	} else {
		if(!is.null(innerCluster)) {
			parallel::stopCluster(innerCluster)	
		}
	}

	logLik = do.call(rbind, lapply(res, function(qq) qq$logL))
	rownames(logLik) = 1:nrow(logLik)

	coarseCells=data$coarseCells

	if(reduce) {
		onlyMax = merge(
			logLik,
			aggregate(logLik[,'logL',drop=FALSE], logLik[,'obs', drop=FALSE], max),
			all=FALSE)
		onlyMax$Drange = match(onlyMax$range, Srange)

		thetaSub = list()
		for(Dobs in 1:nrow(onlyMax)) {
			thetaSub[[as.character(onlyMax[Dobs,'obs'])]] = res[[
			onlyMax[Dobs,'Drange']
			]]$theta[ , as.character(onlyMax[Dobs,'obs'])]
		}
		thetaSub = do.call(cbind, thetaSub)

		thetaArray = array(thetaSub, 
			c(ncol(coarseCells), nrow(coarseCells), 
				ncol(thetaSub)))

		thetaBrick = raster::brick(thetaArray, 
			xmin(coarseCells), xmax(coarseCells),
			ymin(coarseCells), ymax(coarseCells), 
			projection(coarseCells),
			transpose=TRUE)
		names(thetaBrick) = colnames(thetaSub)
	} else {
		thetaArray = lapply(res, function(xx) xx$theta)
		thetaArray = do.call(abind::abind, c(thetaArray, list(along=4)))
#	thetaArray = array(thetaArray, 
#		      c(ncol(coarseCells), nrow(coarseCells), dim(thetaArray)[-1])
#		      )
		thetaBrick = thetaArray
		attributes(thetaBrick)$raster = coarseCells
	}

	logLik$Nlogsd = nrow(data$offsetMatrix)*log(logLik$sd)
#logLik$logL2 = logLik$yLogSum - logLik$sumDenom - logLik$halfLogDet -
#	 logLik$thetaCrossprod/(2*logLik$sd^2) + 
#	 logLik$invCorMatHalfLogDet - logLik$Nlogsd 


	logLprof = merge(
		logLik,
		aggregate(
			logLik[,'logL',drop=FALSE], 
			logLik[,c('sd','range','obs'), drop=FALSE], max),
		all=FALSE)
	Sidvars = c('obs','intercept','sd','range')
	logLikMelt = reshape2::melt(logLik, 
		id.vars = Sidvars,
		measure.vars = setdiff(colnames(logLik), Sidvars))

	logLikArray = reshape2::acast(logLikMelt, obs ~ intercept ~ sd ~ range ~ variable)

	if(verbose){
		cat('\n')
	}

	forMLE = as.data.table(logLik)
	mle = forMLE[,
		list(range = range[which.max(logL)], 
			sd=sd[which.max(logL)], 
			intercept = intercept[which.max(logL)], 
			logL = max(logL)), 
		by=obs]

	res = list(mle = as.data.frame(mle), 
		theta = thetaBrick, logL = logLik, profL = logLprof, 
		array = logLikArray, data=data)
	res
}


emsRgcpPred = function(
	x, 
	param = x$mle,
	theta = x$theta,
	data = x$data,
	cl = NULL
	) {

	if(class(theta) == 'array') {
		warning("haven't written this part yet")
	}

	if(is.null(param$shape)) {
		param$shape = data$shape
	}
	param$variance = param$sd^2

	for2deriv = objectsForLikelihood(
		data$Oijl, 
		data$y, 
		as.matrix(theta)^2) 

	if(is.character(cl)) {
		cl = get(cl)
	}


	if(!is.null(cl)) {
		twoDeriv = parallel::clusterMap(cl, 
			derivDiag,
		param =  lapply(rownames(param), 
			function(i) drop(as.matrix(param[i,c('range','variance','shape')]))),
		diagOf2ndDeriv =  lapply(1:nrow(param), 
			function(i) for2deriv$diagOf2ndDeriv[,i,drop=FALSE]),
		offDiagSecondDerivX = lapply(1:nrow(param), 
			function(i) for2deriv$offDiagSecondDerivX[,i,drop=FALSE]),
		MoreArgs = c(
			data[c('precTemplateMatrix','sparseTemplate','offsetMatrix','cholGmrfCorTemplate')],
			for2deriv['offDiagSecondDerivIJ']
			)
		)
	} else {

	twoDeriv = Map(derivDiag,
	param =  lapply(rownames(param), 
		function(i) drop(as.matrix(param[i,c('range','variance','shape')]))),
	diagOf2ndDeriv =  lapply(1:nrow(param), 
		function(i) for2deriv$diagOf2ndDeriv[,i,drop=FALSE]),
	offDiagSecondDerivX = lapply(1:nrow(param), 
		function(i) for2deriv$offDiagSecondDerivX[,i,drop=FALSE]),
	MoreArgs = c(
		data[c('precTemplateMatrix','sparseTemplate','offsetMatrix','cholGmrfCorTemplate')],
		for2deriv['offDiagSecondDerivIJ']
		)
	)
	}


	twoDeriv = do.call(cbind, twoDeriv)
	twoDeriv = twoDeriv / 2

	coarseCells = data$coarseCells

	seArray = array(twoDeriv, 
			c(ncol(coarseCells), nrow(coarseCells), 
				ncol(twoDeriv)))

	seBrick = raster::brick(seArray, 
			xmin(coarseCells), xmax(coarseCells),
			ymin(coarseCells), ymax(coarseCells), 
			projection(coarseCells),
			transpose=TRUE)
	names(seBrick) = colnames(twoDeriv)
	seBrick

}