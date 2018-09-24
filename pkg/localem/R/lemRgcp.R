
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


if(FALSE) {
nCoresOuter = 2; nCoresInner = 3
outerCluster = parallel::makeCluster(2, type='PSOCK', methods=TRUE)
parallel::clusterEvalQ(outerCluster, bob <- Sys.getpid())
parallel::clusterEvalQ(outerCluster, get('bob'))

parallel::clusterEvalQ(outerCluster, 
	{innerCluster <- parallel::makeCluster(3, type='PSOCK', methods=TRUE)})

parallel::clusterEvalQ(outerCluster, ls())

parallel::clusterEvalQ(outerCluster, parallel::clusterEvalQ(innerCluster, Sys.getpid()))

parallel::clusterEvalQ(outerCluster, parallel::stopCluster(innerCluster))
parallel::stopCluster(outerCluster)
}


allFunctions = c('derivDet', 'eStepFun',
	'objectsForLikelihoodOneMap',
	'objectsForLikelihoodOneDataset',
	'objectsForLikeilhood', 
	'get2ndDeriv',
    'emsOneRange',
    'emsOneSd')


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
	if(nCoresOuter == 1) {
		warning("set outer cores > 1 when inner cores > 1")
	}
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
		    cl = 'innerCluster',
		    reduce = reduce,
		    verbose = 2*verbose
			)) 
}

if(nCoresInner > 1) {
	parallel::clusterEvalQ(outerCluster, parallel::stopCluster(innerCluster))
}
if(nCoresOuter > 1) {
	parallel::stopCluster(outerCluster)
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
	list(logL = logLik, profL = logLprof, theta = thetaBrick, array = logLikArray)
}

