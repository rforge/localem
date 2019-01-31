
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
				unlist(parallel::clusterCall(outerCluster, 
					function() Sys.getpid())), 
				'\n')
		} 
	} else { 
		# noCoresOuter <= 1
		outerCluster = NULL
		require('data.table')
	}

	if(nCoresInner > 1) {	
		if(nCoresOuter <= 1) { 
		# an inner cluster but no outer cluster
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

		 # inner and outer cluster
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
	# nCoresInner <= 1
		innerCluster = NULL
		if(!is.null(outerCluster)) {
			parallel::clusterEvalQ(outerCluster, {
				innerCluster = NULL
			})
		}
	}



	if(!is.null(outerCluster)) {
		res = parallel::clusterMap(outerCluster, 
#			localEM:::emsOneRange,
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

#			parallel::clusterExport(innerCluster, varlist = 'objectsForLikelihoodOneMap', envir = environment())
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


	if(reduce) {
		if(verbose){
			cat('computing quantiles\n')
		}


		if(nCoresOuter > 1) {
			clHere = outerCluster
		} else {
			clHere = innerCluster
		}

		pred = try(rgcpPred(res, cl=clHere, verbose=verbose), silent=TRUE) 
		res$pred = pred
	}

	if(nCoresOuter > 1) {
		if(nCoresInner > 1) {
			try(parallel::clusterEvalQ(outerCluster, parallel::stopCluster(innerCluster)))
		}
		try(parallel::stopCluster(outerCluster))
	} else {
		if(nCoresInner > 1) {
			try(parallel::stopCluster(innerCluster))	
		}
	}

	res
}


rgcpPred = function(
	x,
	sd,
	param = x$mle,
	theta = x$theta,
	data = x$data,
	p= c(0.01,0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975,0.99),
	cl = NULL,
	verbose=FALSE
	) {

	if(is.character(cl)) {
		cl = get(cl)
	}

	if(is.numeric(cl)) {
		cl = parallel::makeCluster(cl, type='PSOCK', methods=TRUE)
		parallel::clusterEvalQ(cl, require('Matrix'))
		parallel::clusterEvalQ(cl, require('data.table'))
		parallel::clusterEvalQ(cl, require('localEM'))

		stopCluster = TRUE

	} else {
		stopCluster = FALSE
	}


	if(!missing(sd) & 'mean' %in% names(x)) {
		varBrick = sd^2
		meanBrick = theta + varBrick
	} else {
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


		if(!is.null(cl) & (nrow(param)<1) ) {

			twoDerivList = parallel::clusterMap(cl,
				derivDiag,
				param =  lapply(rownames(param), 
					function(i) drop(as.matrix(param[i,c('range','variance','shape')]))),
				diagOf2ndDeriv =  lapply(1:nrow(param), 
					function(i) for2deriv$diagOf2ndDeriv[,i,drop=FALSE]),
				offDiagSecondDerivX = lapply(1:nrow(param), 
					function(i) for2deriv$offDiagSecondDerivX[,i,drop=FALSE]),
				MoreArgs = c(
					data[c('precTemplateMatrix','sparseTemplate','offsetMatrix')],
					for2deriv['offDiagSecondDerivIJ']
					)
				)
		} else {

			twoDerivList = Map(derivDiag,
				param =  lapply(rownames(param), 
					function(i) drop(as.matrix(
						param[i,c('range','variance','shape')]))),
				diagOf2ndDeriv =  lapply(1:nrow(param), 
					function(i) for2deriv$diagOf2ndDeriv[,i,drop=FALSE]),
				offDiagSecondDerivX = lapply(1:nrow(param), 
					function(i) for2deriv$offDiagSecondDerivX[,i,drop=FALSE]),
				MoreArgs = c(
					data[c('precTemplateMatrix','sparseTemplate','offsetMatrix')],
					for2deriv['offDiagSecondDerivIJ'],
					verbose=verbose
					)
				)
		}


		diagTwoDeriv = do.call(cbind, twoDerivList)

		coarseCells = data$coarseCells

		varArray = array(diagTwoDeriv/2, 
			c(ncol(coarseCells), nrow(coarseCells), 
				ncol(diagTwoDeriv)))

		varBrick = raster::brick(varArray, 
			xmin(coarseCells), xmax(coarseCells),
			ymin(coarseCells), ymax(coarseCells), 
			projection(coarseCells),
			transpose=TRUE)
		names(varBrick) = names(theta)

		meanBrick = theta + varBrick
	}

	ncp = theta^2/varBrick
	meanNcp = mean(values(ncp), na.rm=TRUE)
	SqForApprox = sort(unique(c(
		seq(0, meanNcp, len=40),
		seq(meanNcp, 2*meanNcp, len=40),
		quantile(values(ncp), seq(0,1,len=201) ))))

	Squant = p

	if(!is.null(cl)) {
		pchisqList = parallel::clusterMap(cl, 
			stats::pchisq,
			q = SqForApprox,
			MoreArgs = list( 
				df= 1, 
				ncp = values(ncp),
				log.p=TRUE, lower.tail=FALSE) 
			)

	} else {

		pchisqList = Map(stats::pchisq,
			q = SqForApprox,
			MoreArgs = list( 
				df= 1, 
				ncp = values(ncp),
				log.p=TRUE, lower.tail=FALSE) 
			)
	}

	if(stopCluster) {
		parallel::stopCluster(cl)
	}

	if(!is.matrix(pchisqList[[1]])) {
		pchisqMat = do.call(cbind, pchisqList)
	} else {
		pchisqMat = do.call(abind::abind, c(pchisqList, list(along=3)))
	}

	pchisqMat[is.nan(pchisqMat)] = -Inf
	pchisqMat = exp(pchisqMat)

	qchisqMat = apply(pchisqMat, 
			seq(1, length(dim(pchisqMat))-1), 
			function(xx) {
				if(length(unique(xx))==1) {
					xx[1] = xx[1]+1
				}
				approx(xx, SqForApprox, 1-Squant, rule=2)$y
			})
	if(is.matrix(qchisqMat)) {
		qchisqMat= t(qchisqMat)
	} else {
		qchisqMat = aperm(qchisqMat, c(2,3,1))
	}

#	dimnames(qchisqMat)[[2]] = paste0('q', Squant)

if(FALSE) {		qchisqBrick = raster::brick(
			array(qchisqMat, c(ncol(meanBrick), nrow(meanBrick),
				prod(dim(qchisqMat)[-1])))
			)
		names(qchisqBrick) = apply(
			do.call(expand.grid, dimnames(qchisqMat)[-1]),
			1, paste, collapse='.')

	qchisqMat2= as.vector(qchisqMat[,,1]) * as.vector(values(varBrick[[1]]))
	qBrick = raster::brick(array(
		qchisqMat2,
		c(ncol(meanBrick), nrow(meanBrick), length(Squant))
		), transpose=TRUE)

}

#	stuff = array(as.vector(values(varBrick)), dim(qchisqMat))
	qchisqMat2 = as.vector(qchisqMat) * as.vector(drop(values(varBrick)))

#	qq2 = array(qchisqMat2, dim(qchisqMat), dimnames=dimnames(qchisqMat))

	qBrick = raster::brick(array(
		qchisqMat2,
		c(ncol(meanBrick), nrow(meanBrick), length(Squant)*nlayers(theta))
		), transpose=TRUE)
	extent(qBrick) = extent(meanBrick)
	projection(qBrick) = projection(meanBrick)


	seBrick = sqrt(varBrick)
	names(seBrick) = paste0(names(varBrick), '.se')

	if(nlayers(theta) == 1) {
		names(qBrick) = paste0('q', Squant)
		names(theta) = 'mode'
		names(meanBrick) = 'mean'
		names(seBrick) = 'se'
	} else {
	names(qBrick) = apply(
		expand.grid(
			names(theta), 
			paste0('q', Squant)),
		1, paste, collapse='.'
		)
		names(theta) = paste0(names(theta), '.mode')
		names(meanBrick) = paste0(names(meanBrick), '.mean')
	}


	resQ = stack(
		meanBrick,
		seBrick,
		qBrick,
		theta
		)

	resQ

}


emsExcProb = function(x, threshold=1) {

	modeLayers = grep('(^|[.])mode$', names(x), value=TRUE)
	seLayers = grep("se$", names(x), value=TRUE)

	if(length(seLayers)>1) {
		seLayers = gsub('mode$', 'se', modeLayers)
	}

	vars = values(x[[seLayers]]^2)

	ncp = values(x[[modeLayers]])^2 / vars

	resE = mapply(
		function(threshold, vars,ncp) {
			stats::pchisq(threshold/vars,1,ncp,lower.tail=FALSE)
		},
		threshold = threshold,
		MoreArgs = list(vars = vars, ncp = ncp), SIMPLIFY=FALSE
		)
	if(length(modeLayers)>1) {
		resE = do.call(abind::abind, c(resE, list(along=3)))
	} else {
		resE = do.call(cbind, resE)
	}

	# in case resE is a vector
	coarseCells = raster(x)
	eArray = array(resE, 
		c(ncol(coarseCells), nrow(coarseCells), 
			prod(dim(resE)[-1])))
	eBrick = raster::brick(eArray, 
		xmin(coarseCells), xmax(coarseCells),
		ymin(coarseCells), ymax(coarseCells), 
		projection(coarseCells),
		transpose=TRUE)
	names(eBrick) = as.vector(outer(
		gsub("[.]?mode$", "", modeLayers),
		threshold, paste, sep='exc'))
	eBrick

}
