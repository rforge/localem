

get2ndDeriv = function(
	diagCombined,
	offDiagSecondDerivIJ,
	offDiagSecondDerivX,		
	precPlusTwoOffset,
	sparseTemplate		
	) {

	# diagonal elemts 2 sum(Y_ij O_ijl)
	cellSeq = seq(0, len=nrow(diagCombined))
	diagCombined = cbind(i=cellSeq, j=cellSeq, 
		x=(-2)*as.matrix(diagCombined))

	offDiagCombined = cbind(
		as.matrix(offDiagSecondDerivIJ),
		4*as.matrix(offDiagSecondDerivX)
		)

	# 4 * sum Y_ij O_ijl theta_ijl O_ijk theta_ijk

	# the correlation matrix and offsets
	precPlusTwoOffsetT = as(precPlusTwoOffset, 'TsparseMatrix')
	precPlusTwoOffsetL =	cbind(
		i=precPlusTwoOffsetT@i, j=precPlusTwoOffsetT@j,
		matrix(precPlusTwoOffsetT@x, 
			length(precPlusTwoOffsetT@x), ncol(	offDiagSecondDerivX),
			dimnames = list(NULL, colnames(	offDiagSecondDerivX)))
		)

	# add up diag, outer offsets, and gmrf mat
	toAgg = data.table::as.data.table(rbind(
		diagCombined, 
		offDiagCombined, 
		precPlusTwoOffsetL
		))

	secondDerivLong = 
		toAgg[,  lapply(.SD, sum), by = .(i,j)]


	if(!missing(sparseTemplate)) {
		st2 = as(sparseTemplate, 'TsparseMatrix')
		secondDerivLong = merge(
			data.table(
				i = st2@i, j=st2@j, index=st2@x
				), secondDerivLong, sort=FALSE, all=TRUE)
		res = secondDerivLong[
				match(secondDerivLong$index, 
					sparseTemplate@x),]
		if(FALSE) {
			# check
			myMat = forceSymmetric(sparseMatrix(
				res[,'i'], res[,'j'], x=res[,4],
				index1=FALSE
				))
			range(diag(expand(Cholesky(myMat))$L))
			myMat2 = sparseTemplate
			myMat2@x = res[,4]
			range(diag(expand(Cholesky(myMat2))$L))
		}
	} else {
		res = secondDerivLong[order(i,j)]
	}
	res
}

derivDet = function(outerOffsetHere, 
	sparseTemplate, 
	cholTemplate,
	verbose) {

	secondDerivHere = sparseTemplate
	secondDerivHere@x = outerOffsetHere	

	cholHere <- try(Matrix::update(
		cholTemplate, secondDerivHere), silent=TRUE)
	detHere = try(Matrix::determinant(cholHere, 
		logarithm=TRUE)$modulus, silent=TRUE)
	if(class(detHere)=='try-error') {
		if(verbose) {
			cat('!')
		}
		bob = try(
			eigen(as.matrix(secondDerivHere), only.values=TRUE)$values,
			silent=TRUE)
		if(class(bob) != 'try-error') {
			if(any(bob < 0)) {
				warning("negative eigen values, setting to a small positive number")
			}
			bob = pmax(0.001, bob)
			detHere = sum(log(bob))/2
		} # eiten not try error
	} # chol not try error
	detHere
}

objectsForLikelihoodOneMap = 
function(Dmap, Oijl, yHere, 
	lambdaHere, thetaHere = lambdaHere^2) {

	# sum_m O_ijm theta_m^2
	offThetaIJ = as.vector(Matrix::crossprod(Oijl[[Dmap]], lambdaHere))

	# sum_ij Y_ij O_ijl / sum_m O_ijm theta_m^2
	diagOf2ndDeriv = apply(Oijl[[Dmap]] %*% Matrix::Diagonal(
		length(offThetaIJ), yHere[[Dmap]]/offThetaIJ),
	1,sum)

	# sum_ij Y_ij O_ijl theta_l O_ijk theta-k / [sum_m O_ijm theta_m^2]^2
	offDiagSecondDeriv = 
	Matrix::tcrossprod(Diagonal(length(lambdaHere), lambdaHere) %*%
		Oijl[[Dmap]] %*% Matrix::Diagonal(length(offThetaIJ), 
			sqrt(yHere[[Dmap]])/offThetaIJ))

	list(
		diagOf2ndDeriv = diagOf2ndDeriv,
		offDiagSecondDeriv = offDiagSecondDeriv,
		dPoisson = sum(yHere[[Dmap]]*log(offThetaIJ) ) -
		sum(offThetaIJ)
		)		
}

objectsForLikelihoodOneDataset = 
function(yHere, 
	lambdaHere, 
	Oijl, 
	thetaHere = sqrt(lambdaHere)) {

	res = mapply(objectsForLikelihoodOneMap,
		Dmap = names(yHere),
		MoreArgs = list(
			yHere = yHere,
			lambdaHere = lambdaHere,
			thetaHere = thetaHere, 
			Oijl = Oijl
			), SIMPLIFY=FALSE)
	res = mapply(
		function(D, yy) lapply(yy, function(yyy) yyy[[D]]),
		D = names(res[[1]]),
		MoreArgs = list(yy=res),
		SIMPLIFY = FALSE
		)

	res$dPoisson = sum(unlist(res$dPoisson))
	res$diagOf2ndDeriv = apply(do.call(cbind,res$diagOf2ndDeriv),1,sum)
	res$offDiagSecondDeriv = do.call(rbind, 
		lapply(res$offDiagSecondDeriv, 
			function(xx) {
				xxx = as(xx, 'TsparseMatrix')
				cbind(i=xxx@i, j=xxx@j, x=xxx@x)
			} 
			))
	res
}

objectsForLikeilhood = function(Oijl, y, lambda, cl = NULL) {

	yExpanded = mapply(
		function(Dobs, yy) lapply(yy, function(yyy) yyy[,Dobs]),
		Dobs = colnames(y[[1]]),
		MoreArgs = list(yy=y),
		SIMPLIFY = FALSE
		)

	if(!is.null(cl)) {
		res = parallel::clusterMap(cl,
			objectsForLikelihoodOneDataset,
			lambdaHere = as.list(as.data.frame(as.matrix(lambda))),
			yHere = yExpanded[gsub('_.*', '', colnames(lambda))],
			MoreArgs = list( Oijl = Oijl))
	} else {

		res = Map(
			objectsForLikelihoodOneDataset,
			lambdaHere = as.list(as.data.frame(as.matrix(lambda))),
			yHere = yExpanded[gsub('_.*', '', colnames(lambda))],
			MoreArgs = list( Oijl = Oijl)
			)
	}
	res2 = list(
		diagOf2ndDeriv = do.call(cbind, lapply(res, function(xx) xx$diagOf2ndDeriv)),
		dPoisson = unlist(lapply(res, function(xx) xx$dPoisson)),
		offDiagSecondDerivIJ = 
		res[[1]]$offDiagSecondDeriv[,c('i','j')],
		offDiagSecondDerivX = do.call(cbind, 
			lapply(res, function(xx) xx$offDiagSecondDeriv[,'x']))
		)
	res2
}

eStepFun = function(
	Yfull, OijlHere, YijOijlHere, 
	mStepMatHere, SobsFull, lambda) {

	denom = Matrix::crossprod(OijlHere, lambda)

	mStepMatHere %*% (
		YijOijlHere[,SobsFull] / denom[YijOijlHere[,'id'],]
		)
}




#' 

#+ functionForEms, include=FALSE
emsOneSd = function(
	sd,
	gmrfCorMatrix,
	data,
	theta,
	maxIter = 2000,
	iterPrint = 10,
	tol = .Machine$double.eps^(1/3),
	cl = NULL,
	reduce = FALSE,
	verbose=FALSE
) {

	Dsd = sd[1]
	DtwoSigmaSquared = Dsd^2*2
	DSigmaSquared = Dsd^2

	sparseTemplate = data$sparseTemplate
	cholTemplate =  data$chol2derivTemplate


	YijOijl = data$YijOijl
	offsetIntercept = data$offsetIntercept
	obsIntercept = data$obsIntercept
	interceptMatrix = data$interceptMatrix
	y = data$y
	Oijl = data$Oijl
	offsetMat = data$offsetMatrix
	coarseCells = data$coarseCells
	twoOffsetMat = 2*offsetMat

	SobsIntercept = colnames(offsetIntercept)
	SobsFull = as.character(obsIntercept[,'obs'])

	Ny = length(SobsIntercept)
	Ncells = nrow(data$offsetMatrix)
	Scells = 1:Ncells
	Smap = names(y)

	StildeInv = DtwoSigmaSquared * offsetMat + gmrfCorMatrix
	StildeChol = Matrix::update(data$cholGmrfCorTemplate, StildeInv)
	
	etaStuff = interceptMatrix - DtwoSigmaSquared * 
	Matrix::solve(StildeChol, offsetIntercept)

	# loop through iterations

	# starting values
	if(missing(theta)) {
	theta = matrix(1, Ncells, Ny,
		dimnames = list(
			values(coarseCells),
			SobsIntercept))
	}
	lambda = theta^2

	if(verbose) {
		cat("sd ", signif(Dsd,2), ' ')
	}


	if(!is.null(cl)) {
		mapFun = function(...) parallel::clusterMap(cl, ...)
	} else {
		mapFun = Map
	}	
	
	Diter = 0
	derivMax = Inf
	while(Diter< maxIter & derivMax > tol) {

		muLambda1 = mapFun(
			eStepFun, 
			YijOijlHere = data$YijOijl,
			OijlHere = data$Oijl,
			mStepMatHere = data$mStepMat,
			MoreArgs = list(
				SobsFull=SobsFull,
				lambda = lambda)
			)
		muLambda = Reduce('+', muLambda1)
		muLambdaTheta = theta * muLambda

		muStuff <- DtwoSigmaSquared * theta * 
			Matrix::solve(StildeChol, muLambdaTheta)

		# update lambda and theta
		lambda = muStuff + theta * etaStuff
		theta = sqrt(as.matrix(lambda))

		firstDerivHere = 2 * muLambdaTheta - 
			twoOffsetMat %*% theta -
			(1/DSigmaSquared) * (
				gmrfCorMatrix %*% (theta  - interceptMatrix))

		derivMax = max(abs(firstDerivHere))

		if(verbose & (!(Diter %% iterPrint))) {
			cat(signif(derivMax, 2), ' ')
		}
		Diter = Diter + 1
	} # Diter
	if(verbose) {
		cat(Diter, derivMax)
	}
	if(Diter >= maxIter) {
		warning("max iterations reached")
	}
	# second derivative
	if(verbose) {
		cat(' det ')
	}

	for2deriv = objectsForLikeilhood(Oijl, y, lambda, cl = cl) 


	if(verbose) {
		cat(".")
	}

	secondDeriv <- get2ndDeriv(
		diagCombined = for2deriv$diagOf2ndDeriv,
		for2deriv$offDiagSecondDerivIJ,
		for2deriv$offDiagSecondDerivX,		
		precPlusTwoOffset = StildeInv/DSigmaSquared,
		sparseTemplate)		
	

	if(verbose) {
		cat('.')
	}


		halfLogDet <- mapFun(
			derivDet,
			# secondDeriv is a data.table, need with=FALSE
			outerOffsetHere = as.list(
					secondDeriv[,SobsIntercept, with=FALSE]),
			MoreArgs = list(
				sparseTemplate = sparseTemplate,
				cholTemplate = cholTemplate,
				verbose = verbose
				)
			)

	if(verbose) {
		cat(".")
	}
	halfLogDet = unlist(halfLogDet)		

	#determinant(gmrfCorMatrix)$modulus/2

	thetaC = theta - interceptMatrix 
	thetaCC = apply(thetaC * (gmrfCorMatrix %*% thetaC),2,sum)


	theL = cbind(
		obsIntercept,
		thetaCrossprod = thetaCC,
		halfLogDet = halfLogDet,
		dPoisson = for2deriv$dPoisson,
		sd = rep(Dsd, Ny)#length(thetaCC))
		)

	theL$logLminusDet = 
		theL[,'dPoisson'] -
		thetaCC/DtwoSigmaSquared - halfLogDet -
		Ncells * log(Dsd)


	if(reduce) {
		# not written yet
	}

	if(verbose) {
		cat("\n")
	}
	list(logLik = theL, theta = theta)
}

emsOneRange = function(
	range,
	Ssd,
	data,
	shape,
	maxIter = 2000,
	tol = .Machine$double.eps^(1/3),
	reduce=TRUE,
	verbose = FALSE,
	cl = NULL
	) {

	if(is.character(cl)) {
		cl = get(cl)
	}

	if(verbose) {
		cat('pid ', Sys.getpid(), ' range ', range, ' ')
	}

	templateMatrix = data$precTemplateMatrix

	if(is.null(templateMatrix)) {
		if(missing(shape)) {
			stop("one of shape and templateMatrix must be specified")
		}
		templateMatrix = geostatsp::NNmat(
			data$coarseCells, 
			nearest = shape+1, 
			adjustEdges=TRUE)
	} else {
		if(!missing(shape)) {
			if(shape != (attributes(templateMatrix)$nearest-1 ) ) {
				warning('ignoring shape argument')
			}
		}
		shape = attributes(templateMatrix)$nearest-1
	}

	gmrfCorMatrix = geostatsp::maternGmrfPrec(
		templateMatrix,
		param = c(
			variance=1, 
			range = as.numeric(range[1]),
			shape = as.numeric(shape[1])),
		adjustEdges=TRUE)



	theta = matrix(1, nrow(data$offsetMatrix), nrow(data$obsIntercept),
		dimnames = list(
			1:nrow(data$offsetMatrix),
			rownames(data$obsIntercept)))
	resAllSd = list()


	for(Dsd in Ssd) {
		resAllSd[[as.character(Dsd)]] = emsOneSd(
			sd = Dsd, gmrfCorMatrix, data,         
			theta = theta, maxIter = maxIter, 
			iterPrint = 10,
			tol=tol, cl = cl, 
			verbose = verbose>1
			)
		if(verbose>0) {
			cat(' ', range[1], '_', Dsd, ' ')
		}
		theta = resAllSd[[as.character(Dsd)]]$theta
		if(reduce) {
			logLik2 = resAllSd[[as.character(Dsd)]]$logLik
			logLik2$name = rownames(logLik2)
			onlyMax = merge(logLik2,
				aggregate(logLik2[,'logLminusDet',drop=FALSE], 
					logLik2[,'obs', drop=FALSE], max),
				all=FALSE)
			resAllSd[[as.character(Dsd)]]$theta = 
			resAllSd[[as.character(Dsd)]]$theta[,onlyMax$name]
			colnames(resAllSd[[as.character(Dsd)]]$theta) = 
			gsub('_.*', '', colnames(resAllSd[[as.character(Dsd)]]$theta))
		}
	}

	theLList = lapply(resAllSd, function(xx) xx$logLik)
	theThetaList = lapply(resAllSd, function(xx) xx$theta)

	logLik = do.call(rbind, theLList)
	rownames(logLik) = 1:nrow(logLik)


	if(reduce) {
		onlyMax = merge(
			logLik,
			aggregate(logLik[,'logLminusDet',drop=FALSE], logLik[,'obs', drop=FALSE], max),
			all=FALSE)
		onlyMax$Dsd = match(onlyMax$sd, Ssd)

		thetaSub = list()
		for(Dobs in 1:nrow(onlyMax)) {
			thetaSub[[as.character(onlyMax[Dobs,'obs'])]] = 
			theThetaList[[
			onlyMax[Dobs,'Dsd']
			]][ , as.character(onlyMax[Dobs,'obs'])]
		}
		thetaAll = do.call(cbind, thetaSub)
		if(ncol(thetaAll) == nrow(onlyMax)) {
			# it's possible the grep command failed and selected zero or 2 columns
			colnames(thetaAll) = onlyMax$obs
		}
	} else {
		thetaAll = do.call(abind::abind, c(theThetaList, list(along=3)))
	}
	# combine results for different SD

	logLik$mean = logLik$intercept^2
	logLik$sdNatrual = logLik$sd/logLik$intercept
	logLik$range = range[1]

	cholGmrfCor = Matrix::update(data$cholGmrfCorTemplate, gmrfCorMatrix)
	logLik$invCorMatHalfLogDet = 
		Matrix::determinant(cholGmrfCor)$modulus
	logLik$logL = logLik$logLminusDet + logLik$invCorMatHalfLogDet

	if(verbose) {
		cat(' ', range, ' ')
	}


	list(logL = logLik, theta = thetaAll)
}
