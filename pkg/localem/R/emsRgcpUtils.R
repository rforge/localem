

get2ndDeriv = function(
	diagCombined,
	offDiagSecondDerivIJ,
	offDiagSecondDerivX,		
	precPlusTwoOffset,
	sparseTemplate		
	) {

	# diagonal elemts 2 sum(Y_ij O_ijl)
	cellSeq = seq(0, len=nrow(diagCombined))
	diagCombined2 = cbind(
		data.frame(i=cellSeq), 
		data.frame(j=cellSeq), 
		(-2)*diagCombined)

	offDiagCombined = cbind(
		as.data.frame(offDiagSecondDerivIJ),
		4*offDiagSecondDerivX
		)

	# 4 * sum Y_ij O_ijl theta_ijl O_ijk theta_ijk

	# the correlation matrix and offsets
	precPlusTwoOffsetT = as(precPlusTwoOffset, 'TsparseMatrix')
	precPlusTwoOffsetL =	cbind(
		data.frame(i=precPlusTwoOffsetT@i), 
		data.frame(j=precPlusTwoOffsetT@j),
		matrix(precPlusTwoOffsetT@x, 
			length(precPlusTwoOffsetT@x), 
			ncol(offDiagSecondDerivX),
			dimnames = list(NULL, 
				colnames(	offDiagSecondDerivX)))
		)

	# add up diag, outer offsets, and gmrf mat
	toAgg = data.table::rbindlist(list(
		diagCombined2, 
		offDiagCombined, 
		precPlusTwoOffsetL
		))

	secondDerivLong = 
		toAgg[ ,  
			lapply(.SD, sum), 
			by = list(i,j)]
	

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

derivDiag = function(
		param, offsetMatrix,precTemplateMatrix,
		diagOf2ndDeriv,offDiagSecondDerivIJ,
		offDiagSecondDerivX, sparseTemplate,
		cholGmrfCorTemplate) {

		precMat = offsetMatrix + 
			geostatsp::maternGmrfPrec(
				precTemplateMatrix,
				param = param, 
				adjustEdges = TRUE
				)

		derivHere = get2ndDeriv(
			diagCombined = diagOf2ndDeriv,
			offDiagSecondDerivIJ=offDiagSecondDerivIJ,
			offDiagSecondDerivX=offDiagSecondDerivX,		
			precPlusTwoOffset = precMat,
			sparseTemplate)
		derivHere = as.data.frame(derivHere)

		derivMat = Matrix::sparseMatrix(
			i = as.vector(derivHere[,'i']),
			j = as.vector(derivHere[,'j']),
			x = as.vector(derivHere[,4]),
			dims = dim(precMat), dimnames = dimnames(precMat),
			index1 = FALSE)

		cholHere <- try(Matrix::update(
			cholGmrfCorTemplate, derivMat), silent=TRUE)

		diag(solve(cholHere))

	}

objectsForLikelihoodOneMap = function(
	OijlHere, yHere, 
	lambdaHere, thetaHere = sqrt(lambdaHere)) {

	# sum_m O_ijm theta_m^2
	offThetaIJ = as.vector(Matrix::crossprod(
		OijlHere, lambdaHere))

	
	# sum_ij Y_ij O_ijl / sum_m O_ijm theta_m^2
	diagOf2ndDeriv = apply(
		OijlHere %*% Matrix::Diagonal(
			length(offThetaIJ), yHere/offThetaIJ),
		1, sum)

	# sum_ij Y_ij O_ijl theta_l O_ijk theta-k / [sum_m O_ijm theta_m^2]^2
	offDiagSecondDeriv = Matrix::tcrossprod(
		Diagonal(
			length(lambdaHere), lambdaHere) %*%
		OijlHere %*% Matrix::Diagonal(
			length(offThetaIJ), 
			sqrt(yHere)/offThetaIJ))

	offDiagSecondDeriv =
		as(offDiagSecondDeriv, 'TsparseMatrix')
	offDiagSecondDeriv = data.frame(
		i = offDiagSecondDeriv@i,
		j = offDiagSecondDeriv@j,
		x = offDiagSecondDeriv@x
		)

	list(
		diagOf2ndDeriv = diagOf2ndDeriv,
		offDiagSecondDeriv = offDiagSecondDeriv,
		dPoisson = sum(yHere*log(offThetaIJ) ) -
		sum(offThetaIJ)
		)		
}

objectsForLikelihoodOneDataset = 
function(yHere, 
	lambdaHere, 
	Oijl, 
	thetaHere = sqrt(lambdaHere)) {

	res = Map(objectsForLikelihoodOneMap,
		OijlHere = Oijl,
		yHere = yHere,
		MoreArgs = list(
			lambdaHere = lambdaHere,
			thetaHere = thetaHere
			))

	res = mapply(
		function(D, yy) lapply(yy, function(yyy) yyy[[D]]),
		D = names(res[[1]]),
		MoreArgs = list(yy=res),
		SIMPLIFY = FALSE
		)

	res$dPoisson = Reduce('+', res$dPoisson)
	res$diagOf2ndDeriv = Reduce('+', 
		res$diagOf2ndDeriv)
	res$offDiagSecondDeriv = as.matrix(
		data.table::rbindlist( 
		res$offDiagSecondDeriv)
	)

	res
}

objectsForLikelihood = function(
	Oijl, y, 
	lambda, cl = NULL) {

	yExpanded = Map(
		function(Dobs, yy) {
			lapply(yy, function(yyy) yyy[,Dobs])
		},
		Dobs = colnames(y[[1]]),
		MoreArgs = list(yy=y)
		)

	if(!is.null(cl)) {
		res = parallel::clusterMap(cl,
			objectsForLikelihoodOneDataset,
			yHere = yExpanded[
				gsub('_.*', '', colnames(lambda))],
			lambdaHere = lapply(
				seq_len(ncol(lambda)), 
				function(D) lambda[, D]),
			MoreArgs = list( Oijl = Oijl))
	} else {

		res = Map(
			objectsForLikelihoodOneDataset,
			yHere = yExpanded[
				gsub('_.*', '', colnames(lambda))],
			lambdaHere = lapply(
				seq_len(ncol(lambda)), 
				function(D) lambda[, D]),
			MoreArgs = list( Oijl = Oijl)
			)
	}
	names(res) = colnames(lambda)
	res2 = list(
		diagOf2ndDeriv = do.call(cbind, 
			lapply(res, 
				function(xx) xx$diagOf2ndDeriv)),
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

	for2deriv = objectsForLikelihood(
		Oijl, y, 
		lambda, cl = cl) 


	if(verbose) {
		cat(".")
	}

	secondDeriv <- get2ndDeriv(
		diagCombined = for2deriv$diagOf2ndDeriv,
		offDiagSecondDerivIJ=for2deriv$offDiagSecondDerivIJ,
		offDiagSecondDerivX=for2deriv$offDiagSecondDerivX,		
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
		cat('pid ', Sys.getpid(), ' range ', range[1], '\n')
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


	verboseHere = verbose==1
	for(Dsd in Ssd) {
		resAllSd[[as.character(Dsd)]] = emsOneSd(
			sd = Dsd, gmrfCorMatrix, data,         
			theta = theta, maxIter = maxIter, 
			iterPrint = 10,
			tol=tol, cl = cl, 
			verbose = verbose>1
			)
		if(verboseHere) {
			cat(' ', range[1], '_', Dsd, '\n')
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
			resAllSd[[as.character(Dsd)]]$theta[,onlyMax$name, drop=FALSE]
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


	list(logL = logLik, theta = thetaAll)
}
