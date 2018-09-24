emsObjects = function(y, x, shape,
	obs = grep("^count", names(y[[1]]), value=TRUE), 
	intercept = 1) {

	Sobs = unique(obs)
	Sintercept = sort(unique(intercept))

	Oijl = x$Oijl
	coarseCells = x$cells$coarse
	Ncells = ncell(coarseCells)

	y = lapply(y, function(xx) {
		res = as.matrix(xx@data[,Sobs])
		rownames(res) = xx$id
		res
	})

	Smap = names(y)

	YijOijl = mStepMat = list()
	for(Dmap in Smap) {
		res = mapply(function(Dobs) {
			res = Oijl[[Dmap]] %*% 
			Matrix::Diagonal(nrow(y[[Dmap]]), y[[Dmap]][,Dobs])
			res = as(res, 'TsparseMatrix')
			cbind(i=res@i, j=res@j, x=res@x)
		}, Dobs = Sobs, SIMPLIFY=FALSE)
		YijOijl[[Dmap]] = cbind(
			res[[1]][,c('i','j')],
			cell = res[[1]][,'i']+1,
			id = res[[1]][,'j']+1,
			do.call(cbind, lapply(res, function(xx) xx[,'x']))
			)

		mStepMat[[Dmap]] = Matrix::Matrix(outer(
			1:ncell(coarseCells), YijOijl[[Dmap]][,'cell'], 
			FUN = '=='
			))
	}

# replicate datasets for use with different intercepts
	obsIntercept = expand.grid(obs = Sobs, intercept = Sintercept)
	rownames(obsIntercept) = apply(
		obsIntercept[,c('obs','intercept')], 1, paste, collapse='_alpha')

	interceptMatrix = t(matrix(
		obsIntercept$intercept,
		ncol=ncell(coarseCells), nrow=nrow(obsIntercept),
		dimnames = list(
			rownames(obsIntercept), values(coarseCells)
			)
		))

	offsetIntercept = x$offset$matrix %*% interceptMatrix


	
# template for the determinant of second derivative
	if(!missing(shape)) {

	precTemplateMatrix = geostatsp::NNmat(
		coarseCells, 
		nearest = shape[1]+1, 
		adjustEdges=TRUE)

	gmrfCorMatrix = geostatsp::maternGmrfPrec(
		precTemplateMatrix, 
		param=c(shape=unname(shape), variance=1, range = 10*mean(res(coarseCells))),
		adjustEdges = TRUE
		)


	Dobs = rownames(obsIntercept)[1]
	Ddataset = as.character(obsIntercept[Dobs,'obs'])
	sumDenomBoth = YlogSumBoth = 
	outerOffsetBoth = diagBoth = list()
	SobsIntercept = Dobs

	for2deriv = objectsForLikeilhood(
		Oijl, 
		y = lapply(y, function(xx) xx[,1, drop=FALSE]), 
		lambda = matrix(1, Ncells, 1,
			dimnames = list(NULL, rownames(obsIntercept)[1] )
		))
	secondDeriv <- get2ndDeriv(
		diagCombined = for2deriv$diagOf2ndDeriv,
		offDiagSecondDerivIJ=for2deriv$offDiagSecondDerivIJ,
		offDiagSecondDerivX=for2deriv$offDiagSecondDerivX,		
		precPlusTwoOffset = gmrfCorMatrix + 2*x$offset$matrix
		)

	secondDeriv = as.matrix(secondDeriv)
	sparseTemplate1 = Matrix::sparseMatrix(
		i = secondDeriv[,'i'],
		j = secondDeriv[,'j'],
		x = 1:nrow(secondDeriv),
		dims = rep(Ncells, 2),
		index1=FALSE
		)
	sparseTemplate = Matrix::forceSymmetric(sparseTemplate1)

	secondDerivHere = sparseTemplate
	secondDerivHere@x = secondDeriv[sparseTemplate@x, Dobs]	

	chol2ndDeriv = try(Matrix::Cholesky(secondDerivHere, LDL=TRUE), silent=TRUE)
	cholGmrfCor = Matrix::Cholesky(gmrfCorMatrix, LDL=TRUE)


	if(class(chol2ndDeriv)=='try-error') {
		if(verbose) {
			cat('!')
		}
		chol2ndDeriv = NULL
		warning("can't cholesky 2nd deriv matrix")
	}
	} else { # shape not provided
		cholGmrfCor = chol2ndDeriv = shape = precTemplateMatrix = NULL
	}


	list(y=y, YijOijl=YijOijl, 
		interceptMatrix = interceptMatrix,
		offsetIntercept = offsetIntercept,
		obsIntercept = obsIntercept,
		Oijl = Oijl, 
		mStepMat = mStepMat,
		offsetMatrix = x$offset$matrix,
		coarseCells = x$cells$coarse,
		sparseTemplate=sparseTemplate,
		shape = shape,
		precTemplateMatrix = precTemplateMatrix,
		chol2derivTemplate = chol2ndDeriv,
		cholGmrfCorTemplate = cholGmrfCor)
}

