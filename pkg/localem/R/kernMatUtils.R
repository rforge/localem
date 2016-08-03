# Sets the size of sigma for Gaussian density kernel
# The default is 3 times sigma
focalWeightWithSize = function(x, bw, size=NULL) {
	
	if(!length(size))
		size = 3*bw
	
	focalWeight(
			x=x, d=c(bw, size),
			type='Gauss'
	)
	
}


# Computes focal weight matrix for a Gaussian density kernel with specified bandwidth and size of sigma
focalFromBw = function(
  bw, fine, ncores=1, minDim = 60, focalSize = NULL
){
 centreCell = 1+dim(fine)[1:2]
	
 theDim = 2*rep(max(centreCell),2)-1
 bigCentreCell = (theDim-1)/2
	
 bigMat = matrix(0,
   theDim[1],theDim[2]
 )
 names(bw) = gsub("[[:space:]]", "",
   format(bw, scientific=FALSE)
 )
	
 focalList = parallel::mcmapply(
   focalWeightWithSize,
   bw=bw,
   MoreArgs=list(
     x=fine,
					size = focalSize
   ),
	 mc.cores=ncores, SIMPLIFY=FALSE
 )
	
 names(focalList) = paste("bw",
   names(bw),
   sep="")
	
 bigList = parallel::mcmapply(
   function(Dbw){
    theDim = dim(focalList[[Dbw]])
    topleft = bigCentreCell - (theDim-1)/2
    xseq = seq(topleft[1], len=theDim[1], by=1)
    xNotZero = (xseq > 0) & (xseq <= dim(bigMat)[1])
    yseq = seq(topleft[2], len=theDim[2], by=1)
    yNotZero = (yseq > 0) & (yseq <= dim(bigMat)[2])
    res = bigMat
    res[
      xseq[xNotZero], yseq[yNotZero]
    ] = focalList[[Dbw]][xNotZero,yNotZero]
    Matrix(res)
   },
   Dbw=names(focalList),
   mc.cores=ncores, SIMPLIFY=FALSE
 )
	
 bw = data.frame(
   bw=bw,
   dim=unlist(lapply(focalList,
       function(qq) min(dim(qq))))
 )
	
 bw$factLog2 = floor(pmax(0,log2(bw$dim) - log2(minDim)))
 bw$fact = 2^bw$factLog2
	
 result = list(focal=focalList,
   extended=bigList,centreCell=bigCentreCell,
   bw=bw)
	
 for(D in sort(setdiff(unique(bw$fact),1))){
  fineAgg = aggregate(fine, fact=D)
  bwHere = sort(bw[bw$fact==D,'bw'])
  focalListD =
    mapply(
      focalWeightWithSize,
      bw=bwHere,
      MoreArgs=list(
        x=fineAgg,
		size = focalSize
      )
    )
  names(focalListD) = paste("bw",
    bwHere,
    sep="")
  result[[paste('focalAgg', D,sep='')]] = focalListD
 }
	
 return(result)
}


# Reorders the focal weight matrix to estimate kernel smoothing function for specified bandwidth
# Utilises the symmetrical design of the focal weight matrix
reorderCellsTranslate = function(Ncell){
	
 Ncell = Ncell[1]
 NcellSq = Ncell^2
 cellSeq = 1:NcellSq
 cellMat = matrix(cellSeq, Ncell,Ncell)
 reorderRows = as.vector(cellMat[,Ncell:1])
 reorderBoth = rev(cellSeq)
 switchRowCol = as.vector(t(cellMat))
	
 list(
   q2=cellSeq[switchRowCol],
   q3=reorderRows[switchRowCol],
   q4=reorderRows
 )
}


# Computes the kernel smoothing function for specified bandwidth
kernMat = function(
  cellDist, focalList,
  fineRelativeRes,
  cell1=c(0,0),
  coarse=NULL,
  reorder=FALSE
) {
 if(length(cell1)==2 & length(cellDist)==2){
  cellDist = cellDist - cell1
 } else if(length(cellDist)== 1 & length(cell1) == 1) {
  if(is.null(coarse)) {
   warning("supply coarse if cellDist is a cell ID")
  }
  cellCoord = rowColFromCell(coarse, c(cell1,cellDist))
  cellDist = apply(cellCoord,2,diff)
 } else {
  warning('cellDist and cel1 should be both either length 1 (cell id) or length 2 (x,y)')
 }
	
 if(reorder){
  cellDist = c(max(abs(cellDist)), min(abs(cellDist)))
 }
	
 Sbw = names(focalList$extended) # names of bandwidths
	
 if(!is.vector(fineRelativeRes)) {
  if(is.null(coarse)) {
   warning("supply coarse if fineRelativeRes is a raster")
  }
  fineRelativeRes = round(res(coarse)/res(fineRelativeRes))
 }
 cellDistFine = cellDist * fineRelativeRes
	
 # construct array to hold result
 Ncell = round(prod(fineRelativeRes))
 kernelArray = array(NA,
   c(Ncell, Ncell, length(Sbw)),
   dimnames=list(NULL, NULL, Sbw)
 )
	
 # loop through fine cells in coarse cell 1
	
 for(Dcol in 1:(fineRelativeRes[2])){
  firstCellInCol = fineRelativeRes[2]*(Dcol-1)
		
  Scol = seq(
    focalList$centreCell[2]+cellDistFine[2]-Dcol+1,
    len= fineRelativeRes[2],by=1)[1:fineRelativeRes[2]]
		
  for(Drow in 1:(fineRelativeRes[1])){
   Dcell = firstCellInCol + Drow
   Srow = seq(
     focalList$centreCell[1]+cellDistFine[1]-Drow+1,
     len= fineRelativeRes[1],by=1)[1:fineRelativeRes[1]]
   for(Dbw in Sbw){
    kernelArray[Dcell,,Dbw] = as.vector(
      focalList$extended[[Dbw]][
        Srow, Scol
      ])
   } # for Dbw
  } #  for Drow
 } # for Dcol
	
 if(reorder){
  resEnv = new.env()
  with(resEnv, kernelArray <- kernelArray)
		
  with(resEnv,
    reorder <- reorderCellsTranslate(fineRelativeRes)
  )
		
  resFun = function(qx) {
   if(qx=='q1x') { #x>0, y=0
    resk = kernelArray
   } else if(qx=='q2x') { #|x|>|y| or |x|=|y|; x>0,y<0
    resk = kernelArray
   } else if(qx=='q2y') { #|y|>|x|; x>0 or x=0, y<0
    resk = kernelArray[reorder$q2,reorder$q2,,drop=FALSE]
   } else if(qx == 'q3y') { #|y|>|x|; x<0 or x=0, y<0
    resk = kernelArray[reorder$q3,reorder$q3,,drop=FALSE]
   } else if(qx == 'q3x') { #|x|>|y| or |x|=|y|; x<0 , y<0
    resk = aperm(kernelArray[reorder$q4,reorder$q4,,drop=FALSE],c(2,1,3))
   } else {
    warning('qx cant be ', qx)
    resk = NULL
   }
   resk
  }
		
  environment(resFun) = resEnv
  resHere = resFun
		
 } else {
  resHere = kernelArray
 }
 resHere
}
