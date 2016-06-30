# Rasterizes the expected counts of the fine polygons over the fine grid cells
offsetSpToRaster = function(
  x,
  rasterFine,
  pattern="expected")
{

  #data showing in raster
  polyFineId = 1:length(x)
  rasterFineId = rasterFine[["idFine"]]

  nCellPerPoly = rep(1, length(x))
  nCellTable = table(values(rasterFineId))
  nCellPerPoly[as.numeric(names(nCellTable))] = nCellTable

  dataHere = x@data[,pattern] / nCellPerPoly
  dataHere[is.na(dataHere)] = 0

  forRasterHere = dataHere[values(rasterFineId)]
  forRasterHere[is.na(forRasterHere)] = 0

  #data not showing in raster
  notInRaster = which(!(polyFineId %in% values(rasterFineId)))

  if(length(notInRaster) > 0) {

    polyCentres = SpatialPoints(x[notInRaster,])
    polyCell = cellFromXY(rasterFineId, polyCentres@coords)

    dataNotInRaster = aggregate(dataHere[notInRaster],
                                list(cell = polyCell),
                                FUN = sum,
                                na.rm = TRUE)
    colnames(dataNotInRaster) = c("cell","expected")

    forRasterHere[dataNotInRaster$cell] = forRasterHere[dataNotInRaster$cell] + dataNotInRaster$expected
  }

  rasterFineArea = prod(res(rasterFine))

  res = raster(rasterFineId)
  values(res) = exp(log(forRasterHere) - log(rasterFineArea))

  return(res)
}
