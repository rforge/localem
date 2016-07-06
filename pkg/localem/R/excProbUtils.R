# Simulates case data with input risk threshold
# Computes the relative risk estimation of bootstrap data
simLemEst = function(bootId,
                     x,
                     lemObjects,
                     threshold,
                     bw,
                     tol,
                     maxIter
) {

  offsetRaster = trim(lemObjects$offset$offset)

  simEvents = geostatsp::simPoissonPP(threshold * offsetRaster)$events
  simEvents$id = over(simEvents, x)[["id"]]
  simEvents$count[!is.na(simEvents$id)] = 1

  simCounts = aggregate(simEvents$count,
                        by = list(id = as.character(simEvents$id)),
                        FUN = length
  )
  colnames(simCounts) = c("id","count")

  x$count = simCounts$count[match(x$id, simCounts$id)]
  x$count[is.na(x$count)] = 0

  result = riskEst(x = x,
                lemObjects = lemObjects,
                bw = bw,
                tol = tol,
                maxIter = maxIter)

  return(values(result))
}
