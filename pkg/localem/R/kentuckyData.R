#' @title Observed cancer case data in Kentucky
#'
#' @description The dataset contains simulated lung cancer cases for the counties of Kentucky.
#'
#' @name kentuckyCounty
#'
#' @docType data
#'
#' @usage data(kentuckyCounty)
#'
#' @format A SpatialPolygonsDataFrame of Kentucky county boundaries and observed cancer cases.
#'
#' @details The observed lung cancer case data are simulated from a log Gaussian Cox process with cancer rates obtained from www.cancer-rates.info for deliberately unspecified years.
#'
#' @keywords datasets
#'
#' @examples
#' data(kentuckyCounty)
#'
#' countyBrks = classInt::classIntervals(kentuckyCounty$count, 
#' 		n = 7, style = 'quantile'
#' )
#' countyCol = classInt::findColours(countyBrks, 
#' 		pal = rev(RColorBrewer::brewer.pal(7, 'RdYlBu'))
#' )
#' 
#' mapmisc::map.new(kentuckyCounty)
#' plot(kentuckyCounty, 
#' 		col = countyCol, 
#' 		axes = FALSE, border = FALSE
#' )
#' 
#' data(kMap)
#' plot(kMap, add = TRUE)
#' 
#' legend('topleft', 
#' 		legend = rev(names(attr(countyCol, 'table'))), 
#' 		fill = rev(attr(countyCol, 'palette')), 
#' 		bg = 'white'
#' )
NULL

#' @title Population and expected cancer case data in Kentucky
#'
#' @description The dataset contains population data and expected lung cancer cases for the census tracts of Kentucky.
#'
#' @name kentuckyTract
#'
#' @docType data
#'
#' @usage data(kentuckyTract)
#'
#' @format A SpatialPolygonsDataFrame of Kentucky census tract boundaries, population and expected cancer cases.
#'
#' @details The expected lung cancer case data are computed with the age-sex population and cancer rates obtained from www.cancer-rates.info for deliberately unspecified years.
#'
#' @keywords datasets
#'
#' @examples
#' data(kentuckyTract)
#' kentuckyTract$expected[is.na(kentuckyTract$expected)] = 0
#'
#' ctBrks = classInt::classIntervals(kentuckyTract$expected, 
#'		n = 7, style = 'fixed', fixedBreaks = c(0, 1, 1.5, 2, 3, 4, 6, 9.2)
#' )
#' ctCol = classInt::findColours(ctBrks, 
#'		pal = rev(RColorBrewer::brewer.pal(7, 'RdYlBu'))
#' )
#'
#' mapmisc::map.new(kentuckyTract)
#' plot(kentuckyTract, 
#'		col = ctCol, 
#'		axes = FALSE, border = FALSE
#' )
#'
#' data(kMap)
#' plot(kMap, add = TRUE)
#'
#' legend('topleft', 
#' 		legend = rev(names(attr(ctCol, 'table'))), 
#'		fill = rev(attr(ctCol, 'palette')), 
#'		bg = 'white'
#' )
NULL

#' @title Map of Kentucky
#'
#' @description Stamen toner map of Kentucky
#'
#' @name kMap
#'
#' @docType data
#'
#' @usage data('kMap')
#'
#' @format A RayerLayer of Kentucky and surrounding areas. 
#'
#' @details The Stamen toner map is a detailed black and white background map to use with colorful overlays. 
#'
#' @keywords datasets
#'
#' @examples
#' data(kMap)
#' 
#' data(kentuckyCounty)
#'
#' countyBrks = classInt::classIntervals(kentuckyCounty$count, 
#' 		n = 7, style = 'quantile'
#' )
#' countyCol = classInt::findColours(countyBrks, 
#' 		pal = rev(RColorBrewer::brewer.pal(7, 'RdYlBu'))
#' )
#' 
#' mapmisc::map.new(kentuckyCounty)
#' plot(kentuckyCounty, 
#' 		col = countyCol, 
#' 		axes = FALSE, border = FALSE
#' )
#' plot(kMap, add = TRUE)
#' 
#' legend('topleft', 
#' 		legend = rev(names(attr(countyCol, 'table'))), 
#' 		fill = rev(attr(countyCol, 'palette')), 
#' 		bg = 'white'
#' )
#'
#' data(kentuckyTract)
#' kentuckyTract$expected[is.na(kentuckyTract$expected)] = 0
#'
#' ctBrks = classInt::classIntervals(kentuckyTract$expected, 
#'		n = 7, style = 'fixed', fixedBreaks = c(0, 1, 1.5, 2, 3, 4, 6, 9.2)
#' )
#' ctCol = classInt::findColours(ctBrks, 
#'		pal = rev(RColorBrewer::brewer.pal(7, 'RdYlBu'))
#' )
#'
#' mapmisc::map.new(kentuckyTract)
#' plot(kentuckyTract, 
#'		col = ctCol, 
#'		axes = FALSE, border = FALSE
#' )
#' plot(kMap, add = TRUE)
#'
#' legend('topleft', 
#' 		legend = rev(names(attr(ctCol, 'table'))), 
#'		fill = rev(attr(ctCol, 'palette')), 
#'		bg = 'white'
#' )
NULL
